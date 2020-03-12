glue.command(["multi-unset", "link-target", "variation", "cov_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_deletion_sequence", "cov_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_deletion_sequence", "sequence", "-a"]);
glue.command(["multi-delete", "cov_deletion", "-a"]);
glue.command(["multi-delete", "cov_deletion_sequence", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del%'"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del_detect%'"]);

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));

// create "detection" variation objects to detect deletions in each of the features respectively.
var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
_.each(featuresList, function(featureObj) {
	glue.inMode("reference/"+comparisonRefName+"/feature-location/"+featureObj.name, function() {
		var codonLabels= glue.getTableColumn(glue.command(["list", "labeled-codon"]), "codonLabel");
		var codonLabel2 = codonLabels[1];
		var codonLabelPenultimate = codonLabels[codonLabels.length - 2];
		glue.command(["create", "variation", "cov_aa_del_detect:"+featureObj.name, 
			"-t", "aminoAcidDeletion", "-c", codonLabel2, codonLabelPenultimate]);
	});
});

var deletionsSet = {};
var orf1aDeletions = {}; 
var orf1abDeletions = {}; 


/*
if non-codon-aligned deletions are detected at this step, we try to fix them if the deletion size is a multiple of 3.
This is because they are probably due to deletions which MAFFT did not properly codon-align.
Some examples:
nsp15:227:del in EPI_ISL_408670. 
nsp1:82-86:del in EPI_ISL_410044, EPI_ISL_413623
nsp2:268del in EPI_ISL_413564, EPI_ISL_413568, EPI_ISL_413573, EPI_ISL_413577, 
               EPI_ISL_413580, EPI_ISL_413581, EPI_ISL_413582, EPI_ISL_413583
S:145del in EPI_ISL_413522
the fix involves 
(a) removing the segments of the sequence from AL_GISAID_UNCONSTRAINED 
(b) adding the sequence to AL_GISAID_CONSTRAINED 
(c) recomputing this row of AL_GISAID_CONSTRAINED using the codon-aware compound aligner.
(d) deriving the segments back to AL_GISAID_UNCONSTRAINED.
(e) attempting to generate deletions again.
*/

var realignedMemberObjs = [];

glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
	var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", "sequence.analyse_aa_deletions = true"]));
	deletionScan(almtMemberObjs, true);
});

_.each(realignedMemberObjs, function(realignedMemberObj) {
	var sourceName = realignedMemberObj["sequence.source.name"];
	var sequenceID = realignedMemberObj["sequence.sequenceID"];
	
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED/member/"+sourceName+"/"+sequenceID, function() {
		glue.command(["remove", "segment", "-a"]);
	});

	glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
		glue.command(["add", "member", sourceName, sequenceID]);
	});

	glue.command(["compute", "alignment", "AL_GISAID_CONSTRAINED", "covCompoundAligner", 
		"-w", "sequence.source.name = '"+sourceName+"' and sequence.sequenceID = '"+sequenceID+"'"]);

	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
		glue.command(["derive", "segments", "AL_GISAID_CONSTRAINED", 
			"-w", "sequence.source.name = '"+sourceName+"' and sequence.sequenceID = '"+sequenceID+"'"]);
	});
});

glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
	deletionScan(realignedMemberObjs, false);
});


function deletionScan(almtMemberObjs, attemptRealign) {
	_.each(almtMemberObjs, function(almtMemberObj) {
		glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
			
			var allFeatureDelObjs = [];
			
			_.each(featuresList, function(featureObj) {
				var memberDelObjs = glue.tableToObjects(glue.command(["variation", "scan", 
					"-r", comparisonRefName, "-f", featureObj.name, 
					"--whereClause", "name = 'cov_aa_del_detect:"+featureObj.name+"'", 
					"--excludeAbsent", "--showMatchesAsTable"]));
				
				_.each(memberDelObjs, function(memberDelObj) {
					memberDelObj.featureName = featureObj.name;
					memberDelObj.parentFeatureName = featureObj["parent.name"];
				});
				
				allFeatureDelObjs = allFeatureDelObjs.concat(memberDelObjs);
			});
			var realign = false;
			_.each(allFeatureDelObjs, function(memberDelObj) {
				var codonAligned = memberDelObj.deletionIsCodonAligned;
				if(!codonAligned) {
					if(attemptRealign && memberDelObj.deletedRefNts.length % 3 == 0) { // deletion is multiple of 3, so may attempt realignment.
						realign = true;
					} else {
						glue.log("INFO", "alignment member", almtMemberObj);
						glue.log("INFO", "non-codon-aligned deletion", memberDelObj);
						throw new Error("Non-codon-aligned deletion, either length is not a mulitple of 3 or realignment failed! "+
								"We could add this but what are the implications for protein translation elsewhere in the feature? "+
								"Investigate, then consider setting analyse_aa_deletions false in covLoadSequenceData.glue");
					}
					
				}
			});
			if(realign) {
				realignedMemberObjs.push(almtMemberObj);
				return;
			}
			
			_.each(allFeatureDelObjs, function(memberDelObj) {
				var codonAligned = memberDelObj.deletionIsCodonAligned;
				var deletionID = memberDelObj.featureName+":ca:"+memberDelObj.refFirstCodonDeleted+":"+memberDelObj.refLastCodonDeleted;
				var deletionObj = deletionsSet[deletionID];
				if(deletionObj == null) {
					deletionObj = {
						id: deletionID,
						feature: memberDelObj.featureName,
						parentFeature: memberDelObj.parentFeatureName,
						codonAligned: memberDelObj.deletionIsCodonAligned,
						startCodon: memberDelObj.refFirstCodonDeleted,
						endCodon: memberDelObj.refLastCodonDeleted,
						refNtStart: memberDelObj.refFirstNtDeleted,
						refNtEnd: memberDelObj.refLastNtDeleted,
						memberSeqs: []
					};
					deletionsSet[deletionID] = deletionObj;
					if(memberDelObj.featureName == "ORF_1a") {
						orf1aDeletions[memberDelObj.refFirstNtDeleted+":"+memberDelObj.refLastNtDeleted] = deletionObj; 
					}
					if(memberDelObj.featureName == "ORF_1ab") {
						orf1abDeletions[memberDelObj.refFirstNtDeleted+":"+memberDelObj.refLastNtDeleted] = deletionObj; 
					}

				}
				deletionObj.memberSeqs.push(almtMemberObj);
			});
			
		});
	}); 
}

_.each(_.values(deletionsSet), function(deletionObj) {
	if(deletionObj.feature == "ORF_1a" || deletionObj.feature == "ORF_1ab") {
		return;
	}
	createDeletion(deletionObj);
});

// create any ORF1a/ORF1ab deletions which are not already represented by NSP deletions
// eg if they span cleavage locations.
_.each(_.values(orf1aDeletions), function(deletionObj) {
	if(deletionObj.skipCreation) {
		return;
	}
	createDeletion(deletionObj);
});


_.each(_.values(orf1abDeletions), function(deletionObj) {
	if(deletionObj.skipCreation) {
		return;
	}
	createDeletion(deletionObj);
});

function createDeletion(deletionObj) {
	glue.log("FINEST", "Creating deletion object", deletionObj);
	var variationName = "cov_aa_del:"+deletionObj.id;
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+deletionObj.feature, function() {
		glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidDeletion", 
			"--labeledCodon", deletionObj.startCodon, deletionObj.endCodon]);
	});
	
	glue.command(["create", "custom-table-row", "cov_deletion", deletionObj.id]);
	glue.inMode("custom-table-row/cov_deletion/"+deletionObj.id, function() {
		var displayName;
		if(deletionObj.startCodon == deletionObj.endCodon) {
			displayName = deletionObj.startCodon+":del";	
		} else {
			displayName = deletionObj.startCodon+"-"+deletionObj.endCodon+":del";	
		}
		glue.command(["set", "field", "display_name", displayName]);
		glue.command(["set", "field", "start_codon", deletionObj.startCodon]);		
		glue.command(["set", "field", "start_codon_int", parseInt(deletionObj.startCodon)]);		
		glue.command(["set", "field", "end_codon", deletionObj.endCodon]);		
		glue.command(["set", "field", "end_codon_int", parseInt(deletionObj.endCodon)]);		
		glue.command(["set", "field", "codon_aligned", deletionObj.codonAligned]);		
		glue.command(["set", "field", "reference_nt_start", deletionObj.refNtStart]);		
		glue.command(["set", "field", "reference_nt_end", deletionObj.refNtEnd]);		
		glue.command(["set", "field", "num_seqs", deletionObj.memberSeqs.length]);
		if(deletionObj.parentFeature == "ORF_1a") {
			glue.command(["set", "field", "parent_feature", "ORF_1a"]);
			var parent1aDelObj = orf1aDeletions[deletionObj.refNtStart+":"+deletionObj.refNtEnd];
			parent1aDelObj.skipCreation = true;
			var parent1abDelObj = orf1abDeletions[deletionObj.refNtStart+":"+deletionObj.refNtEnd];
			parent1abDelObj.skipCreation = true;
			glue.command(["set", "field", "parent_start_codon", parent1aDelObj.startCodon]);
			glue.command(["set", "field", "parent_end_codon", parent1aDelObj.endCodon]);
		}
		if(deletionObj.parentFeature == "ORF_1ab") {
			glue.command(["set", "field", "parent_feature", "ORF_1ab"]);
			var parentDelObj = orf1abDeletions[deletionObj.refNtStart+":"+deletionObj.refNtEnd];
			parentDelObj.skipCreation = true;
			glue.command(["set", "field", "parent_start_codon", parentDelObj.startCodon]);
			glue.command(["set", "field", "parent_end_codon", parentDelObj.endCodon]);
		}
		glue.command(["set", "link-target", "variation", 
			"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+deletionObj.feature+
			"/variation/"+variationName]);
	});
	
	_.each(deletionObj.memberSeqs, function(memberObj) {
		var sourceName = memberObj["sequence.source.name"];
		var sequenceID = memberObj["sequence.sequenceID"];
		var linkObjId = deletionObj.id+":"+sequenceID;
		glue.command(["create", "custom-table-row", "cov_deletion_sequence", linkObjId]);
		glue.inMode("custom-table-row/cov_deletion_sequence/"+linkObjId, function() {
			glue.command(["set", "link-target", "cov_deletion", "custom-table-row/cov_deletion/"+deletionObj.id]);
			glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
		});
	});
}