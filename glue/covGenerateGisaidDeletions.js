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


glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
	var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", "sequence.analyse_aa_deletions = true"]));
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
});



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