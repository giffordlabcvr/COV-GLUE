glue.command(["multi-unset", "link-target", "variation", "cov_insertion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_insertion_sequence", "cov_insertion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_insertion_sequence", "sequence", "-a"]);
glue.command(["multi-delete", "cov_insertion", "-a"]);
glue.command(["multi-delete", "cov_insertion_sequence", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_ins%'"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_ins_detect%'"]);

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));

// create "detection" variation objects to detect insertions in each of the features respectively.
var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
_.each(featuresList, function(featureObj) {
	glue.inMode("reference/"+comparisonRefName+"/feature-location/"+featureObj.name, function() {
		var codonLabels= glue.getTableColumn(glue.command(["list", "labeled-codon"]), "codonLabel");
		var codonLabel1 = codonLabels[0];
		var codonLabelLast = codonLabels[codonLabels.length - 1];
		glue.command(["create", "variation", "cov_aa_ins_detect:"+featureObj.name, 
			"-t", "aminoAcidInsertion", "-c", codonLabel1, codonLabelLast]);
	});
});

var insertionsSet = {};
var orf1aInsertions = {}; 
var orf1abInsertions = {}; 
var numInsertions = 0;

_.each(featuresList, function(featureObj) {
	glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
		var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", "sequence.analyse_aa_insertions = true"]));
		_.each(almtMemberObjs, function(almtMemberObj) {
			glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
				
				var memberInsObjs = glue.tableToObjects(glue.command(["variation", "scan", 
					"-r", comparisonRefName, "-f", featureObj.name, 
					"--whereClause", "name = 'cov_aa_ins_detect:"+featureObj.name+"'", 
					"--excludeAbsent", "--showMatchesAsTable"]));
				
				_.each(memberInsObjs, function(memberInsObj) {
					var codonAligned = memberInsObj.insertionIsCodonAligned;
					var insertionID;
					if(codonAligned) {
						insertionID = featureObj.name+":ca:"+memberInsObj.refLastCodonBeforeIns+":"+memberInsObj.insertedQryAas.length+":"+memberInsObj.refFirstCodonAfterIns;
					} else {
						glue.log("INFO", "alignment member", almtMemberObj);
						glue.log("INFO", "non-codon-aligned insertion", memberInsObj);
						throw new Error("Non-codon-aligned insertion! We could add this but what are the implications for protein translation elsewhere in the feature");
					}
					var insertionObj = insertionsSet[insertionID];
					if(insertionObj == null) {
						insertionObj = {
							id: insertionID,
							feature: featureObj.name,
							parentFeature: featureObj["parent.name"],
							codonAligned: memberInsObj.insertionIsCodonAligned,
							lastCodonBefore: memberInsObj.refLastCodonBeforeIns,
							firstCodonAfter: memberInsObj.refFirstCodonAfterIns,
							lastRefNtBefore: memberInsObj.refLastNtBeforeIns,
							firstRefNtAfter: memberInsObj.refFirstNtAfterIns,
							insertedAas: memberInsObj.insertedQryAas,
							memberSeqs: []
						};
						insertionsSet[insertionID] = insertionObj;
						if(featureObj.name == "ORF_1a") {
							orf1aInsertions[memberInsObj.refLastNtBeforeIns+":"+memberInsObj.insertedQryAas.length+":"+memberInsObj.refFirstNtAfterIns] = insertionObj; 
						}
						if(featureObj.name == "ORF_1ab") {
							orf1abInsertions[memberInsObj.refLastNtBeforeIns+":"+memberInsObj.insertedQryAas.length+":"+memberInsObj.refFirstNtAfterIns] = insertionObj; 
						}

					}
					insertionObj.memberSeqs.push(almtMemberObj);
				});
			});
		}); 
	});
});

_.each(_.values(insertionsSet), function(insertionObj) {
	if(insertionObj.feature == "ORF_1a" || insertionObj.feature == "ORF_1ab") {
		return;
	}
	createInsertion(insertionObj);
});

// create any ORF1a/ORF1ab insertions which are not already represented by NSP insertions
// eg if they span cleavage locations.
_.each(_.values(orf1aInsertions), function(insertionObj) {
	if(insertionObj.skipCreation) {
		return;
	}
	createInsertion(insertionObj);
});


_.each(_.values(orf1abInsertions), function(insertionObj) {
	if(insertionObj.skipCreation) {
		return;
	}
	createInsertion(insertionObj);
});

function createInsertion(insertionObj) {
	numInsertions++;
	glue.log("FINEST", "Creating insertion object", insertionObj);
	var variationName = "cov_aa_ins:"+insertionObj.id;
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+insertionObj.feature, function() {
		glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidInsertion", 
			"--labeledCodon", insertionObj.lastCodonBefore, insertionObj.firstCodonAfter]);
		glue.inMode("variation/"+variationName, function() {
			glue.command(["set", "metatag", "MIN_INSERTION_LENGTH_AAS", insertionObj.insertedAas.length]);
			glue.command(["set", "metatag", "MAX_INSERTION_LENGTH_AAS", insertionObj.insertedAas.length]);
		});
	});
	
	glue.command(["create", "custom-table-row", "cov_insertion", insertionObj.id]);
	glue.inMode("custom-table-row/cov_insertion/"+insertionObj.id, function() {
		var displayName;
		displayName = insertionObj.lastCodonBefore+"-"+insertionObj.insertedAas.length+"-"+insertionObj.firstCodonAfter+":ins";	
		glue.command(["set", "field", "display_name", displayName]);
		glue.command(["set", "field", "last_codon_before", insertionObj.lastCodonBefore]);		
		glue.command(["set", "field", "last_codon_before_int", parseInt(insertionObj.lastCodonBefore)]);		
		glue.command(["set", "field", "first_codon_after", insertionObj.firstCodonAfter]);		
		glue.command(["set", "field", "first_codon_after_int", parseInt(insertionObj.firstCodonAfter)]);		
		glue.command(["set", "field", "codon_aligned", insertionObj.codonAligned]);		
		glue.command(["set", "field", "last_ref_nt_before", insertionObj.lastRefNtBefore]);		
		glue.command(["set", "field", "first_ref_nt_after", insertionObj.firstRefNtAfter]);		
		glue.command(["set", "field", "num_seqs", insertionObj.memberSeqs.length]);
		if(insertionObj.parentFeature == "ORF_1a") {
			glue.command(["set", "field", "parent_feature", "ORF_1a"]);
			var parent1aInsObj = orf1aInsertions[insertionObj.lastRefNtBefore+":"+insertionObj.insertedAas.length+":"+insertionObj.firstRefNtAfter];
			parent1aInsObj.skipCreation = true;
			var parent1abInsObj = orf1abInsertions[insertionObj.lastRefNtBefore+":"+insertionObj.insertedAas.length+":"+insertionObj.firstRefNtAfter];
			parent1abInsObj.skipCreation = true;
			glue.command(["set", "field", "parent_last_codon_before", parent1aInsObj.lastCodonBefore]);
			glue.command(["set", "field", "parent_first_codon_after", parent1aInsObj.firstCodonAfter]);
		}
		if(insertionObj.parentFeature == "ORF_1ab") {
			glue.command(["set", "field", "parent_feature", "ORF_1ab"]);
			var parentInsObj = orf1abInsertions[insertionObj.lastRefNtBefore+":"+insertionObj.insertedAas.length+":"+insertionObj.firstRefNtAfter];
			parentInsObj.skipCreation = true;
			glue.command(["set", "field", "parent_last_codon_before", parentInsObj.lastCodonBefore]);
			glue.command(["set", "field", "parent_first_codon_after", parentInsObj.firstCodonAfter]);
		}
		glue.command(["set", "link-target", "variation", 
			"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+insertionObj.feature+
			"/variation/"+variationName]);
	});
	
	_.each(insertionObj.memberSeqs, function(memberObj) {
		var sourceName = memberObj["sequence.source.name"];
		var sequenceID = memberObj["sequence.sequenceID"];
		var linkObjId = insertionObj.id+":"+sequenceID;
		glue.command(["create", "custom-table-row", "cov_insertion_sequence", linkObjId]);
		glue.inMode("custom-table-row/cov_insertion_sequence/"+linkObjId, function() {
			glue.command(["set", "link-target", "cov_insertion", "custom-table-row/cov_insertion/"+insertionObj.id]);
			glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
		});
	});
}


if(numInsertions > 1) {
	throw new Error("Expected single NSP6 insertion in a Swiss sequence, please check.");
}