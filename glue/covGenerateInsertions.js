
var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));

glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_ins_detect%'"]);

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

var ntInsertionsSet = {};
var orf1aNtInsertions = {}; 
var orf1abNtInsertions = {}; 

var insertionsSet = {};
var orf1aInsertions = {}; 
var orf1abInsertions = {}; 

var processed = 0;

// variant cache test 
// var whereClause = "sequence.sequenceID in ('EPI_ISL_500981', 'EPI_ISL_465549')";

//just delete everything
//var whereClause = "false";

//sequences containing (a) codon-aligned insertion in NSP16, (b) no insertions
//var whereClause = "sequence.analyse_variation = true and sequence.sequenceID in ('EPI_ISL_414588', 'EPI_ISL_402125')"

// production
var whereClause = "sequence.analyse_variation = true and cg_insertions_from_cache = false";

glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
	var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", whereClause]));
	_.each(almtMemberObjs, function(almtMemberObj) {
		glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
			
			var allFeatureInsObjs = [];

			_.each(featuresList, function(featureObj) {
				var memberInsObjs = glue.tableToObjects(glue.command(["variation", "scan", 
					"-r", comparisonRefName, "-f", featureObj.name, 
					"--whereClause", "name = 'cov_aa_ins_detect:"+featureObj.name+"'", 
					"--excludeAbsent", "--showMatchesAsTable"]));
				
				_.each(memberInsObjs, function(memberInsObj) {
					memberInsObj.featureName = featureObj.name;
					memberInsObj.parentFeatureName = featureObj["parent.name"];
				});
				
				allFeatureInsObjs = allFeatureInsObjs.concat(memberInsObjs);
			});
		
			_.each(allFeatureInsObjs, function(memberInsObj) {
				// there seem to be spurious very long insertions 
				if(memberInsObj.insertedQryNts.length > 100) {
					return;
				}
				// skip insertions consisting purely of Ns
				if(memberInsObj.insertedQryNts.match("^N+$") != null) {
					return;
				}
				var ntHash = stringHash(memberInsObj.insertedQryNts);
				var ntInsertionID = memberInsObj.featureName+":nca:"+memberInsObj.refLastNtBeforeIns+":"+ntHash+":"+memberInsObj.refFirstNtAfterIns;
				var ntInsertionObj = ntInsertionsSet[ntInsertionID];
				if(ntInsertionObj == null) {
					ntInsertionObj = {
						id: ntInsertionID,
						feature: memberInsObj.featureName,
						parentFeature: memberInsObj.parentFeatureName,
						lastRefNtBefore: memberInsObj.refLastNtBeforeIns,
						firstRefNtAfter: memberInsObj.refFirstNtAfterIns,
						insertedNts: memberInsObj.insertedQryNts,
						memberSeqs: []
					};
					ntInsertionsSet[ntInsertionID] = ntInsertionObj;
					if(memberInsObj.featureName == "ORF_1a") {
						orf1aNtInsertions[memberInsObj.refLastNtBeforeIns+":"+ntHash+":"+memberInsObj.refFirstNtAfterIns] = ntInsertionObj; 
					}
					if(memberInsObj.featureName == "ORF_1ab") {
						orf1abNtInsertions[memberInsObj.refLastNtBeforeIns+":"+ntHash+":"+memberInsObj.refFirstNtAfterIns] = ntInsertionObj; 
					}
	
				}
				ntInsertionObj.memberSeqs.push(almtMemberObj);
				
				if(memberInsObj.insertionIsCodonAligned) {
					var aaHash = stringHash(memberInsObj.insertedQryAas);

					var insertionID =  memberInsObj.featureName+":ca:"+memberInsObj.refLastCodonBeforeIns+":"+aaHash+":"+memberInsObj.refFirstCodonAfterIns;
					var insertionObj = insertionsSet[insertionID];
					if(insertionObj == null) {
						insertionObj = {
							id: insertionID,
							feature: memberInsObj.featureName,
							parentFeature:  memberInsObj.parentFeatureName,
							lastCodonBefore: memberInsObj.refLastCodonBeforeIns,
							firstCodonAfter: memberInsObj.refFirstCodonAfterIns,
							lastRefNtBefore: memberInsObj.refLastNtBeforeIns,
							firstRefNtAfter: memberInsObj.refFirstNtAfterIns,
							insertedAas: memberInsObj.insertedQryAas,
							ntInsertionID: ntInsertionID,
							memberSeqs: []
						};
						insertionsSet[insertionID] = insertionObj;
						if(memberInsObj.featureName == "ORF_1a") {
							orf1aInsertions[memberInsObj.refLastNtBeforeIns+":"+aaHash+":"+memberInsObj.refFirstNtAfterIns] = insertionObj; 
						}
						if(memberInsObj.featureName == "ORF_1ab") {
							orf1abInsertions[memberInsObj.refLastNtBeforeIns+":"+aaHash+":"+memberInsObj.refFirstNtAfterIns] = insertionObj; 
						}
	
					}
					insertionObj.memberSeqs.push(almtMemberObj);
				}
			});
			processed++;
			if(processed % 250 == 0) {
				glue.log("FINEST", "Processed "+processed+" alignment members for insertions");
				glue.command(["new-context"]);
			}
		});
	}); 
});

//create NT (non codon aligned) insertions

_.each(_.values(ntInsertionsSet), function(ntInsertionObj) {
	if(ntInsertionObj.feature == "ORF_1a" || ntInsertionObj.feature == "ORF_1ab") {
		return;
	}
	createNtInsertion(ntInsertionObj);
});
// create any ORF1a/ORF1ab ntInsertions which are not already represented by NSP ntInsertions
// eg if they span cleavage locations.
_.each(_.values(orf1aNtInsertions), function(ntInsertionObj) {
	if(ntInsertionObj.skipCreation) {
		return;
	}
	createNtInsertion(ntInsertionObj);
});
_.each(_.values(orf1abNtInsertions), function(ntInsertionObj) {
	if(ntInsertionObj.skipCreation) {
		return;
	}
	createNtInsertion(ntInsertionObj);
});

// ---------

//create codon aligned insertions
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


function createNtInsertion(ntInsertionObj) {
	glue.log("FINEST", "Creating insertion object", ntInsertionObj);
	var ntHash = stringHash(ntInsertionObj.insertedNts);

	var variationName = "cov_nt_ins:"+ntInsertionObj.id;
	var variationExists = false;

	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+ntInsertionObj.feature, function() {
		var existing = glue.tableToObjects(glue.command(["list", "variation", "-w", "name = '"+variationName+"'"]));
		if(existing.length > 0) {
			variationExists = true;
		}
		if(!variationExists) {
			glue.command(["create", "variation", variationName, 
				"-t", "nucleotideInsertion", 
				"--nucleotide", ntInsertionObj.lastRefNtBefore, ntInsertionObj.firstRefNtAfter]);
			glue.inMode("variation/"+variationName, function() {
				glue.command(["set", "metatag", "MIN_INSERTION_LENGTH_NTS", ntInsertionObj.insertedNts.length]);
				glue.command(["set", "metatag", "MAX_INSERTION_LENGTH_NTS", ntInsertionObj.insertedNts.length]);
			});
		}
	});

	if(ntInsertionObj.parentFeature == "ORF_1a") {
		var parent1aInsObj = orf1aNtInsertions[ntInsertionObj.lastRefNtBefore+":"+ntHash+":"+ntInsertionObj.firstRefNtAfter];
		if(parent1aInsObj != null) {
			parent1aInsObj.skipCreation = true;
		}
		var parent1abInsObj = orf1abNtInsertions[ntInsertionObj.lastRefNtBefore+":"+ntHash+":"+ntInsertionObj.firstRefNtAfter];
		if(parent1abInsObj != null) {
			parent1abInsObj.skipCreation = true;
		}
	}
	if(ntInsertionObj.parentFeature == "ORF_1ab") {
		var parentInsObj = orf1abNtInsertions[ntInsertionObj.lastRefNtBefore+":"+ntHash+":"+ntInsertionObj.firstRefNtAfter];
		if(parentInsObj != null) {
			parentInsObj.skipCreation = true;
		}
	}

	
	if(!variationExists) {
		glue.command(["create", "custom-table-row", "cov_nt_insertion", ntInsertionObj.id]);
		glue.inMode("custom-table-row/cov_nt_insertion/"+ntInsertionObj.id, function() {
			var displayName;
			displayName = ntInsertionObj.lastRefNtBefore+"-"+ntInsertionObj.insertedNts+"-"+ntInsertionObj.firstRefNtAfter;	
			glue.command(["set", "field", "display_name", displayName]);
			glue.command(["set", "field", "last_ref_nt_before", ntInsertionObj.lastRefNtBefore]);		
			glue.command(["set", "field", "first_ref_nt_after", ntInsertionObj.firstRefNtAfter]);		
			glue.command(["set", "field", "inserted_nts_length", ntInsertionObj.insertedNts.length]);
			glue.command(["set", "field", "inserted_nts", ntInsertionObj.insertedNts]);
			
			if(ntInsertionObj.parentFeature == "ORF_1a") {
				glue.command(["set", "field", "parent_feature", "ORF_1a"]);
			}
			if(ntInsertionObj.parentFeature == "ORF_1ab") {
				glue.command(["set", "field", "parent_feature", "ORF_1ab"]);
			}
			glue.command(["set", "link-target", "variation", 
				"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+ntInsertionObj.feature+
				"/variation/"+variationName]);
		});
	}
	
	_.each(ntInsertionObj.memberSeqs, function(memberObj) {
		var sourceName = memberObj["sequence.source.name"];
		var sequenceID = memberObj["sequence.sequenceID"];
		var linkObjId = ntInsertionObj.id+":"+sourceName+":"+sequenceID;
		glue.command(["create", "custom-table-row", "cov_nt_insertion_sequence", linkObjId]);
		glue.inMode("custom-table-row/cov_nt_insertion_sequence/"+linkObjId, function() {
			glue.command(["set", "link-target", "cov_nt_insertion", "custom-table-row/cov_nt_insertion/"+ntInsertionObj.id]);
			glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
		});
	});
}


function createInsertion(insertionObj) {
	glue.log("FINEST", "Creating insertion object", insertionObj);
	var aaHash = stringHash(insertionObj.insertedAas);

	var variationName = "cov_aa_ins:"+insertionObj.id;
	var variationExists = false;

	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+insertionObj.feature, function() {
		var existing = glue.tableToObjects(glue.command(["list", "variation", "-w", "name = '"+variationName+"'"]));
		if(existing.length > 0) {
			variationExists = true;
		}
		if(!variationExists) {
			glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidInsertion", 
			"--labeledCodon", insertionObj.lastCodonBefore, insertionObj.firstCodonAfter]);
			glue.inMode("variation/"+variationName, function() {
				glue.command(["set", "metatag", "MIN_INSERTION_LENGTH_AAS", insertionObj.insertedAas.length]);
				glue.command(["set", "metatag", "MAX_INSERTION_LENGTH_AAS", insertionObj.insertedAas.length]);
			});
		}
	});
	
	if(insertionObj.parentFeature == "ORF_1a") {
		var parent1aInsObj = orf1aInsertions[insertionObj.lastRefNtBefore+":"+aaHash+":"+insertionObj.firstRefNtAfter];
		if(parent1aInsObj != null) {
			parent1aInsObj.skipCreation = true;
		}
		var parent1abInsObj = orf1abInsertions[insertionObj.lastRefNtBefore+":"+aaHash+":"+insertionObj.firstRefNtAfter];
		if(parent1abInsObj != null) {
			parent1abInsObj.skipCreation = true;
		}
	}
	if(insertionObj.parentFeature == "ORF_1ab") {
		var parentInsObj = orf1abInsertions[insertionObj.lastRefNtBefore+":"+aaHash+":"+insertionObj.firstRefNtAfter];
		if(parentInsObj != null) {
			parentInsObj.skipCreation = true;
		}
	}
	
	if(!variationExists) {
		glue.command(["create", "custom-table-row", "cov_insertion", insertionObj.id]);
		glue.inMode("custom-table-row/cov_insertion/"+insertionObj.id, function() {
			var displayName;
			displayName = insertionObj.lastCodonBefore+"-"+insertionObj.insertedAas+"-"+insertionObj.firstCodonAfter;	
			glue.command(["set", "field", "display_name", displayName]);
			glue.command(["set", "field", "last_codon_before", insertionObj.lastCodonBefore]);		
			glue.command(["set", "field", "last_codon_before_int", parseInt(insertionObj.lastCodonBefore)]);		
			glue.command(["set", "field", "first_codon_after", insertionObj.firstCodonAfter]);		
			glue.command(["set", "field", "first_codon_after_int", parseInt(insertionObj.firstCodonAfter)]);		
			glue.command(["set", "field", "last_ref_nt_before", insertionObj.lastRefNtBefore]);		
			glue.command(["set", "field", "first_ref_nt_after", insertionObj.firstRefNtAfter]);		
			glue.command(["set", "field", "inserted_aas_length", insertionObj.insertedAas.length]);
			glue.command(["set", "field", "inserted_aas", insertionObj.insertedAas]);
			
			if(insertionObj.parentFeature == "ORF_1a") {
				glue.command(["set", "field", "parent_feature", "ORF_1a"]);
				var parent1aInsObj = orf1aInsertions[insertionObj.lastRefNtBefore+":"+aaHash+":"+insertionObj.firstRefNtAfter];
				if(parent1aInsObj != null) {
					glue.command(["set", "field", "parent_last_codon_before", parent1aInsObj.lastCodonBefore]);
					glue.command(["set", "field", "parent_first_codon_after", parent1aInsObj.firstCodonAfter]);
				}
			}
			if(insertionObj.parentFeature == "ORF_1ab") {
				glue.command(["set", "field", "parent_feature", "ORF_1ab"]);
				var parentInsObj = orf1abInsertions[insertionObj.lastRefNtBefore+":"+aaHash+":"+insertionObj.firstRefNtAfter];
				if(parentInsObj != null) {
					glue.command(["set", "field", "parent_last_codon_before", parentInsObj.lastCodonBefore]);
					glue.command(["set", "field", "parent_first_codon_after", parentInsObj.firstCodonAfter]);
				}
			}
			glue.command(["set", "link-target", "variation", 
				"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+insertionObj.feature+
				"/variation/"+variationName]);
			glue.command(["set", "link-target", "cov_nt_insertion", 
				"custom-table-row/cov_nt_insertion/"+insertionObj.ntInsertionID]);
	
		});
	}	
	_.each(insertionObj.memberSeqs, function(memberObj) {
		var sourceName = memberObj["sequence.source.name"];
		var sequenceID = memberObj["sequence.sequenceID"];
		var linkObjId = insertionObj.id+":"+sourceName+":"+sequenceID;
		glue.command(["create", "custom-table-row", "cov_insertion_sequence", linkObjId]);
		glue.inMode("custom-table-row/cov_insertion_sequence/"+linkObjId, function() {
			glue.command(["set", "link-target", "cov_insertion", "custom-table-row/cov_insertion/"+insertionObj.id]);
			glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
		});
	});
}







function stringHash(string) {
	var hash = 0, i, chr;
    for (i = 0; i < string.length; i++) {
      chr   = string.charCodeAt(i);
      hash  = ((hash << 5) - hash) + chr;
      hash |= 0; // Convert to 32bit integer
    }
    return hash;
}
