glue.command(["multi-unset", "link-target", "variation", "cov_replacement", "-a"]);
glue.command(["multi-unset", "link-target", "cov_replacement_sequence", "cov_replacement", "-a"]);
glue.command(["multi-unset", "link-target", "cov_replacement_sequence", "sequence", "-a"]);
glue.command(["multi-delete", "cov_replacement", "-a"]);
glue.command(["multi-delete", "cov_replacement_sequence", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_rpl%'"]);

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));


var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
var replacementsSet = {};
var orf1aReplacements = {}; 
var orf1abReplacements = {}; 

_.each(featuresList, function(featureObj) {
	var refAaObjsMap = {};
	glue.inMode("reference/"+comparisonRefName+"/feature-location/"+featureObj.name, function() {
		var refAaObjs = glue.tableToObjects(glue.command(["amino-acid"]));
		_.each(refAaObjs, function(refAaObj) {
			refAaObjsMap[refAaObj.codonLabel] = refAaObj;
		});
	});
	glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
		var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", "sequence.analyse_aa_replacements = true"]));
		_.each(almtMemberObjs, function(almtMemberObj) {
			glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
				var memberAaObjs = glue.tableToObjects(glue.command(["amino-acid", "-r", comparisonRefName, "-f", featureObj.name]));
				_.each(memberAaObjs, function(memberAaObj) {
					// Require no Ns in the codonNts in order to generate a replacement,
					// unless the replacement is unambiguously a single AA residue.
					// This means we are interpreting N as 'unable to sequence' rather than
					// 'equal proportion A, C, G, T' 
					if(memberAaObj.definiteAas != null && memberAaObj.definiteAas != "" &&
							(memberAaObj.definiteAas.length == 1 || memberAaObj.codonNts.indexOf('N') < 0)) {
						refAaObj = refAaObjsMap[memberAaObj.codonLabel];
						if(refAaObj != null && refAaObj.definiteAas != null && refAaObj.definiteAas != "" && 
								refAaObj.definiteAas != memberAaObj.definiteAas) {
							
							var refAas = refAaObj.definiteAas.split('');
							var memberAas = memberAaObj.definiteAas.split('');
							_.each(refAas, function(refAa) {
								_.each(memberAas, function(memberAa) {
									if(refAa != memberAa) {
										var replacementID = featureObj.name+":"+refAa+":"+memberAaObj.codonLabel+":"+memberAa;
										var replacementObj = replacementsSet[replacementID];
										if(replacementObj == null) {
											replacementObj = {
												id: replacementID,
												feature: featureObj.name,
												parentFeature: featureObj["parent.name"],
												codonLabel: memberAaObj.codonLabel,
												refNt: memberAaObj.relRefNt,
												refAa: refAa,
												replacementAa: memberAa,
												memberSeqs: []
											};
											replacementsSet[replacementID] = replacementObj;
											if(featureObj.name == "ORF_1a") {
												orf1aReplacements[memberAaObj.relRefNt+":"+memberAa] = replacementObj; 
											}
											if(featureObj.name == "ORF_1ab") {
												orf1abReplacements[memberAaObj.relRefNt+":"+memberAa] = replacementObj; 
											}
										}
										replacementObj.memberSeqs.push(almtMemberObj);
									}
								});
							});
						}
					}
				});
			});
			
		}); 
	});
	
});

_.each(_.values(replacementsSet), function(replacementObj) {
	if(replacementObj.feature == "ORF_1a" || replacementObj.feature == "ORF_1ab") {
		return;
	}
	glue.log("FINEST", "Creating replacement object", replacementObj);
	var variationName = "cov_aa_rpl:"+replacementObj.id;
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+replacementObj.feature, function() {
		glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidSimplePolymorphism", 
			"--labeledCodon", replacementObj.codonLabel, replacementObj.codonLabel]);
		glue.inMode("variation/"+variationName, function() {
			glue.command(["set", "metatag", "SIMPLE_AA_PATTERN", replacementObj.replacementAa]);
			glue.command(["set", "metatag", "MIN_COMBINED_TRIPLET_FRACTION", 0.25]);
		});
	});

	var hanada_radical_I;
	var hanada_radical_II;
	var hanada_radical_III;
	var grantham_distance_double;
	var grantham_distance_int;
	var miyata_distance;
	var classifyReplacement = false;
	
	if(replacementObj.refAa != '*' && replacementObj.refAa != 'X'
		 && replacementObj.replacementAa != '*' && replacementObj.replacementAa != 'X') {
		classifyReplacement = true;
		glue.inMode("module/covHanada2006ReplacementClassifier", function() {
			var classifierResults = glue.tableToObjects(glue.command(["classify", "replacement", replacementObj.refAa, replacementObj.replacementAa]));
			hanada_radical_I = classifierResults[0].radical;
			hanada_radical_II = classifierResults[1].radical;
			hanada_radical_III = classifierResults[2].radical;
		});
		glue.inMode("module/covGrantham1974DistanceCalculator", function() {
			var granthamResult = glue.command(["distance", replacementObj.refAa, replacementObj.replacementAa]).grantham1974DistanceResult;
			grantham_distance_double = granthamResult.distanceDouble;
			grantham_distance_int = granthamResult.distanceInt;
		});
		glue.inMode("module/covMiyata1979DistanceCalculator", function() {
			miyataDistance = glue.command(["distance", replacementObj.refAa, replacementObj.replacementAa]).miyata1979DistanceResult.distance;
		});
	}
	
	glue.command(["create", "custom-table-row", "cov_replacement", replacementObj.id]);
	glue.inMode("custom-table-row/cov_replacement/"+replacementObj.id, function() {
		var displayName = replacementObj.refAa+replacementObj.codonLabel+replacementObj.replacementAa;
		glue.command(["set", "field", "display_name", displayName]);
		glue.command(["set", "field", "codon_label", replacementObj.codonLabel]);		
		glue.command(["set", "field", "codon_label_int", parseInt(replacementObj.codonLabel)]);		
		glue.command(["set", "field", "reference_nt", replacementObj.refNt]);
		glue.command(["set", "field", "reference_aa", replacementObj.refAa]);
		glue.command(["set", "field", "replacement_aa", replacementObj.replacementAa]);
		glue.command(["set", "field", "num_seqs", replacementObj.memberSeqs.length]);
		if(replacementObj.parentFeature == "ORF_1a") {
			glue.command(["set", "field", "parent_feature", "ORF_1a"]);
			glue.command(["set", "field", "parent_codon_label", orf1aReplacements[replacementObj.refNt+":"+replacementObj.replacementAa].codonLabel]);
		}
		if(replacementObj.parentFeature == "ORF_1ab") {
			glue.command(["set", "field", "parent_feature", "ORF_1ab"]);
			glue.command(["set", "field", "parent_codon_label", orf1abReplacements[replacementObj.refNt+":"+replacementObj.replacementAa].codonLabel]);
		}
		
		glue.command(["set", "link-target", "variation", 
			"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+replacementObj.feature+
			"/variation/"+variationName]);
		if(classifyReplacement) {
			glue.command(["set", "field", "radical_hanada_category_i", hanada_radical_I]);
			glue.command(["set", "field", "radical_hanada_category_ii", hanada_radical_II]);
			glue.command(["set", "field", "radical_hanada_category_iii", hanada_radical_III]);
			glue.command(["set", "field", "grantham_distance_double", grantham_distance_double]);
			glue.command(["set", "field", "grantham_distance_int", grantham_distance_int]);
			glue.command(["set", "field", "miyata_distance", miyataDistance]);
		}
	});
	
	_.each(replacementObj.memberSeqs, function(memberObj) {
		var sourceName = memberObj["sequence.source.name"];
		var sequenceID = memberObj["sequence.sequenceID"];
		var linkObjId = replacementObj.id+":"+sequenceID;
		glue.command(["create", "custom-table-row", "cov_replacement_sequence", linkObjId]);
		glue.inMode("custom-table-row/cov_replacement_sequence/"+linkObjId, function() {
			glue.command(["set", "link-target", "cov_replacement", "custom-table-row/cov_replacement/"+replacementObj.id]);
			glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
		});
	});
});