glue.command(["multi-unset", "link-target", "variation", "cov_replacement", "-a"]);
glue.command(["multi-unset", "link-target", "cov_replacement_sequence", "cov_replacement", "-a"]);
glue.command(["multi-unset", "link-target", "cov_replacement_sequence", "sequence", "-a"]);
glue.command(["multi-delete", "cov_replacement", "-a"]);
glue.command(["multi-delete", "cov_replacement_sequence", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_rpl%'"]);

var featuresList = [
    { name: "E",
      displayName: "E" },
    { name: "M",
      displayName: "M" },
    { name: "N",
      displayName: "N" },
    { name: "ORF_10",
      displayName: "ORF 10" },
    { name: "ORF_1ab",
      displayName: "ORF 1ab" },
    { name: "ORF_3a",
      displayName: "ORF_3a" },
    { name: "ORF_6",
      displayName: "ORF 6" },
    { name: "ORF_7a",
      displayName: "ORF 7a" },
    { name: "ORF_8",
      displayName: "ORF 8" },
    { name: "S",
      displayName: "S" }];


var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
var replacementsSet = {};
_.each(featuresList, function(featureObj) {
	var refAaObjsMap = {};
	glue.inMode("reference/"+comparisonRefName+"/feature-location/"+featureObj.name, function() {
		var refAaObjs = glue.tableToObjects(glue.command(["amino-acid"]));
		_.each(refAaObjs, function(refAaObj) {
			refAaObjsMap[refAaObj.codonLabel] = refAaObj;
		});
	});
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
		var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member"]));
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
												codonLabel: memberAaObj.codonLabel,
												refNt: memberAaObj.relRefNt,
												refAa: refAa,
												replacementAa: memberAa,
												memberSeqs: []
											};
											replacementsSet[replacementID] = replacementObj;
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
		glue.command(["set", "link-target", "variation", 
			"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+replacementObj.feature+
			"/variation/"+variationName]);
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