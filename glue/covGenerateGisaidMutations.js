glue.command(["multi-unset", "link-target", "variation", "cov_mutation", "-a"]);
glue.command(["multi-unset", "link-target", "cov_mutation_sequence", "cov_mutation", "-a"]);
glue.command(["multi-unset", "link-target", "cov_mutation_sequence", "sequence", "-a"]);
glue.command(["multi-delete", "cov_mutation", "-a"]);
glue.command(["multi-delete", "cov_mutation_sequence", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like 'covmut%'"]);

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
var mutationsSet = {};
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
					if(memberAaObj.definiteAas != null && memberAaObj.definiteAas != "") {
						refAaObj = refAaObjsMap[memberAaObj.codonLabel];
						if(refAaObj != null && refAaObj.definiteAas != null && refAaObj.definiteAas != "" && 
								refAaObj.definiteAas != memberAaObj.definiteAas) {
							var refAas = refAaObj.definiteAas.split('');
							var memberAas = memberAaObj.definiteAas.split('');
							_.each(refAas, function(refAa) {
								_.each(memberAas, function(memberAa) {
									if(refAa != memberAa) {
										var mutationID = featureObj.name+":"+refAa+":"+memberAaObj.codonLabel+":"+memberAa;
										var mutationObj = mutationsSet[mutationID];
										if(mutationObj == null) {
											mutationObj = {
												id: mutationID,
												feature: featureObj.name,
												codonLabel: memberAaObj.codonLabel,
												refAa: refAa,
												mutationAa: memberAa,
												memberSeqs: []
											};
											mutationsSet[mutationID] = mutationObj;
										}
										mutationObj.memberSeqs.push(almtMemberObj);
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

_.each(_.values(mutationsSet), function(mutationObj) {
	glue.log("FINEST", "Creating mutation object", mutationObj);
	var variationName = "covmut:"+mutationObj.id;
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+mutationObj.feature, function() {
		glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidSimplePolymorphism", 
			"--labeledCodon", mutationObj.codonLabel, mutationObj.codonLabel]);
		glue.inMode("variation/"+variationName, function() {
			glue.command(["set", "metatag", "SIMPLE_AA_PATTERN", mutationObj.mutationAa]);
			glue.command(["set", "metatag", "MIN_COMBINED_TRIPLET_FRACTION", 0.25]);
		});
	});
	
	glue.command(["create", "custom-table-row", "cov_mutation", mutationObj.id]);
	glue.inMode("custom-table-row/cov_mutation/"+mutationObj.id, function() {
		var displayName = mutationObj.refAa+mutationObj.codonLabel+mutationObj.mutationAa;
		glue.command(["set", "field", "display_name", displayName]);
		glue.command(["set", "field", "codon_label", mutationObj.codonLabel]);		
		glue.command(["set", "field", "reference_aa", mutationObj.refAa]);
		glue.command(["set", "field", "mutation_aa", mutationObj.mutationAa]);
		glue.command(["set", "link-target", "variation", 
			"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+mutationObj.feature+
			"/variation/"+variationName]);
	});
	
	_.each(mutationObj.memberSeqs, function(memberObj) {
		var sourceName = memberObj["sequence.source.name"];
		var sequenceID = memberObj["sequence.sequenceID"];
		var linkObjId = mutationObj.id+":"+sequenceID;
		glue.command(["create", "custom-table-row", "cov_mutation_sequence", linkObjId]);
		glue.inMode("custom-table-row/cov_mutation_sequence/"+linkObjId, function() {
			glue.command(["set", "link-target", "cov_mutation", "custom-table-row/cov_mutation/"+mutationObj.id]);
			glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
		});
	});
});