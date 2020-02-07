glue.command(["multi-unset", "link-target", "variation", "cov_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_deletion_sequence", "cov_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_deletion_sequence", "sequence", "-a"]);
glue.command(["multi-delete", "cov_deletion", "-a"]);
glue.command(["multi-delete", "cov_deletion_sequence", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del%'"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del_detect%'"]);

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
_.each(featuresList, function(featureObj) {
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
		var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member"]));
		_.each(almtMemberObjs, function(almtMemberObj) {
			glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
				
				var memberDelObjs = glue.tableToObjects(glue.command(["variation", "scan", 
					"-r", comparisonRefName, "-f", featureObj.name, 
					"--whereClause", "name = 'cov_aa_del_detect:"+featureObj.name+"'", 
					"--excludeAbsent", "--showMatchesAsTable"]));
				
				_.each(memberDelObjs, function(memberDelObj) {
					var codonAligned = memberDelObj.deletionIsCodonAligned;
					var deletionID;
					if(codonAligned) {
						deletionID = featureObj.name+":ca:"+memberDelObj.refFirstCodonDeleted+":"+memberDelObj.refLastCodonDeleted;
					} else {
						throw new Error("Non-codon-aligned deletion! We could add this but what are the implications for protein translation elsewhere in the feature");
					}
					var deletionObj = deletionsSet[deletionID];
					if(deletionObj == null) {
						deletionObj = {
							id: deletionID,
							feature: featureObj.name,
							codonAligned: memberDelObj.deletionIsCodonAligned,
							startCodon: memberDelObj.refFirstCodonDeleted,
							endCodon: memberDelObj.refLastCodonDeleted,
							refNtStart: memberDelObj.refFirstNtDeleted,
							refNtEnd: memberDelObj.refLastNtDeleted,
							memberSeqs: []
						};
						deletionsSet[deletionID] = deletionObj;
					}
					deletionObj.memberSeqs.push(almtMemberObj);
				});
			});
		}); 
	});
});

_.each(_.values(deletionsSet), function(deletionObj) {
	glue.log("FINEST", "Creating deletion object", deletionObj);
	var variationName = "cov_aa_del:"+deletionObj.id;
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+deletionObj.feature, function() {
		glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidDeletion", 
			"--labeledCodon", deletionObj.startCodon, deletionObj.endCodon]);
	});
	
	glue.command(["create", "custom-table-row", "cov_deletion", deletionObj.id]);
	glue.inMode("custom-table-row/cov_deletion/"+deletionObj.id, function() {
		var displayName = deletionObj.startCodon+"-"+deletionObj.endCodon+"del";
		glue.command(["set", "field", "display_name", displayName]);
		glue.command(["set", "field", "start_codon", deletionObj.startCodon]);		
		glue.command(["set", "field", "start_codon_int", parseInt(deletionObj.startCodon)]);		
		glue.command(["set", "field", "end_codon", deletionObj.endCodon]);		
		glue.command(["set", "field", "end_codon_int", parseInt(deletionObj.endCodon)]);		
		glue.command(["set", "field", "codon_aligned", deletionObj.codonAligned]);		
		glue.command(["set", "field", "reference_nt_start", deletionObj.refNtStart]);		
		glue.command(["set", "field", "reference_nt_end", deletionObj.refNtEnd]);		
		glue.command(["set", "field", "num_seqs", deletionObj.memberSeqs.length]);
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
});