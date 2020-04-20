glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del_detect%'"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_ins_detect%'"]);

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));

// create "detection" variation objects to detect insertions and deletions in each of the features respectively.
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
_.each(featuresList, function(featureObj) {
	glue.inMode("reference/"+comparisonRefName+"/feature-location/"+featureObj.name, function() {
		var codonLabels= glue.getTableColumn(glue.command(["list", "labeled-codon"]), "codonLabel");
		var codonLabel1 = codonLabels[0];
		var codonLabelLast = codonLabels[codonLabels.length - 1];
		glue.command(["create", "variation", "cov_aa_ins_detect:"+featureObj.name, 
			"-t", "aminoAcidInsertion", "-c", codonLabel1, codonLabelLast]);
	});
});


var implausibleDeletionAlmtMembers = [];
var implausibleInsertionAlmtMembers = [];

glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
	var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member"]));
	_.each(almtMemberObjs, function(almtMemberObj) {
		glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
			
			var containsImplausibleDeletion = false;
			
			_.each(featuresList, function(featureObj) {
				var memberDelObjs = glue.tableToObjects(glue.command(["variation", "scan", 
					"-r", comparisonRefName, "-f", featureObj.name, 
					"--whereClause", "name = 'cov_aa_del_detect:"+featureObj.name+"'", 
					"--excludeAbsent", "--showMatchesAsTable"]));
				
				_.each(memberDelObjs, function(memberDelObj) {
					if(memberDelObj.deletedRefNts.length % 3 != 0) {
						containsImplausibleDeletion = true;
					}
				});
			});

			if(containsImplausibleDeletion) {
				implausibleDeletionAlmtMembers.push(almtMemberObj);
			}
			
			var containsImplausibleInsertion = false;
			
			_.each(featuresList, function(featureObj) {
				var memberInsObjs = glue.tableToObjects(glue.command(["variation", "scan", 
					"-r", comparisonRefName, "-f", featureObj.name, 
					"--whereClause", "name = 'cov_aa_ins_detect:"+featureObj.name+"'", 
					"--excludeAbsent", "--showMatchesAsTable"]));
				
				_.each(memberInsObjs, function(memberInsObj) {
					if(memberInsObj.insertedQryNts.length % 3 != 0) {
						containsImplausibleInsertion = true;
					}
				});
			});

			if(containsImplausibleInsertion) {
				implausibleInsertionAlmtMembers.push(almtMemberObj);
			}

		});
	}); 
});

_.each(implausibleDeletionAlmtMembers, function(almtMemberObj) {
	glue.log("WARNING", "implausible deletion detected for sequence "+almtMemberObj["sequence.sequenceID"]+", it will be excluded from variation analysis");
	glue.inMode("sequence/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
		glue.command(["set", "field", "analyse_variation", "false"]);
	});
});

_.each(implausibleInsertionAlmtMembers, function(almtMemberObj) {
	glue.log("WARNING", "implausible insertion detected for sequence "+almtMemberObj["sequence.sequenceID"]+", it will be excluded from variation analysis");
	glue.inMode("sequence/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
		glue.command(["set", "field", "analyse_variation", "false"]);
	});
});

