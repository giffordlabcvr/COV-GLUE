glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_ins_detect%'"]);

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName"]));

var refName = "REF_MASTER_WUHAN_HU_1";
_.each(featuresList, function(featureObj) {
	glue.inMode("reference/"+refName+"/feature-location/"+featureObj.name, function() {
		var codonLabels= glue.getTableColumn(glue.command(["list", "labeled-codon"]), "codonLabel");
		var codonLabel1 = codonLabels[0];
		var codonLabelLast = codonLabels[codonLabels.length - 1];
		glue.command(["create", "variation", "cov_aa_ins_detect:"+featureObj.name, 
			"-t", "aminoAcidInsertion", "-c", codonLabel1, codonLabelLast]);
	});
});

var someInsertionsFound = false;

glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
	_.each(featuresList, function(featureObj) {
		var insertionsFound = glue.tableToObjects(
				glue.command(["variation", "member", "scan", 
					"-r", refName, "-f", featureObj.name, "-v", "cov_aa_ins_detect:"+featureObj.name, 
					"--excludeAbsent"]));
		featureObj.insertionsFound = insertionsFound;
		glue.logInfo("Insertions found for "+featureObj.displayName, insertionsFound);
		if(insertionsFound.length > 0) {
			someInsertionsFound = true;
		}
	});
});

if(someInsertionsFound) {
	throw new Error("Some insertions found!");
}
