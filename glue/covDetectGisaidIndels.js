glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del_detect%'"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_ins_detect%'"]);

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


var refName = "REF_MASTER_WUHAN_HU_1";
_.each(featuresList, function(featureObj) {
	glue.inMode("reference/"+refName+"/feature-location/"+featureObj.name, function() {
		var codonLabels= glue.getTableColumn(glue.command(["list", "labeled-codon"]), "codonLabel");
		var codonLabel1 = codonLabels[0];
		var codonLabel2 = codonLabels[1];
		var codonLabelPenultimate = codonLabels[codonLabels.length - 2];
		var codonLabelLast = codonLabels[codonLabels.length - 1];
		glue.command(["create", "variation", "cov_aa_del_detect:"+featureObj.name, 
			"-t", "aminoAcidDeletion", "-c", codonLabel2, codonLabelPenultimate]);
		glue.command(["create", "variation", "cov_aa_ins_detect:"+featureObj.name, 
			"-t", "aminoAcidInsertion", "-c", codonLabel1, codonLabelLast]);
	});
});

var someDeletionsFound = false;
var someInsertionsFound = false;

glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
	_.each(featuresList, function(featureObj) {
		var deletionsFound = glue.tableToObjects(
				glue.command(["variation", "member", "scan", 
					"-r", refName, "-f", featureObj.name, "-v", "cov_aa_del_detect:"+featureObj.name, 
					"--excludeAbsent"]));
		featureObj.deletionsFound = deletionsFound;
		glue.logInfo("Deletions found for "+featureObj.displayName, deletionsFound);
		if(deletionsFound.length > 0) {
			someDeletionsFound = true;
		}
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

if(someDeletionsFound || someInsertionsFound) {
	throw new Error("Some indels found!");
}
