glue.command(["delete", "module", "-w", "name like 'cov_region_selector_%'"]);

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", 
			"featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "max_codon_number"]));

_.each(featuresList, function(featureObj) {
	var selectorModuleName = "cov_region_selector_"+featureObj.name;
	glue.command(["create", "module", "--moduleType", "alignmentColumnsSelector", selectorModuleName]);
	glue.inMode("module/"+selectorModuleName, function() {
		glue.command(["set", "property", "relRefName", "REF_MASTER_WUHAN_HU_1"]);
		if(featureObj.name == 'ORF_1ab') {
			// for ORF 1ab, remove codons 4401 and 4402 which are involved in the ribosomal slippage
			glue.command(["add", "region-selector", "-f", featureObj.name, "-a", "-l", "1", "4400"]);
			glue.command(["add", "region-selector", "-f", featureObj.name, "-a", "-l", "4403", featureObj.max_codon_number]);
		} else if(featureObj.name == 'NSP12') {
			// for nsp12 (RdRp), remove codons 9 and 10 which are involved in the ribosomal slippage
			glue.command(["add", "region-selector", "-f", featureObj.name, "-a", "-l", "1", "8"]);
			glue.command(["add", "region-selector", "-f", featureObj.name, "-a", "-l", "11", featureObj.max_codon_number]);
		} else {
			glue.command(["add", "region-selector", "-f", featureObj.name, "-a"]);
		}
	});
	glue.inMode("module/covFastaProteinAlignmentExporterSeqIdOnly", function() {
		glue.command(["export", "AL_GISAID_UNCONSTRAINED", "-s", selectorModuleName, "-a", "-o", "alignments/cov_"+featureObj.name+".faa"]);
	});
	glue.inMode("module/covFastaAlignmentExporterSeqIdOnly", function() {
		glue.command(["export", "AL_GISAID_UNCONSTRAINED", "-s", selectorModuleName, "-a", "-o", "alignments/cov_"+featureObj.name+".fna"]);
	});
});