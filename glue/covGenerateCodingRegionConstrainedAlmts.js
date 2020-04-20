// script to generate constrained nucleotide and protein alignments for coding regions.

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", 
			"featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name"]));

_.each(featuresList, function(featureObj) {
	glue.inMode("module/covFastaProteinAlignmentExporterSeqIdOnly", function() {
		glue.command(["export", "AL_GISAID_CONSTRAINED", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", featureObj.name,
			"-a", 
			"-o", "alignments/covConstrainedAlignment_"+featureObj.name+".faa"]);
	});
	glue.inMode("module/covFastaAlignmentExporterSeqIdOnly", function() {
		glue.command(["export", "AL_GISAID_CONSTRAINED", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", featureObj.name,
			"-a",
			"-o", "alignments/covConstrainedAlignment_"+featureObj.name+".fna"]);
	});
});