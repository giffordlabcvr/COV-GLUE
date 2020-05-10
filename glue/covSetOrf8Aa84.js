glue.command(["multi-unset", "field", "sequence", "-a", "orf_8_aa_84"]);

var proteinAlmt;

glue.inMode("module/covFastaProteinAlignmentExporterSourcePlusSeqId", function() {
	proteinAlmt = glue.command(["export", "AL_GISAID_CONSTRAINED", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", "ORF_8", "-l", "84", "84", 
			"-w", "sequence.include_in_ref_tree = true", 
			"-p"]);
});

_.each(proteinAlmt.aminoAcidFasta.sequences, function(membObj) {
	if(membObj.sequence != "-") {
		glue.inMode("sequence/"+membObj.id, function() {
			glue.command(["set", "field", "orf_8_aa_84", membObj.sequence]);
		});
	}
});
