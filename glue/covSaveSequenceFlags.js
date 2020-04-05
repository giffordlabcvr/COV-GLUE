

var sequenceFlags = glue.command(["list", "sequence", "source.name", "sequenceID", "include_in_ref_tree", "ref_tree_candidate",
	"analyse_variation", "cdhit_cluster", "orf_8_aa_84"]);

glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":sequenceFlags, 
		"fileName":"tabular/sequenceFlags.txt"
	}});
});