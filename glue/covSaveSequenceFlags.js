

var sequenceFlags = glue.command(["list", "sequence", "source.name", "sequenceID", "include_in_ref_tree", "ref_tree_candidate",
	"analyse_aa_replacements", "analyse_aa_deletions", "analyse_aa_insertions", "cdhit_cluster", "is_l_lineage"]);

glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":sequenceFlags, 
		"fileName":"tabular/sequenceFlags.txt"
	}});
});