

var tabularData = glue.command(["list", "sequence", 
	"source.name", 
	"sequenceID", 
	"length", 
	// isolate == "virus name" in GISAID
    "isolate",
    "host_species",
    "collection_year",
    "collection_month",
    "collection_month_day",
    "submission_year",
    "submission_month",
    "submission_month_day",
	"place_sampled",
	"m49_country.display_name",
	"m49_country.m49_sub_region.display_name",
	"m49_country.m49_region.display_name",
	"gisaid_passage_details",
	// GISAID QC mechanism: warn / info / normal
	"gisaid_status_code", 
	// GISAID QC comment
	"gisaid_status_comment", 
	"gisaid_originating_lab",
	"gisaid_submitting_lab",
	"gisaid_authors",
	// true if various QC and metadata criteria are met
	// if true, CoV-GLUE includes sequence in its analysis of AA replacements etc.
	"analyse_variation", 
	// ref_tree_candidate is a subset of analyse_variation, 
	// if true, sequence may be selected to appear in the CoV-GLUE reference tree
	"ref_tree_candidate", 
	// if true, appears in the CoV-GLUE reference tree
	"include_in_ref_tree", 
	// for those sequences where include_in_ref_tree == true, 
	// number of cluster which CD-HIT-EST assigned the sequence to
	"cdhit_cluster", 
	// most coarse-grained lineage assignment for the pandemic.
	// if true, belongs to the clade characterised by amino acid L at ORF8 codon 84
	// if false, belongs to the clade characterised by amino acid S at ORF8 codon 84
	"is_l_lineage",
	
	]);

glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":tabularData, 
		"fileName":"tabular/covGlueSanitisedMetadata.txt"
	}});
});