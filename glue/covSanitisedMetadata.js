

var tabularData = glue.command(["list", "sequence", 
	"-w", "source.name = 'cov-gisaid'",
	"sequenceID", 
	"source.name", 
	"length", 
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
	"gisaid_status_code", 
	"gisaid_status_comment", 
	"gisaid_originating_lab",
	"gisaid_submitting_lab",
	"gisaid_authors",
	"analyse_variation", 
	"include_in_ref_tree", 
    "num_ns",
    "longest_n_run",
    "num_bivalent_ambigs",
    "num_trivalent_ambigs",
    "num_hyphens",
    "num_whitespace",
    "num_other",
    "uncorrected_illegals",
    "num_unique_snps"
	]);

glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":tabularData, 
		"fileName":"tabular/covGlueSanitisedMetadata.txt"
	}});
});