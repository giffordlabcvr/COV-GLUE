

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
	"cov_glue_lineage",
	"cov_glue_lw_ratio",
	"gisaid_passage_details",
    "gisaid_specimen",
    "gisaid_patient_sex",
    "gisaid_patient_age",
    "gisaid_sequencing_technology",
    "gisaid_assembly_method",
    "gisaid_patient_status",
    "gisaid_lineage",
    "gisaid_clade",
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
    "num_unique_snps",
    "csr_coverage"
	]);

glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":tabularData, 
		"fileName":"tabular/covGlueSanitisedMetadata.txt"
	}});
});