// This is not part of any build but is to be run on an ad-hoc basis.
// Possible reasons for running it: 
// (1) the global lineages have been updated
// (2) a number of COGUK sequences marked as lineage representatives, previously only on COGUK,
//     are now on GISAID.
// It assumes the COGUK-COV-GLUE extension is loaded.

// It ouputs these updates:
// (a) The cov-coguk-refs source, containing those lineage reps only found in COG-UK
//
// (b) A file gisaidLineageRepresentatives.tsv, with one row per GISAID lineage rep, containing:
//
//     sequenceID,
//     pang_lineage,
//     pang_representative,
//     pang_masked_snps,
//     gisaid_authors
//     gisaid_originating_lab
//     gisaid_submitting_lab
//
// (c) A file cogukLineageRepresentatives.tsv, with one row per COG-UK-only lineage reps, 
//     containing:
//     
//     sequenceID,
//     pang_lineage,
//     pang_representative,
//     pang_masked_snps,
//     isolate,
//     place_sampled,
//     m49_country.id,
//     collection_month_day,
//     collection_month,
//     collection_year

// Steps:
// 1. Delete cov-coguk-refs files on disk.
// 2. Load the files lineages.metadata.csv and singleton.csv for the snps.
// 3. For each row that is marked as a lineage rep:
//    If row does not have a GISAID ID
//    -- Try to establish GISAID ID by looking up virus name, making substitutions where necessary
//    If row now has a GISAID ID
//    -- lookup authorship data from GISAID
//    -- output a line to gisaidLineageRepresentatives.tsv
//    Else
//    -- Find COG-UK sequence based on virus name
//    -- If it can't be found, throw an error.
//    -- output a line to cogukLineageRepresentatives.tsv
//    -- create a sequence in source cov-coguk-refs
// 4. Save cov-coguk-refs FASTA files
// 5. Save file gisaidLineageRepresentatives.tsv
// 6. Save file cogukLineageRepresentatives.tsv

function handleNull(value) {
	if(value == null) {
		return "";
	}
	return value;
}

glue.command(["file-util", "make-directories", "sources/cov-coguk-refs"]);

var fileObjs = glue.tableToObjects(
		glue.command(["file-util", "list-files", "-d", "sources/cov-coguk-refs"]));
_.each(fileObjs, function(fileObj) {
	if(fileObj.fileName.indexOf(".fasta") >= 0) {
		glue.command(["file-util", "delete-file", "sources/cov-coguk-refs/"+fileObj.fileName]);
	}
});

var gisaidRefsMetadata = { "listResult": {
	"column": ["sequenceID", "pang_lineage", "pang_representative", "pang_masked_snps", "gisaid_authors", "gisaid_originating_lab", "gisaid_submitting_lab"],
	"row" : []
} }

var cogukRefsMetadata = { "listResult": {
	"column": ["sequenceID", "pang_lineage", "pang_representative", "pang_masked_snps", "isolate", "place_sampled", "m49_country.id", 
		 "collection_month_day", "collection_month", "collection_year"],
	"row" : []
} }


var lineageMetadataObjs;
glue.inMode("module/tabularUtilityCsv", function() {
	lineageMetadataObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/lineages.metadata.csv"]));
});

var isolateToSnps = {};

//load data from PANGOLIN project lineages table.
glue.inMode("module/tabularUtilityCsv", function() {
	var singletonRows = glue.tableToObjects(glue.command(["load-tabular", "tabular/singletons.csv"]));
	_.each(singletonRows, function(singletonRow) {
		var taxon = singletonRow.taxon.trim();
		var snps = isolateToSnps[taxon];
		if(snps == null) {
			snps = [];
			isolateToSnps[taxon] = snps;
		}
		var oldSnp = singletonRow.snp.trim();
		var newSnp = oldSnp.substring(oldSnp.length-2, oldSnp.length-1)+oldSnp.substring(0, oldSnp.length-2)+oldSnp.substring(oldSnp.length-1, oldSnp.length);
		if(!newSnp.endsWith("-")) { 
			// no point recording these snps as we would replace a "-" with a "-" effectively.
			// also there are some misaligned deletions here
			snps.push(newSnp);
		}
	});
});


var processed = 0;

_.each(lineageMetadataObjs, function(lineageMetadataObj) {
	var representative = lineageMetadataObj.representative != null && lineageMetadataObj.representative.trim() == "1";
	if(!representative) {
		return
	}
	var name = lineageMetadataObj["name"];
	var lineage = lineageMetadataObj.lineage.trim();
	var snps = isolateToSnps[name];
	if(snps == null) {
		snps = [];
	}
	var seqID;
	var sourceName;
	var gisaidID = lineageMetadataObj["GISAID ID"];
	if(gisaidID != null) {
		// first try cov-gisaid
		var gisaidIdMatches = 
			glue.tableToObjects(glue.command(["list", "sequence", 
				"-w", "source.name = 'cov-gisaid' and sequenceID = '"+gisaidID+"'"]));
		if(gisaidIdMatches.length > 0) {
			seqID = gisaidIdMatches[0]["sequenceID"]
			sourceName = gisaidIdMatches[0]["source.name"]
		} else {
			throw new Error("lineages.metadata.csv representative row has GISAID ID = '"+gisaidID+"' but sequence could not be found in cov-gisaid");
		}
	} else {
		// no GISAID ID
		var cogukMatches;
		var bits = name.split("/");
		if(bits.length >= 3 && ["England", "Scotland", "Wales", "Northern_Ireland"].indexOf(bits[0]) >= 0) {
			// UK sequence
			var name2 = name.replace("_", " ");
			// try GISAID lookup based on isolate
			var gisaidIsolateMatches = 
				glue.tableToObjects(glue.command(["list", "sequence", 
					"-w", "source.name = 'cov-gisaid' and isolate in ('"+name2+"', '"+name+"')"]));
			if(gisaidIsolateMatches.length > 0) {
				seqID = gisaidIsolateMatches[0]["sequenceID"]
				sourceName = gisaidIsolateMatches[0]["source.name"]
			} else {
				// try COGUK lookup based on COGUK ID
				var cogUkId = bits[1];
				var cogUkIdMatches = 
					glue.tableToObjects(glue.command(["list", "sequence", 
						"-w", "source.name = 'cov-coguk' and sequenceID = '"+cogUkId+"'"]));
				if(cogUkIdMatches.length > 0) {
					seqID = cogUkIdMatches[0]["sequenceID"]
					sourceName = cogUkIdMatches[0]["source.name"]
				} else {
					// try COGUK lookup based on isolate
					var cogUkIsolateMatches = 
						glue.tableToObjects(glue.command(["list", "sequence", 
							"-w", "source.name = 'cov-coguk' and isolate in ('"+name2+"', '"+name+"')"]));
					if(cogUkIsolateMatches.length > 0) {
						seqID = cogUkIsolateMatches[0]["sequenceID"]
						sourceName = cogUkIsolateMatches[0]["source.name"]
					} else {
						// try secondary identifier
						var cogUkSIMatches = glue.tableToObjects(glue.command(["list", "sequence", 
							"-w", "source.name = 'cov-coguk' and coguk_metadata.secondary_identifier in ('hCoV-19/"+name+"', 'hCoV-19/"+name2+"')"]));
						if(cogUkSIMatches.length > 0) {
							seqID = cogUkSIMatches[0]["sequenceID"]
							sourceName = cogUkSIMatches[0]["source.name"]
						} else {
							throw new Error("UK lineages.metadata.csv representative row with name = '"+name+"' not found in COG-UK dataset");						
						}
					}
				}

			}
		} else {
			throw new Error("Non-UK lineages.metadata.csv representative row with name = '"+name+"' has no GISAID ID");
		}
	}

	if(sourceName == 'cov-coguk') {
		var fullSeqObj = glue.tableToObjects(glue.command(["list", "sequence", 
			"-w", "source.name = '"+sourceName+"' and sequenceID = '"+seqID+"'", 
			"isolate",
			"place_sampled",
			"m49_country.id",
			"collection_month_day",
			"collection_month",
			"collection_year"]))[0];
		cogukRefsMetadata.listResult.row.push(
			{ "value": [
				seqID, 
				lineage,
				true,
				snps.join(","), 
				handleNull(fullSeqObj["isolate"]),
				handleNull(fullSeqObj["place_sampled"]),
				handleNull(fullSeqObj["m49_country.id"]),
				handleNull(fullSeqObj["collection_month_day"]),
				handleNull(fullSeqObj["collection_month"]),
				handleNull(fullSeqObj["collection_year"]),
			] }
		);
		var seqNts;
		glue.inMode("sequence/"+sourceName+"/"+seqID, function() {
			seqNts = glue.command(["show", "nucleotides"]).nucleotidesResult.nucleotides;
		});
		
		glue.inMode("module/covFastaUtility", function() {
			seqNts = glue.command({
				"save-nucleotide-fasta": {
					"fastaCommandDocument": {
						"nucleotideFasta": {
							"sequences": [
								{
									"id": seqID,
									"sequence": seqNts
								}
							]
						}
					},
					"outputFile": "sources/cov-coguk-refs/"+seqID+".fasta"
				}
			});
		});
	} else if(sourceName == 'cov-gisaid') {
		var urlPath = "/acknowledgement/"+
			seqID.substring(10,12)+"/"+seqID.substring(12,14)+"/"+seqID+".json";
		var authorship;
		glue.inMode("module/covEpicovAcknowledgementHttpRunner", function() {
			authorship = JSON.parse(glue.command(["get", urlPath]).httpRunnerResult.entityAsString);
		});
		
		var gisaid_authors = authorship["covv_authors"];
		var gisaid_originating_lab = authorship["covv_orig_lab"];
		var gisaid_submitting_lab = authorship["covv_subm_lab"];
		
		gisaidRefsMetadata.listResult.row.push({ "value": [
			seqID, 
			lineage, 
			true,
			snps.join(","), 
			gisaid_authors,
			gisaid_originating_lab,
			gisaid_submitting_lab,
		] } );

	} else {
		throw new Error("Don't know what to do with sourceName "+sourceName);
	}
	processed++;
	if(processed % 10 == 0) {
		glue.logInfo("Processed "+processed+" lineage representatives.");
	}
});


glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":gisaidRefsMetadata, 
		"fileName":"tabular/gisaidLineageRepresentatives.tsv"
	}});
});

glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":cogukRefsMetadata, 
		"fileName":"tabular/cogukLineageRepresentatives.tsv"
	}});
});

