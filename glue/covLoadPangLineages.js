glue.command(["multi-unset", "field", "sequence", "-w", "pang_lineage != null", "pang_lineage"]);
glue.command(["multi-unset", "field", "sequence", "-w", "pang_representative != null", "pang_representative"]);
glue.command(["multi-unset", "field", "sequence", "-w", "pang_masked_snps != null", "pang_masked_snps"]);

// load data from PANGOLIN project lineages table.
// this will populate pang_... fields for
// (1) Any GISAID sequences in the lineages table, if they can be found in COV-GLUE
//     Although note, UK sequences in the lineages table do not have GISAID IDs even if they are on GISAID.
// (2) Any COGUK sequences which are also marked as lineage representatives.
var pangSeqObjs;
glue.inMode("module/tabularUtilityCsv", function() {
	pangSeqObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/lineages.metadata.csv"]));
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

_.each(pangSeqObjs, function(pangSeqObj) {
	var representative = pangSeqObj.representative != null && pangSeqObj.representative.trim() == "1";
	var name = pangSeqObj["name"];
	var lineage = pangSeqObj.lineage.trim();
	var snps = isolateToSnps[name];

	var seqID;
	var sourceName;

	var gisaidID = pangSeqObj["GISAID ID"];
	if(gisaidID != null) {
		// first try cov-gisaid
		var gisaidIdMatches = 
			glue.tableToObjects(glue.command(["list", "sequence", 
				"-w", "source.name = 'cov-gisaid' and sequenceID = '"+gisaidID+"'"]));
		if(gisaidIdMatches.length > 0) {
			seqID = gisaidIdMatches[0]["sequenceID"]
			sourceName = gisaidIdMatches[0]["source.name"]
		}
	} else {
		if(!representative) { // no GISAID ID, not representative
			if(pangSeqObj["country"] == "UK") {
				var transformedName = name.replace("Northern_Ireland", "Northern Ireland");
				var gisaidIsolateMatches = 
					glue.tableToObjects(glue.command(["list", "sequence", 
						"-w", "source.name = 'cov-gisaid' and m49_country.id = 'GBR' and isolate = '"+transformedName+"'"]));
				if(gisaidIsolateMatches.length > 0) {
					seqID = gisaidIsolateMatches[0]["sequenceID"]
					sourceName = gisaidIsolateMatches[0]["source.name"]
				} else {
					return;
				}
			} else {
				return; 
			}
		} else {
			// no GISAID ID, representative
			var cogukMatches;
			var bits = name.split("/");
			if(bits.length >= 3 && ["England", "Scotland", "Wales", "Northern_Ireland"].indexOf(bits[0]) >= 0) {
				var cogUkId = bits[1];
				// first try 'cov-coguk' source. This will exist when we are updating refs but not in the main build
				var cogUkIdMatches = glue.tableToObjects(glue.command(["list", "sequence", 
					"-w", "source.name = 'cov-coguk' and (sequenceID = '"+cogUkId+"' or isolate = '"+name+"')"]));
				if(cogUkIdMatches.length > 0) {
					seqID = cogUkIdMatches[0]["sequenceID"]
					sourceName = cogUkIdMatches[0]["source.name"]
				}
				// As fallback use 'cov-coguk-refs' source. This will exist when in the main build
				if(seqID == null && sourceName == null) {
					var cogUkRefsIdMatches = glue.tableToObjects(glue.command(["list", "sequence", 
						"-w", "source.name = 'cov-coguk-refs' and (sequenceID = '"+cogUkId+"' or isolate = '"+name+"')"]));
					if(cogUkRefsIdMatches.length > 0) {
						seqID = cogUkRefsIdMatches[0]["sequenceID"]
						sourceName = cogUkRefsIdMatches[0]["source.name"]
					}
				}
			}
		}
	}
	

	if(seqID == null && sourceName == null) {
		if(!representative) {
			return; // not a GISAID sequence we could find, not representative
		}
		// should only drop through to this bit when we are updating to a new version of the lineages CSV.
		// if the COGUK-COV-GLUE extension is loaded, it will associate sequences. 
		// but we then need to export them back out again using covExportCogUkRefs.js
		var cogukExtendedMatches;
		var bits = name.split("/");
		var cogUkIdMatches = glue.tableToObjects(glue.command(["list", "sequence", 
			"-w", "source.name = 'cov-coguk' and coguk_metadata.secondary_identifier = 'hCoV-19/"+name+"'"]));
		if(cogUkIdMatches.length > 0) {
			seqID = cogUkIdMatches[0]["sequenceID"]
			sourceName = cogUkIdMatches[0]["source.name"]
		} else {
			glue.logInfo("bits", bits);
		}
	}
	if(seqID != null && sourceName != null) {
		glue.inMode("sequence/"+sourceName+"/"+seqID, function() {
			glue.command(["set", "field", "pang_lineage", lineage]);
			glue.command(["set", "field", "pang_representative", representative]);
			if(snps != null && snps.length > 0) {
				glue.command(["set", "field", "pang_masked_snps", snps.join(",")]);
			}
		});
	} else {
		glue.log("WARNING", "No GISAID/COGUK sequence found for "+name+", representative of lineage "+lineage);
	} 
	processed++;
	if(processed % 250 == 0) {
		glue.log("FINEST", "Loaded lineage for "+processed+" sequences");
		glue.command(["new-context"]);
	}
});

glue.log("FINEST", "Loaded lineage for "+processed+" sequences");
glue.command(["new-context"]);
