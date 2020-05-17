glue.command(["multi-unset", "field", "sequence", "-w", "pang_lineage != null", "pang_lineage"]);
glue.command(["multi-unset", "field", "sequence", "-w", "pang_representative != null", "pang_representative"]);
glue.command(["multi-unset", "field", "sequence", "-w", "pang_masked_snps != null", "pang_masked_snps"]);

//load data from PANGOLIN project lineages table.
var pangSeqObjs;
glue.inMode("module/tabularUtilityCsv", function() {
	pangSeqObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/lineages.2020-05-07.csv"]));
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
		snps.push(newSnp);
	});
});

var processed = 0;

_.each(pangSeqObjs, function(pangSeqObj) {
	var name = pangSeqObj["name"];
	var representative = pangSeqObj.representative != null && pangSeqObj.representative.trim() == "1";
	if(!representative) {
		return;
	}
	var snps = isolateToSnps[name];
	var lineage = pangSeqObj.lineage.trim();
	// various rules to undo the transformation of isolate names that the lineages repo has done
	// within lineage representatives
	name = name.replace("Hong_Kong", "Hong Kong");
	name = name.replace("USA/VA-DCLS-00", "USA/VA-DCLS-0");
	name = name.replace("USA/UT-000", "USA/UT-0");
	
	var seqID;
	var sourceName;
	var gisaidIsolateMatches = 
		glue.tableToObjects(glue.command(["list", "sequence", 
			"-w", "source.name = 'cov-gisaid' and isolate = '"+name+"'"]));
	if(gisaidIsolateMatches.length > 0) {
		seqID = gisaidIsolateMatches[0]["sequenceID"]
		sourceName = gisaidIsolateMatches[0]["source.name"]
	}
	if(seqID == null && sourceName == null) {
		var cogukMatches;
		var bits = name.split("/");
		if(bits.length >= 3 && ["England", "Scotland", "Wales", "Northern_Ireland"].indexOf(bits[0]) >= 0) {
			var cogUkId = bits[1];
			var cogUkIdMatches = glue.tableToObjects(glue.command(["list", "sequence", 
				"-w", "source.name in ('cov-coguk', 'cov-coguk-refs') and (sequenceID = '"+cogUkId+"' or isolate = '"+name+"')"]));
			if(cogUkIdMatches.length > 0) {
				seqID = cogUkIdMatches[0]["sequenceID"]
				sourceName = cogUkIdMatches[0]["source.name"]
			}
		}
	}
	if(seqID == null && sourceName == null) {
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