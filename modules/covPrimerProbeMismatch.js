function reportSingleFasta(fastaFilePath) {
	
	var htmlFilePath = fastaFilePath.substring(0, fastaFilePath.lastIndexOf("."))+".html";
	var fastaDoc;
	glue.inMode("module/covFastaUtility", function() {
		fastaDoc = glue.command(["load-nucleotide-fasta", fastaFilePath]);
	});
	if(fastaDoc.nucleotideFasta.sequences.length > 1) {
		throw new Error("FASTA file can have only one sequence");
	}
	if(fastaDoc.nucleotideFasta.sequences.length == 0) {
		throw new Error("FASTA file has no sequences");
	}
	var queryName = fastaDoc.nucleotideFasta.sequences[0].id;
	var queryNucleotides = fastaDoc.nucleotideFasta.sequences[0].sequence;
	var targetRefName = "REF_MASTER_WUHAN_HU_1";
	var queryToTargetRefSegs = generateQueryToTargetRefSegs(targetRefName, queryNucleotides);
	
	var projectVersion = 
		glue.command(["show","setting","project-version"]).projectShowSettingResult.settingValue;
	var glueVersion = 
		glue.command(["glue-engine","show-version"]).glueEngineShowVersionResult.glueEngineVersion;
	
	var rasVariationMatchDocument;
	
	var matchResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_match:%'");
	//var matchResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_match:RdRP_SARSr-P1'");
	var insertionResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_insertion:%'");
	var deletionResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_deletion:%'");

	//glue.log("FINEST", "matchResults", matchResults);
	//glue.log("FINEST", "insertionResults", insertionResults);
	//glue.log("FINEST", "deletionResults", deletionResults);

	var vNameToMatchResult = {};
	_.each(matchResults, function(mr) {vNameToMatchResult[mr.variationName] = mr});
	var vNameToInsertionResult = {};
	_.each(insertionResults, function(ir) {vNameToInsertionResult[ir.variationName] = ir});
	var vNameToDeletionResult = {};
	_.each(deletionResults, function(dr) {vNameToDeletionResult[dr.variationName] = dr});
	
	var assayObjs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_primer_probe_assay", "id", "display_name", "organisation", "url"]));
	
	_.each(assayObjs, function(assayObj) {
		glue.inMode("custom-table-row/cov_primer_probe_assay/"+assayObj.id, function() {
			assayObj.ppObjs = glue.tableToObjects(glue.command(["list", "link-target", "cov_primer_probe", 
				"id", "display_name", "sequence_fwd", "sequence_fwd_regex", "sequence_rev", "sequence_rev_regex", 
				"ref_hits", "sequence_to_scan", "fwd_orientation", "ref_start", "ref_end", "length", 
				"seq_match.name", "seq_insertion.name", "seq_deletion.name"]));
		});
		_.each(assayObj.ppObjs, function(ppObj) {
			ppObj.seqMatchResult = vNameToMatchResult[ppObj["seq_match.name"]];
			ppObj.seqInsertionResult = vNameToInsertionResult[ppObj["seq_insertion.name"]];
			ppObj.seqDeletionResult = vNameToDeletionResult[ppObj["seq_deletion.name"]];
		});
		assayObj.anyProblem = false;
		// map of primer/probe display name to list of string issues.
		assayObj.ppDnToIssues = {};
		_.each(assayObj.ppObjs, function(ppObj) {
			var issues = [];
			assayObj.ppDnToIssues[ppObj.display_name] = issues;
			if(ppObj.seqMatchResult.sufficientCoverage == false) {
				
			}
		});		
	});
	
	var resultDoc = {
		fastaFilePath: fastaFilePath,
		sequenceID: queryName,
		assays: assayObjs,
		projectVersion: projectVersion,
		glueVersion: glueVersion,
		date: todaysDate(),
	};
	
	glue.log("FINEST", "resultDoc", resultDoc);
}

function todaysDate() {
	var today = new Date();
	var dd = today.getDate();
	var mm = today.getMonth()+1; // January is 0!
	var yyyy = today.getFullYear();
	if(dd < 10) {
	    dd = '0'+dd
	} 
	if(mm < 10) {
	    mm = '0'+mm
	} 
	return dd + '/' + mm + '/' + yyyy;
}

function variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, whereClause) {
	var results;
	glue.inMode("module/covFastaSequenceReporter", function() {
		results = glue.command({
			"string-plus-alignment" :{
				"variation":{
					"scan":{
						"fastaString":queryNucleotides,
						"queryToTargetSegs": {
							queryToTargetSegs: {
								alignedSegment: queryToTargetRefSegs
							}
						},
						"whereClause": whereClause,
						"targetRefName":targetRefName,
						"relRefName":targetRefName,
						"linkingAlmtName":"AL_GISAID_UNCONSTRAINED",
						"featureName":"whole_genome",
						"descendentFeatures":"false",
						"excludeAbsent":"false",
						"excludeInsufficientCoverage":"false",
						"showMatchesAsDocument":"true",
						"showMatchesAsTable":"false"
					}
				}
			}
		}).variationScanMatchCommandResult.variations;
	});
	return results;
}

function generateQueryToTargetRefSegs(targetRefName, nucleotides) {
	var alignerModule;
	glue.inMode("module/covFastaSequenceReporter", function() {
		alignerModule = glue.command(["show", "property", "alignerModuleName"]).moduleShowPropertyResult.propertyValue;
	});
	var alignResult;
	glue.inMode("module/"+alignerModule, function() {
		alignResult = glue.command({align: {
				referenceName: targetRefName,
				sequence: [
				    { 
				    	queryId: "query", 
				    	nucleotides: nucleotides
				    }
				]
			}
		});
		//glue.log("FINE", "covPrimerProbeMismatch.generateQueryToTargetRefSegs, alignResult", alignResult);
	});
	return alignResult.compoundAlignerResult.sequence[0].alignedSegment;
}


function reportMultiFasta(fastaFilePath) {
	throw new Error("reportMultiFasta is unimplemented");
}

function reportFastaDirectory(fastaDir) {
	throw new Error("reportFastaDirectory is unimplemented");
}