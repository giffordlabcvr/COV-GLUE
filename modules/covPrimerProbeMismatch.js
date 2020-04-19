// static map of concrete / ambiguous NT chars to the set of concrete chars which match them.
var ntCharToSubChars = {
    "A": ["A"],
	"C": ["C"],
	"G": ["G"],
	"T": ["T"],
	"R": ["A", "G"],
	"Y": ["C", "T"],
	"K": ["G", "T"],
	"M": ["A", "C"],
	"S": ["C", "G"],
	"W": ["A", "T"],
	"B": ["C", "G", "T"],
	"D": ["A", "G", "T"],
	"H": ["A", "C", "T"],
	"V": ["A", "C", "G"],
	"N": ["A", "C", "G", "T"]
}

// some mismatches that commonly occur, due to the primer/probe design.
// so shouldn't really be flagged up as a "red" problem.
var knownMismatchIssues = [
  { "ppId": "RdRP_SARSr-P1", "mismatch": "R15480C" },
  { "ppId": "RdRP_SARSr-P1", "mismatch": "T15489A" },
  { "ppId": "RdRP_SARSr-R1", "mismatch": "S15519T" },
  { "ppId": "NIID_2019-nCOV_N_R2", "mismatch": "C29277G" },
];

// build map concrete / ambiguous NT chars to a list of concrete / ambig chars which match them.
var ntCharToAllowed = {};

_.each(_.pairs(ntCharToSubChars), function(pair) {
	var ntChar = pair[0];
	var subChars = pair[1];
	var allowed;
	if(subChars.length == 1) {
		allowed = [ntChar];
	} else {
		var fullSubChars = subChars.slice(); // copy array
		_.each(_.pairs(ntCharToSubChars), function(pair2) {
			if(_.difference(pair2[1], subChars).length == 0 && fullSubChars.indexOf(pair2[0]) < 0) {
				fullSubChars.push(pair2[0]);
			}
		});
		allowed = fullSubChars;
	}
	ntCharToAllowed[ntChar] = allowed;
});

function reportSingleFasta(fastaFilePath) {
	
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
	
	var resultDoc = singleSequenceReportAux(fastaFilePath, queryName, queryNucleotides);
	glue.log("FINEST", "resultDoc", resultDoc);
	var htmlFilePath = fastaFilePath.substring(0, fastaFilePath.lastIndexOf("."))+"_ppReport.html";

	glue.inMode("module/covPrimerProbeReportTransformer", function() {
		glue.command({"transform-to-file" : {
			commandDocument: resultDoc,
			outputFile: htmlFilePath
		}});
	});
	glue.log("FINEST", "Wrote primer/probe report: "+htmlFilePath);
}


function reportSingleFastaFromDocument(document) {
	var fastaFilePath = document.inputDocument.fastaFilePath;
	var queryName = document.inputDocument.queryName;
	var queryNucleotides = document.inputDocument.queryNucleotides;
	var targetRefName = document.inputDocument.targetRefName;
	var queryToTargetRefSegs = document.inputDocument.queryToTargetRefSegs;
	return singleSequenceReport(fastaFilePath, queryName, queryNucleotides, targetRefName, queryToTargetRefSegs);
}

function singleSequenceReportAux(fastaFilePath, queryName, queryNucleotides) {
	var targetRefName = "REF_MASTER_WUHAN_HU_1";
	var queryToTargetRefSegs = generateQueryToTargetRefSegs(targetRefName, queryNucleotides);
	return singleSequenceReport(fastaFilePath, queryName, queryNucleotides, targetRefName, queryToTargetRefSegs);
}

function singleSequenceReport(fastaFilePath, queryName, queryNucleotides, targetRefName, queryToTargetRefSegs) {
	
	var projectVersion = 
		glue.command(["show","setting","project-version"]).projectShowSettingResult.settingValue;
	var glueVersion = 
		glue.command(["glue-engine","show-version"]).glueEngineShowVersionResult.glueEngineVersion;
	
	var rasVariationMatchDocument;
	
	var matchResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_match:%'");
	
	var vNameToMatchResult = {};
	_.each(matchResults, function(mr) {vNameToMatchResult[mr.variationName] = mr});
	
	var anyDelResult = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name = 'cov_pp_any_deletion'")[0];
	
	var assayObjs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_primer_probe_assay", 
				"-s", "organisation",
				//"-w", "id = 'NIID_2019-nCOV_N'", // Limited test 
				//"-w", "id = 'E_sarbeco'", // Limited test 
				"id", "display_name", "organisation", "url", "assay_type", "assay_type_display"]));
	
	_.each(assayObjs, function(assayObj) {
		glue.inMode("custom-table-row/cov_primer_probe_assay/"+assayObj.id, function() {
			assayObj.ppObjs = glue.tableToObjects(glue.command(["list", "link-target", "cov_primer_probe", 
				"id", "display_name", "sequence_fwd", "sequence_fwd_regex", "sequence_rev", "sequence_rev_regex", 
				"ref_hits", "sequence_to_scan", "fwd_orientation", "ref_start", "ref_end", "length", 
				"seq_match.name", "seq_insertion.name"]));
		});
		assayObj.primersWithIssues = 0;
		_.each(assayObj.ppObjs, function(ppObj) {
			ppObj.numIssues = 0;
			ppObj.numKnownIssues = 0;
			ppObj.issues = [];
			ppObj.displayIssues = [];
			ppObj.unknownDisplayIssues = [];
			// use the seq_match variation as a first pass. This will normally get
			// sufficientCoverage == true and present == true, in which case, no issues.
			// Otherwise run the insertion, deletion and single mismatch variations and
			// generate issues as necessary.
			ppObj.seqMatchResult = vNameToMatchResult[ppObj["seq_match.name"]];
			
			if( (!ppObj.seqMatchResult.sufficientCoverage) ||
				(!ppObj.seqMatchResult.present) ) {
				
				ppObj.seqInsertionResult = 
					variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, 
							"name = '"+ppObj["seq_insertion.name"]+"'")[0];
				ppObj.seqMismatchResults = 
					variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, 
							"cov_pp_seq_mismatch.id = '"+ppObj.id+"'");

				ppObj.seqMismatchResults = _.sortBy(ppObj.seqMismatchResults, function(smr) {return smr.variationName; });
			
				var deletedLocs = [];
				if(anyDelResult.present == true) {
					var deletions = [];
					_.each(anyDelResult.matches, function(delMatch) {
						if(delMatch.refFirstNtDeleted <= ppObj.ref_end &&
								delMatch.refLastNtDeleted >= ppObj.ref_start) {
							var refFirstNtDeleted = delMatch.refFirstNtDeleted;
							var refLastNtDeleted = delMatch.refLastNtDeleted;
							var deletedRefNts = delMatch.deletedRefNts;
							ppObj.issues.push({ type: "deletion", 
								refFirstNtDeleted: refFirstNtDeleted,
								refLastNtDeleted: refLastNtDeleted,
								qryLastNtBeforeDel: delMatch.qryLastNtBeforeDel,
								qryFirstNtAfterDel: delMatch.qryFirstNtAfterDel,
								deletedRefNts: deletedRefNts });
							var delString;
							if(refFirstNtDeleted == refLastNtDeleted) {
								delString = refFirstNtDeleted;
							} else {
								delString = refFirstNtDeleted+"-"+refLastNtDeleted;
							}
							deletions.push(delString);
							for(var loc = refFirstNtDeleted; loc <= refLastNtDeleted; loc++) {
								if(loc >= ppObj.ref_start && loc <= ppObj.ref_end) {
									deletedLocs.push(loc);
								}
							}
						}
					});
					if(deletions.length > 0) {
						var displayString = deletions.length + " " +
							(deletions.length == 1 ? "deletion" : "deletions") + ": " +
								deletions.join(", ");
						ppObj.displayIssues.push(displayString);
						ppObj.unknownDisplayIssues.push(displayString);

					}
				}

				
				var insufficientCoverageRegions = [];
				var currentInsufficientCoverage = null;
				var mismatches = [];
				var unknownMismatches = [];
				
				var sequenceToScan = ppObj.sequence_to_scan;
				for(var i = 0; i < sequenceToScan.length; i++) {
					var refCoord = ppObj.ref_start+i;
					if(deletedLocs.indexOf(refCoord) >= 0) {
						continue;
					}
					var mismatchResult = ppObj.seqMismatchResults[i];
					var sufficientCoverageLoc = true;
					if(mismatchResult.sufficientCoverage == false) {
						sufficientCoverageLoc = false;
					} else {
						if(mismatchResult.present == true) {
							var match = mismatchResult.matches[0];
							// offset accounts for flanking
							var offset = refCoord - match.refNtStart;
							var ppSequenceChar = sequenceToScan[i];
							var allowedNts = ntCharToAllowed[ppSequenceChar];
							var queryNt = match.queryNts[offset];
							var queryCoord = match.queryNtStart+offset;
							if(queryNt == "N") {
								sufficientCoverageLoc = false;
							} else {
								currentInsufficientCoverage = null;
								var mismatchString = ppSequenceChar+refCoord+queryNt;
								var knownIssue = _.find(knownMismatchIssues, function(kmi) {
									return kmi.ppId == ppObj.id && kmi.mismatch == mismatchString;
								}) != null;

								ppObj.issues.push({ type: "mismatch", 
									ppSequenceChar: ppSequenceChar,
									allowedNts: allowedNts,
									queryNt: queryNt,
									refCoord: refCoord,
									queryCoord: queryCoord, 
									knownIssue: knownIssue});
								
								if(!knownIssue) {
									unknownMismatches.push(mismatchString);
								}
								mismatches.push(mismatchString+(knownIssue?"*":""));
							}
						} else {
							currentInsufficientCoverage = null;
						}
					}
					if(sufficientCoverageLoc == false) {
						if(currentInsufficientCoverage == null) {
							currentInsufficientCoverage = {
									ref_start: refCoord,
									ref_end: refCoord
							};
							insufficientCoverageRegions.push(currentInsufficientCoverage);
						} else {
							currentInsufficientCoverage.ref_end = refCoord;
						}
					}
				}
				
				if(insufficientCoverageRegions.length > 0) {
					ppObj.issues.push({ type: "coverageAlmtIssues", 
						regions: insufficientCoverageRegions });
					var displayRegions = _.map(insufficientCoverageRegions, function(region) {
						if(region.ref_start == region.ref_end) { 
							return region.ref_start; 
						} else {
							return region.ref_start + "-" + region.ref_end;
						}
					});
					var message = "Coverage/alignment issues at "+displayRegions.join(", ");
					ppObj.displayIssues.push(message);
					ppObj.unknownDisplayIssues.push(message);
				}
				if(mismatches.length > 0) {
					var displayString = mismatches.length + " " +
						(mismatches.length == 1 ? "mismatch" : "mismatches") + ": " +
						mismatches.join(", ");
					ppObj.displayIssues.push(displayString);
				} 
				if(unknownMismatches.length > 0) {
					var displayString = unknownMismatches.length + " " +
						(unknownMismatches.length == 1 ? "mismatch" : "mismatches") + ": " +
						unknownMismatches.join(", ");
					ppObj.unknownDisplayIssues.push(displayString);
				} 
				if(ppObj.seqInsertionResult.sufficientCoverage == true && ppObj.seqInsertionResult.present == true) {
					var insertions = [];
					_.each(ppObj.seqInsertionResult.matches, function(insMatch) {
						var refLastNtBeforeIns = insMatch.refLastNtBeforeIns;
						var refFirstNtAfterIns = insMatch.refFirstNtAfterIns;
						var insertedQryNts = insMatch.insertedQryNts;
						ppObj.issues.push({ type: "insertion", 
							refLastNtBeforeIns: refLastNtBeforeIns,
							refFirstNtAfterIns: refFirstNtAfterIns,
							qryFirstInsertedNt: insMatch.qryFirstInsertedNt,
							qryLastInsertedNt: insMatch.qryLastInsertedNt,
							insertedQryNts: insertedQryNts });
						var displayInsertedNts = insertedQryNts;
						if(insertedQryNts.length > 6) {
							displayInsertedNts = insertedQryNts.substring(0, 3)+"..."+insertedQryNts.substring(insertedQryNts.length-3, insertedQryNts.length);
						}
						insertions.push(refLastNtBeforeIns+"-"+displayInsertedNts+"-"+refFirstNtAfterIns);
					});
					if(insertions.length > 0) {
						var displayString = insertions.length + " " +
							(insertions.length == 1 ? "insertion" : "insertions") + ": " +
								insertions.join(", ");
						ppObj.displayIssues.push(displayString);
						ppObj.unknownDisplayIssues.push(displayString);
					}
				}
				ppObj.issues = _.sortBy(ppObj.issues, function(iss) {
					if(iss.type == "insertion") {
						return iss.refLastNtBeforeIns;
					} else if(iss.type == "deletion") {
						return iss.refFirstNtDeleted;
					} else if(iss.type == "mismatch") {
						return iss.refCoord;
					} else if(iss.type == "coverageAlmtIssues") {
						return iss.regions[0].ref_start;
					} else {
						throw new Error("Unexpected else");
					}
				});
				if(ppObj.issues.length > 0) {
					ppObj.numIssues = ppObj.issues.length;
					ppObj.numKnownIssues = _.filter(ppObj.issues, function(iss) {return iss.knownIssue;} ).length;
					assayObj.primersWithIssues ++;
				}
			}

		});		
	});
	
	var resultDoc = {
		covPPReport: {
			fastaFilePath: fastaFilePath,
			sequenceID: queryName,
			assays: assayObjs,
			projectVersion: projectVersion,
			glueVersion: glueVersion,
			reportGenerationDate: todaysDate()
		}
	};
	
	return resultDoc;
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
	return alignResult.mafftAlignerResult.sequence[0].alignedSegment;
}

function runTests() {
	runTest("vanilla", [], []);
	runTest("E_sarbeco_F1, T26274G", [{type:"replaceSingle",loc:26274, replacement:"G"}], 
			[
				{
				    "type": "mismatch",
				    "ppSequenceChar": "T",
				    "allowedNts": [
				      "T"
				    ],
				    "queryNt": "G",
				    "refCoord": 26274,
				    "queryCoord": 26274,
				    "knownIssue": false,
				    "ppId": "E_Sarbeco_F1"
				}
			]);
	runTest("E_sarbeco_F1, 26274del", [{type:"replaceSingle",loc:26274, replacement:""}], 
			[
				{
				    "type": "deletion",
				    "refFirstNtDeleted": 26274,
				    "refLastNtDeleted": 26274,
				    "qryLastNtBeforeDel": 26273,
				    "qryFirstNtAfterDel": 26274,
				    "deletedRefNts": "T",
				    "ppId": "E_Sarbeco_F1"
				}
			]);
	runTest("E_sarbeco_F1, 26275-26277del", [{type:"replaceRegion",locStart:26275, locEnd:26277, replacement:""}], 
			[
				{
				    "type": "deletion",
				    "refFirstNtDeleted": 26275,
				    "refLastNtDeleted": 26277,
				    "qryLastNtBeforeDel": 26274,
				    "qryFirstNtAfterDel": 26275,
				    "deletedRefNts": "ACG",
				    "ppId": "E_Sarbeco_F1"
				  }
			]);
	// deletion of 4 nts at end of E_sarbeco_F1, this is also the first 4 of nCoV-2019_86_RIGHT
	runTest("E_sarbeco_F1, nCoV-2019_86_RIGHT, 26291-26294del", [{type:"replaceRegion",locStart:26291, locEnd:26294, replacement:""}], 
			[
				{
				    "type": "deletion",
				    "refFirstNtDeleted": 26291,
				    "refLastNtDeleted": 26294,
				    "qryLastNtBeforeDel": 26290,
				    "qryFirstNtAfterDel": 26291,
				    "deletedRefNts": "GCGT",
				    "ppId": "E_Sarbeco_F1"
				  },
				  {
				    "type": "deletion",
				    "refFirstNtDeleted": 26291,
				    "refLastNtDeleted": 26294,
				    "qryLastNtBeforeDel": 26290,
				    "qryFirstNtAfterDel": 26291,
				    "deletedRefNts": "GCGT",
				    "ppId": "nCoV-2019_86_RIGHT"
				  }
			]);
	// deletion of 13 nts overlapping the start of E_sarbeco_F1
	runTest("E_sarbeco_F1 26260-26272del", [{type:"replaceRegion",locStart:26260, locEnd:26272, replacement:""}], 
			[
				{
				    "type": "deletion",
				    "refFirstNtDeleted": 26260,
				    "refLastNtDeleted": 26272,
				    "qryLastNtBeforeDel": 26259,
				    "qryFirstNtAfterDel": 26260,
				    "deletedRefNts": "TCGGAAGAGACAG",
				    "ppId": "E_Sarbeco_F1"
				  }

			]);
	// deletion of 10 nts overlapping the end of E_sarbeco_F1 and start of nCoV-2019_86_RIGHT
	runTest("E_sarbeco_F1, nCoV-2019_86_RIGHT 26286-26295del", [{type:"replaceRegion",locStart:26286, locEnd:26295, replacement:""}], 
			[
				{
				    "type": "deletion",
				    "refFirstNtDeleted": 26286,
				    "refLastNtDeleted": 26295,
				    "qryLastNtBeforeDel": 26285,
				    "qryFirstNtAfterDel": 26286,
				    "deletedRefNts": "TAATAGCGTA",
				    "ppId": "E_Sarbeco_F1"
				  },
				  {
				    "type": "deletion",
				    "refFirstNtDeleted": 26286,
				    "refLastNtDeleted": 26295,
				    "qryLastNtBeforeDel": 26285,
				    "qryFirstNtAfterDel": 26286,
				    "deletedRefNts": "TAATAGCGTA",
				    "ppId": "nCoV-2019_86_RIGHT"
				  }
			]);
			runTest("E_sarbeco_F1, A26271G, A26275G, T26275G", [{type:"replaceSingle",loc:26271, replacement:"G"},
				{type:"replaceSingle",loc:26275, replacement:"G"},
				{type:"replaceSingle",loc:26278, replacement:"G"}], 
			[
				  {
					    "type": "mismatch",
					    "ppSequenceChar": "A",
					    "allowedNts": [
					      "A"
					    ],
					    "queryNt": "G",
					    "refCoord": 26271,
					    "queryCoord": 26271,
					    "knownIssue": false,
					    "ppId": "E_Sarbeco_F1"
					  },
					  {
					    "type": "mismatch",
					    "ppSequenceChar": "A",
					    "allowedNts": [
					      "A"
					    ],
					    "queryNt": "G",
					    "refCoord": 26275,
					    "queryCoord": 26275,
					    "knownIssue": false,
					    "ppId": "E_Sarbeco_F1"
					  },
					  {
					    "type": "mismatch",
					    "ppSequenceChar": "T",
					    "allowedNts": [
					      "T"
					    ],
					    "queryNt": "G",
					    "refCoord": 26278,
					    "queryCoord": 26278,
					    "knownIssue": false,
					    "ppId": "E_Sarbeco_F1"
					  }
			]);
			runTest("E_sarbeco_F1 insertion 26270-TAGAT-26271", [{type:"insertRegion",afterLoc:26270,insertion:"TAGAT"}], 
			[
				{
				    "type": "insertion",
				    "refLastNtBeforeIns": 26270,
				    "refFirstNtAfterIns": 26271,
				    "qryFirstInsertedNt": 26271,
				    "qryLastInsertedNt": 26275,
				    "insertedQryNts": "TAGAT",
				    "ppId": "E_Sarbeco_F1"
				  }
			]);
			runTest("E_sarbeco_F1 various", [
				{type:"insertRegion",afterLoc:26270,insertion:"TAGTAGATAT"},
				// note, locations shifted after insertion
				{type:"replaceRegion",locStart:26296, locEnd:26305, replacement:""},
				{type:"replaceSingle",loc:26288, replacement:"G"}, 
				{type:"replaceSingle",loc:26290, replacement:"N"}], 
			[
				{
				    "type": "insertion",
				    "refLastNtBeforeIns": 26270,
				    "refFirstNtAfterIns": 26271,
				    "qryFirstInsertedNt": 26271,
				    "qryLastInsertedNt": 26280,
				    "insertedQryNts": "TAGTAGATAT",
				    "ppId": "E_Sarbeco_F1"
				  },
				  {
				    "type": "mismatch",
				    "ppSequenceChar": "T",
				    "allowedNts": [
				      "T"
				    ],
				    "queryNt": "G",
				    "refCoord": 26278,
				    "queryCoord": 26288,
				    "knownIssue": false,
				    "ppId": "E_Sarbeco_F1"
				  },
				  {
				    "type": "coverageAlmtIssues",
				    "regions": [
				      {
				        "ref_start": 26280,
				        "ref_end": 26280
				      }
				    ],
				    "ppId": "E_Sarbeco_F1"
				  },
				  {
				    "type": "deletion",
				    "refFirstNtDeleted": 26286,
				    "refLastNtDeleted": 26295,
				    "qryLastNtBeforeDel": 26295,
				    "qryFirstNtAfterDel": 26296,
				    "deletedRefNts": "TAATAGCGTA",
				    "ppId": "E_Sarbeco_F1"
				  },
				  {
				    "type": "deletion",
				    "refFirstNtDeleted": 26286,
				    "refLastNtDeleted": 26295,
				    "qryLastNtBeforeDel": 26295,
				    "qryFirstNtAfterDel": 26296,
				    "deletedRefNts": "TAATAGCGTA",
				    "ppId": "nCoV-2019_86_RIGHT"
				  }
		]);
	// this primer location has a Y (C/T) ambiguity character.
	// the reference has a T, so replacing with a C or Y should give no issues
	runTest("2019-nCoV_N3-P ambig 1", [{type:"replaceSingle",loc:28705, replacement:"C"}], []);
	runTest("2019-nCoV_N3-P ambig 2", [{type:"replaceSingle",loc:28705, replacement:"T"}], []);
	runTest("2019-nCoV_N3-P ambig 3", [{type:"replaceSingle",loc:28705, replacement:"Y"}], []);
	// but replacing with a G should give a mismatch.
	runTest("2019-nCoV_N3-P ambig 4", [{type:"replaceSingle",loc:28705, replacement:"G"}], [
		{
		    "type": "mismatch",
		    "ppSequenceChar": "Y",
		    "allowedNts": [
		      "C",
		      "T",
		      "Y"
		    ],
		    "queryNt": "G",
		    "refCoord": 28705,
		    "queryCoord": 28705,
		    "knownIssue": false,
		    "ppId": "2019-nCoV_N3-P"
		  }
	]);
	runTest("nCoV-2019_86_RIGHT coverage 26306-26308", [{type:"replaceSingle",loc:26306, replacement:"N"},
		{type:"replaceSingle",loc:26307, replacement:"N"},
		{type:"replaceSingle",loc:26308, replacement:"N"}], 
			[
				{
				    "type": "coverageAlmtIssues",
				    "regions": [
				      {
				        "ref_start": 26306,
				        "ref_end": 26308
				      }
				    ],
				    "ppId": "nCoV-2019_86_RIGHT"
				  }
			]);

	
		
}

function runTest(testName, mods, expectedIssues) {
	var testSeq = buildTestSequence(mods);
	var report = singleSequenceReportAux(testName, testName, testSeq);
	var actualIssues = [];
	_.each(report.covPPReport.assays, function(assay) {
		_.each(assay.ppObjs, function(ppObj) {
			var issues = ppObj.issues;
			_.each(issues, function(issue) {
				issue.ppId = ppObj.id;
				if(!issue.knownIssue) {
					actualIssues.push(issue);
				}
			});
		});
	});
	actualIssues = _.sortBy(actualIssues, function(ai) {return ai.ppId} );
	if(_.isEqual(expectedIssues, actualIssues)) {
		glue.log("INFO", "Test '"+testName+"' passed");
	} else {
		glue.logInfo("report", report);
		glue.log("INFO", "Test '"+testName+"' failed");
		glue.log("FINEST", "Test '"+testName+"' expectedIssues", expectedIssues);
		glue.log("FINEST", "Test '"+testName+"' actualIssues", actualIssues);
		throw new Error("Test failed");
	}
	
} 
	
function buildTestSequence(mods) {
	var wuhan_hu_1 = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTATGTCGATAACAACTTCTGTGGCCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGCACTTTGTCCGAACAACTGGACTTTATTGACACTAAGAGGGGTGTATACTGCTGCCGTGAACATGAGCATGAAATTGCTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAACCAAGGGTTGAAAAGAAAAAGCTTGATGGCTTTATGGGTAGAATTCGATCTGTCTATCCAGTTGCGTCACCAAATGAATGCAACCAAATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCAGACGGGCGATTTTGTTAAAGCCACTTGCGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAATTTATTGTCCAGCATGTCACAATTCAGAAGTAGGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAAACCATTCTTCGTAAGGGTGGTCGCACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGCCTATTGGGTTCCACGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATGACAACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGAGATCGCCATTATTTTGGCATCTTTTTCTGCTTCCACAAGTGCTTTTGTGGAAACTGTGAAAGGTTTGGATTATAAAGCATTCAAACAAATTGTTGAATCCTGTGGTAATTTTAAAGTTACAAAAGGAAAAGCTAAAAAAGGTGCCTGGAATATTGGTGAACAGAAATCAATACTGAGTCCTCTTTATGCATTTGCATCAGAGGCTGCTCGTGTTGTACGATCAATTTTCTCCCGCACTCTTGAAACTGCTCAAAATTCTGTGCGTGTTTTACAGAAGGCCGCTATAACAATACTAGATGGAATTTCACAGTATTCACTGAGACTCATTGATGCTATGATGTTCACATCTGATTTGGCTACTAACAATCTAGTTGTAATGGCCTACATTACAGGTGGTGTTGTTCAGTTGACTTCGCAGTGGCTAACTAACATCTTTGGCACTGTTTATGAAAAACTCAAACCCGTCCTTGATTGGCTTGAAGAGAAGTTTAAGGAAGGTGTAGAGTTTCTTAGAGACGGTTGGGAAATTGTTAAATTTATCTCAACCTGTGCTTGTGAAATTGTCGGTGGACAAATTGTCACCTGTGCAAAGGAAATTAAGGAGAGTGTTCAGACATTCTTTAAGCTTGTAAATAAATTTTTGGCTTTGTGTGCTGACTCTATCATTATTGGTGGAGCTAAACTTAAAGCCTTGAATTTAGGTGAAACATTTGTCACGCACTCAAAGGGATTGTACAGAAAGTGTGTTAAATCCAGAGAAGAAACTGGCCTACTCATGCCTCTAAAAGCCCCAAAAGAAATTATCTTCTTAGAGGGAGAAACACTTCCCACAGAAGTGTTAACAGAGGAAGTTGTCTTGAAAACTGGTGATTTACAACCATTAGAACAACCTACTAGTGAAGCTGTTGAAGCTCCATTGGTTGGTACACCAGTTTGTATTAACGGGCTTATGTTGCTCGAAATCAAAGACACAGAAAAGTACTGTGCCCTTGCACCTAATATGATGGTAACAAACAATACCTTCACACTCAAAGGCGGTGCACCAACAAAGGTTACTTTTGGTGATGACACTGTGATAGAAGTGCAAGGTTACAAGAGTGTGAATATCACTTTTGAACTTGATGAAAGGATTGATAAAGTACTTAATGAGAAGTGCTCTGCCTATACAGTTGAACTCGGTACAGAAGTAAATGAGTTCGCCTGTGTTGTGGCAGATGCTGTCATAAAAACTTTGCAACCAGTATCTGAATTACTTACACCACTGGGCATTGATTTAGATGAGTGGAGTATGGCTACATACTACTTATTTGATGAGTCTGGTGAGTTTAAATTGGCTTCACATATGTATTGTTCTTTCTACCCTCCAGATGAGGATGAAGAAGAAGGTGATTGTGAAGAAGAAGAGTTTGAGCCATCAACTCAATATGAGTATGGTACTGAAGATGATTACCAAGGTAAACCTTTGGAATTTGGTGCCACTTCTGCTGCTCTTCAACCTGAAGAAGAGCAAGAAGAAGATTGGTTAGATGATGATAGTCAACAAACTGTTGGTCAACAAGACGGCAGTGAGGACAATCAGACAACTACTATTCAAACAATTGTTGAGGTTCAACCTCAATTAGAGATGGAACTTACACCAGTTGTTCAGACTATTGAAGTGAATAGTTTTAGTGGTTATTTAAAACTTACTGACAATGTATACATTAAAAATGCAGACATTGTGGAAGAAGCTAAAAAGGTAAAACCAACAGTGGTTGTTAATGCAGCCAATGTTTACCTTAAACATGGAGGAGGTGTTGCAGGAGCCTTAAATAAGGCTACTAACAATGCCATGCAAGTTGAATCTGATGATTACATAGCTACTAATGGACCACTTAAAGTGGGTGGTAGTTGTGTTTTAAGCGGACACAATCTTGCTAAACACTGTCTTCATGTTGTCGGCCCAAATGTTAACAAAGGTGAAGACATTCAACTTCTTAAGAGTGCTTATGAAAATTTTAATCAGCACGAAGTTCTACTTGCACCATTATTATCAGCTGGTATTTTTGGTGCTGACCCTATACATTCTTTAAGAGTTTGTGTAGATACTGTTCGCACAAATGTCTACTTAGCTGTCTTTGATAAAAATCTCTATGACAAACTTGTTTCAAGCTTTTTGGAAATGAAGAGTGAAAAGCAAGTTGAACAAAAGATCGCTGAGATTCCTAAAGAGGAAGTTAAGCCATTTATAACTGAAAGTAAACCTTCAGTTGAACAGAGAAAACAAGATGATAAGAAAATCAAAGCTTGTGTTGAAGAAGTTACAACAACTCTGGAAGAAACTAAGTTCCTCACAGAAAACTTGTTACTTTATATTGACATTAATGGCAATCTTCATCCAGATTCTGCCACTCTTGTTAGTGACATTGACATCACTTTCTTAAAGAAAGATGCTCCATATATAGTGGGTGATGTTGTTCAAGAGGGTGTTTTAACTGCTGTGGTTATACCTACTAAAAAGGCTGGTGGCACTACTGAAATGCTAGCGAAAGCTTTGAGAAAAGTGCCAACAGACAATTATATAACCACTTACCCGGGTCAGGGTTTAAATGGTTACACTGTAGAGGAGGCAAAGACAGTGCTTAAAAAGTGTAAAAGTGCCTTTTACATTCTACCATCTATTATCTCTAATGAGAAGCAAGAAATTCTTGGAACTGTTTCTTGGAATTTGCGAGAAATGCTTGCACATGCAGAAGAAACACGCAAATTAATGCCTGTCTGTGTGGAAACTAAAGCCATAGTTTCAACTATACAGCGTAAATATAAGGGTATTAAAATACAAGAGGGTGTGGTTGATTATGGTGCTAGATTTTACTTTTACACCAGTAAAACAACTGTAGCGTCACTTATCAACACACTTAACGATCTAAATGAAACTCTTGTTACAATGCCACTTGGCTATGTAACACATGGCTTAAATTTGGAAGAAGCTGCTCGGTATATGAGATCTCTCAAAGTGCCAGCTACAGTTTCTGTTTCTTCACCTGATGCTGTTACAGCGTATAATGGTTATCTTACTTCTTCTTCTAAAACACCTGAAGAACATTTTATTGAAACCATCTCACTTGCTGGTTCCTATAAAGATTGGTCCTATTCTGGACAATCTACACAACTAGGTATAGAATTTCTTAAGAGAGGTGATAAAAGTGTATATTACACTAGTAATCCTACCACATTCCACCTAGATGGTGAAGTTATCACCTTTGACAATCTTAAGACACTTCTTTCTTTGAGAGAAGTGAGGACTATTAAGGTGTTTACAACAGTAGACAACATTAACCTCCACACGCAAGTTGTGGACATGTCAATGACATATGGACAACAGTTTGGTCCAACTTATTTGGATGGAGCTGATGTTACTAAAATAAAACCTCATAATTCACATGAAGGTAAAACATTTTATGTTTTACCTAATGATGACACTCTACGTGTTGAGGCTTTTGAGTACTACCACACAACTGATCCTAGTTTTCTGGGTAGGTACATGTCAGCATTAAATCACACTAAAAAGTGGAAATACCCACAAGTTAATGGTTTAACTTCTATTAAATGGGCAGATAACAACTGTTATCTTGCCACTGCATTGTTAACACTCCAACAAATAGAGTTGAAGTTTAATCCACCTGCTCTACAAGATGCTTATTACAGAGCAAGGGCTGGTGAAGCTGCTAACTTTTGTGCACTTATCTTAGCCTACTGTAATAAGACAGTAGGTGAGTTAGGTGATGTTAGAGAAACAATGAGTTACTTGTTTCAACATGCCAATTTAGATTCTTGCAAAAGAGTCTTGAACGTGGTGTGTAAAACTTGTGGACAACAGCAGACAACCCTTAAGGGTGTAGAAGCTGTTATGTACATGGGCACACTTTCTTATGAACAATTTAAGAAAGGTGTTCAGATACCTTGTACGTGTGGTAAACAAGCTACAAAATATCTAGTACAACAGGAGTCACCTTTTGTTATGATGTCAGCACCACCTGCTCAGTATGAACTTAAGCATGGTACATTTACTTGTGCTAGTGAGTACACTGGTAATTACCAGTGTGGTCACTATAAACATATAACTTCTAAAGAAACTTTGTATTGCATAGACGGTGCTTTACTTACAAAGTCCTCAGAATACAAAGGTCCTATTACGGATGTTTTCTACAAAGAAAACAGTTACACAACAACCATAAAACCAGTTACTTATAAATTGGATGGTGTTGTTTGTACAGAAATTGACCCTAAGTTGGACAATTATTATAAGAAAGACAATTCTTATTTCACAGAGCAACCAATTGATCTTGTACCAAACCAACCATATCCAAACGCAAGCTTCGATAATTTTAAGTTTGTATGTGATAATATCAAATTTGCTGATGATTTAAACCAGTTAACTGGTTATAAGAAACCTGCTTCAAGAGAGCTTAAAGTTACATTTTTCCCTGACTTAAATGGTGATGTGGTGGCTATTGATTATAAACACTACACACCCTCTTTTAAGAAAGGAGCTAAATTGTTACATAAACCTATTGTTTGGCATGTTAACAATGCAACTAATAAAGCCACGTATAAACCAAATACCTGGTGTATACGTTGTCTTTGGAGCACAAAACCAGTTGAAACATCAAATTCGTTTGATGTACTGAAGTCAGAGGACGCGCAGGGAATGGATAATCTTGCCTGCGAAGATCTAAAACCAGTCTCTGAAGAAGTAGTGGAAAATCCTACCATACAGAAAGACGTTCTTGAGTGTAATGTGAAAACTACCGAAGTTGTAGGAGACATTATACTTAAACCAGCAAATAATAGTTTAAAAATTACAGAAGAGGTTGGCCACACAGATCTAATGGCTGCTTATGTAGACAATTCTAGTCTTACTATTAAGAAACCTAATGAATTATCTAGAGTATTAGGTTTGAAAACCCTTGCTACTCATGGTTTAGCTGCTGTTAATAGTGTCCCTTGGGATACTATAGCTAATTATGCTAAGCCTTTTCTTAACAAAGTTGTTAGTACAACTACTAACATAGTTACACGGTGTTTAAACCGTGTTTGTACTAATTATATGCCTTATTTCTTTACTTTATTGCTACAATTGTGTACTTTTACTAGAAGTACAAATTCTAGAATTAAAGCATCTATGCCGACTACTATAGCAAAGAATACTGTTAAGAGTGTCGGTAAATTTTGTCTAGAGGCTTCATTTAATTATTTGAAGTCACCTAATTTTTCTAAACTGATAAATATTATAATTTGGTTTTTACTATTAAGTGTTTGCCTAGGTTCTTTAATCTACTCAACCGCTGCTTTAGGTGTTTTAATGTCTAATTTAGGCATGCCTTCTTACTGTACTGGTTACAGAGAAGGCTATTTGAACTCTACTAATGTCACTATTGCAACCTACTGTACTGGTTCTATACCTTGTAGTGTTTGTCTTAGTGGTTTAGATTCTTTAGACACCTATCCTTCTTTAGAAACTATACAAATTACCATTTCATCTTTTAAATGGGATTTAACTGCTTTTGGCTTAGTTGCAGAGTGGTTTTTGGCATATATTCTTTTCACTAGGTTTTTCTATGTACTTGGATTGGCTGCAATCATGCAATTGTTTTTCAGCTATTTTGCAGTACATTTTATTAGTAATTCTTGGCTTATGTGGTTAATAATTAATCTTGTACAAATGGCCCCGATTTCAGCTATGGTTAGAATGTACATCTTCTTTGCATCATTTTATTATGTATGGAAAAGTTATGTGCATGTTGTAGACGGTTGTAATTCATCAACTTGTATGATGTGTTACAAACGTAATAGAGCAACAAGAGTCGAATGTACAACTATTGTTAATGGTGTTAGAAGGTCCTTTTATGTCTATGCTAATGGAGGTAAAGGCTTTTGCAAACTACACAATTGGAATTGTGTTAATTGTGATACATTCTGTGCTGGTAGTACATTTATTAGTGATGAAGTTGCGAGAGACTTGTCACTACAGTTTAAAAGACCAATAAATCCTACTGACCAGTCTTCTTACATCGTTGATAGTGTTACAGTGAAGAATGGTTCCATCCATCTTTACTTTGATAAAGCTGGTCAAAAGACTTATGAAAGACATTCTCTCTCTCATTTTGTTAACTTAGACAACCTGAGAGCTAATAACACTAAAGGTTCATTGCCTATTAATGTTATAGTTTTTGATGGTAAATCAAAATGTGAAGAATCATCTGCAAAATCAGCGTCTGTTTACTACAGTCAGCTTATGTGTCAACCTATACTGTTACTAGATCAGGCATTAGTGTCTGATGTTGGTGATAGTGCGGAAGTTGCAGTTAAAATGTTTGATGCTTACGTTAATACGTTTTCATCAACTTTTAACGTACCAATGGAAAAACTCAAAACACTAGTTGCAACTGCAGAAGCTGAACTTGCAAAGAATGTGTCCTTAGACAATGTCTTATCTACTTTTATTTCAGCAGCTCGGCAAGGGTTTGTTGATTCAGATGTAGAAACTAAAGATGTTGTTGAATGTCTTAAATTGTCACATCAATCTGACATAGAAGTTACTGGCGATAGTTGTAATAACTATATGCTCACCTATAACAAAGTTGAAAACATGACACCCCGTGACCTTGGTGCTTGTATTGACTGTAGTGCGCGTCATATTAATGCGCAGGTAGCAAAAAGTCACAACATTGCTTTGATATGGAACGTTAAAGATTTCATGTCATTGTCTGAACAACTACGAAAACAAATACGTAGTGCTGCTAAAAAGAATAACTTACCTTTTAAGTTGACATGTGCAACTACTAGACAAGTTGTTAATGTTGTAACAACAAAGATAGCACTTAAGGGTGGTAAAATTGTTAATAATTGGTTGAAGCAGTTAATTAAAGTTACACTTGTGTTCCTTTTTGTTGCTGCTATTTTCTATTTAATAACACCTGTTCATGTCATGTCTAAACATACTGACTTTTCAAGTGAAATCATAGGATACAAGGCTATTGATGGTGGTGTCACTCGTGACATAGCATCTACAGATACTTGTTTTGCTAACAAACATGCTGATTTTGACACATGGTTTAGCCAGCGTGGTGGTAGTTATACTAATGACAAAGCTTGCCCATTGATTGCTGCAGTCATAACAAGAGAAGTGGGTTTTGTCGTGCCTGGTTTGCCTGGCACGATATTACGCACAACTAATGGTGACTTTTTGCATTTCTTACCTAGAGTTTTTAGTGCAGTTGGTAACATCTGTTACACACCATCAAAACTTATAGAGTACACTGACTTTGCAACATCAGCTTGTGTTTTGGCTGCTGAATGTACAATTTTTAAAGATGCTTCTGGTAAGCCAGTACCATATTGTTATGATACCAATGTACTAGAAGGTTCTGTTGCTTATGAAAGTTTACGCCCTGACACACGTTATGTGCTCATGGATGGCTCTATTATTCAATTTCCTAACACCTACCTTGAAGGTTCTGTTAGAGTGGTAACAACTTTTGATTCTGAGTACTGTAGGCACGGCACTTGTGAAAGATCAGAAGCTGGTGTTTGTGTATCTACTAGTGGTAGATGGGTACTTAACAATGATTATTACAGATCTTTACCAGGAGTTTTCTGTGGTGTAGATGCTGTAAATTTACTTACTAATATGTTTACACCACTAATTCAACCTATTGGTGCTTTGGACATATCAGCATCTATAGTAGCTGGTGGTATTGTAGCTATCGTAGTAACATGCCTTGCCTACTATTTTATGAGGTTTAGAAGAGCTTTTGGTGAATACAGTCATGTAGTTGCCTTTAATACTTTACTATTCCTTATGTCATTCACTGTACTCTGTTTAACACCAGTTTACTCATTCTTACCTGGTGTTTATTCTGTTATTTACTTGTACTTGACATTTTATCTTACTAATGATGTTTCTTTTTTAGCACATATTCAGTGGATGGTTATGTTCACACCTTTAGTACCTTTCTGGATAACAATTGCTTATATCATTTGTATTTCCACAAAGCATTTCTATTGGTTCTTTAGTAATTACCTAAAGAGACGTGTAGTCTTTAATGGTGTTTCCTTTAGTACTTTTGAAGAAGCTGCGCTGTGCACCTTTTTGTTAAATAAAGAAATGTATCTAAAGTTGCGTAGTGATGTGCTATTACCTCTTACGCAATATAATAGATACTTAGCTCTTTATAATAAGTACAAGTATTTTAGTGGAGCAATGGATACAACTAGCTACAGAGAAGCTGCTTGTTGTCATCTCGCAAAGGCTCTCAATGACTTCAGTAACTCAGGTTCTGATGTTCTTTACCAACCACCACAAACCTCTATCACCTCAGCTGTTTTGCAGAGTGGTTTTAGAAAAATGGCATTCCCATCTGGTAAAGTTGAGGGTTGTATGGTACAAGTAACTTGTGGTACAACTACACTTAACGGTCTTTGGCTTGATGACGTAGTTTACTGTCCAAGACATGTGATCTGCACCTCTGAAGACATGCTTAACCCTAATTATGAAGATTTACTCATTCGTAAGTCTAATCATAATTTCTTGGTACAGGCTGGTAATGTTCAACTCAGGGTTATTGGACATTCTATGCAAAATTGTGTACTTAAGCTTAAGGTTGATACAGCCAATCCTAAGACACCTAAGTATAAGTTTGTTCGCATTCAACCAGGACAGACTTTTTCAGTGTTAGCTTGTTACAATGGTTCACCATCTGGTGTTTACCAATGTGCTATGAGGCCCAATTTCACTATTAAGGGTTCATTCCTTAATGGTTCATGTGGTAGTGTTGGTTTTAACATAGATTATGACTGTGTCTCTTTTTGTTACATGCACCATATGGAATTACCAACTGGAGTTCATGCTGGCACAGACTTAGAAGGTAACTTTTATGGACCTTTTGTTGACAGGCAAACAGCACAAGCAGCTGGTACGGACACAACTATTACAGTTAATGTTTTAGCTTGGTTGTACGCTGCTGTTATAAATGGAGACAGGTGGTTTCTCAATCGATTTACCACAACTCTTAATGACTTTAACCTTGTGGCTATGAAGTACAATTATGAACCTCTAACACAAGACCATGTTGACATACTAGGACCTCTTTCTGCTCAAACTGGAATTGCCGTTTTAGATATGTGTGCTTCATTAAAAGAATTACTGCAAAATGGTATGAATGGACGTACCATATTGGGTAGTGCTTTATTAGAAGATGAATTTACACCTTTTGATGTTGTTAGACAATGCTCAGGTGTTACTTTCCAAAGTGCAGTGAAAAGAACAATCAAGGGTACACACCACTGGTTGTTACTCACAATTTTGACTTCACTTTTAGTTTTAGTCCAGAGTACTCAATGGTCTTTGTTCTTTTTTTTGTATGAAAATGCCTTTTTACCTTTTGCTATGGGTATTATTGCTATGTCTGCTTTTGCAATGATGTTTGTCAAACATAAGCATGCATTTCTCTGTTTGTTTTTGTTACCTTCTCTTGCCACTGTAGCTTATTTTAATATGGTCTATATGCCTGCTAGTTGGGTGATGCGTATTATGACATGGTTGGATATGGTTGATACTAGTTTGTCTGGTTTTAAGCTAAAAGACTGTGTTATGTATGCATCAGCTGTAGTGTTACTAATCCTTATGACAGCAAGAACTGTGTATGATGATGGTGCTAGGAGAGTGTGGACACTTATGAATGTCTTGACACTCGTTTATAAAGTTTATTATGGTAATGCTTTAGATCAAGCCATTTCCATGTGGGCTCTTATAATCTCTGTTACTTCTAACTACTCAGGTGTAGTTACAACTGTCATGTTTTTGGCCAGAGGTATTGTTTTTATGTGTGTTGAGTATTGCCCTATTTTCTTCATAACTGGTAATACACTTCAGTGTATAATGCTAGTTTATTGTTTCTTAGGCTATTTTTGTACTTGTTACTTTGGCCTCTTTTGTTTACTCAACCGCTACTTTAGACTGACTCTTGGTGTTTATGATTACTTAGTTTCTACACAGGAGTTTAGATATATGAATTCACAGGGACTACTCCCACCCAAGAATAGCATAGATGCCTTCAAACTCAACATTAAATTGTTGGGTGTTGGTGGCAAACCTTGTATCAAAGTAGCCACTGTACAGTCTAAAATGTCAGATGTAAAGTGCACATCAGTAGTCTTACTCTCAGTTTTGCAACAACTCAGAGTAGAATCATCATCTAAATTGTGGGCTCAATGTGTCCAGTTACACAATGACATTCTCTTAGCTAAAGATACTACTGAAGCCTTTGAAAAAATGGTTTCACTACTTTCTGTTTTGCTTTCCATGCAGGGTGCTGTAGACATAAACAAGCTTTGTGAAGAAATGCTGGACAACAGGGCAACCTTACAAGCTATAGCCTCAGAGTTTAGTTCCCTTCCATCATATGCAGCTTTTGCTACTGCTCAAGAAGCTTATGAGCAGGCTGTTGCTAATGGTGATTCTGAAGTTGTTCTTAAAAAGTTGAAGAAGTCTTTGAATGTGGCTAAATCTGAATTTGACCGTGATGCAGCCATGCAACGTAAGTTGGAAAAGATGGCTGATCAAGCTATGACCCAAATGTATAAACAGGCTAGATCTGAGGACAAGAGGGCAAAAGTTACTAGTGCTATGCAGACAATGCTTTTCACTATGCTTAGAAAGTTGGATAATGATGCACTCAACAACATTATCAACAATGCAAGAGATGGTTGTGTTCCCTTGAACATAATACCTCTTACAACAGCAGCCAAACTAATGGTTGTCATACCAGACTATAACACATATAAAAATACGTGTGATGGTACAACATTTACTTATGCATCAGCATTGTGGGAAATCCAACAGGTTGTAGATGCAGATAGTAAAATTGTTCAACTTAGTGAAATTAGTATGGACAATTCACCTAATTTAGCATGGCCTCTTATTGTAACAGCTTTAAGGGCCAATTCTGCTGTCAAATTACAGAATAATGAGCTTAGTCCTGTTGCACTACGACAGATGTCTTGTGCTGCCGGTACTACACAAACTGCTTGCACTGATGACAATGCGTTAGCTTACTACAACACAACAAAGGGAGGTAGGTTTGTACTTGCACTGTTATCCGATTTACAGGATTTGAAATGGGCTAGATTCCCTAAGAGTGATGGAACTGGTACTATCTATACAGAACTGGAACCACCTTGTAGGTTTGTTACAGACACACCTAAAGGTCCTAAAGTGAAGTATTTATACTTTATTAAAGGATTAAACAACCTAAATAGAGGTATGGTACTTGGTAGTTTAGCTGCCACAGTACGTCTACAAGCTGGTAATGCAACAGAAGTGCCTGCCAATTCAACTGTATTATCTTTCTGTGCTTTTGCTGTAGATGCTGCTAAAGCTTACAAAGATTATCTAGCTAGTGGGGGACAACCAATCACTAATTGTGTTAAGATGTTGTGTACACACACTGGTACTGGTCAGGCAATAACAGTTACACCGGAAGCCAATATGGATCAAGAATCCTTTGGTGGTGCATCGTGTTGTCTGTACTGCCGTTGCCACATAGATCATCCAAATCCTAAAGGATTTTGTGACTTAAAAGGTAAGTATGTACAAATACCTACAACTTGTGCTAATGACCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCGTGCGGCACAGGCACTAGTACTGATGTCGTATACAGGGCTTTTGACATCTACAATGATAAAGTAGCTGGTTTTGCTAAATTCCTAAAAACTAATTGTTGTCGCTTCCAAGAAAAGGACGAAGATGACAATTTAATTGATTCTTACTTTGTAGTTAAGAGACACACTTTCTCTAACTACCAACATGAAGAAACAATTTATAATTTACTTAAGGATTGTCCAGCTGTTGCTAAACATGACTTCTTTAAGTTTAGAATAGACGGTGACATGGTACCACATATATCACGTCAACGTCTTACTAAATACACAATGGCAGACCTCGTCTATGCTTTAAGGCATTTTGATGAAGGTAATTGTGACACATTAAAAGAAATACTTGTCACATACAATTGTTGTGATGATGATTATTTCAATAAAAAGGACTGGTATGATTTTGTAGAAAACCCAGATATATTACGCGTATACGCCAACTTAGGTGAACGTGTACGCCAAGCTTTGTTAAAAACAGTACAATTCTGTGATGCCATGCGAAATGCTGGTATTGTTGGTGTACTGACATTAGATAATCAAGATCTCAATGGTAACTGGTATGATTTCGGTGATTTCATACAAACCACGCCAGGTAGTGGAGTTCCTGTTGTAGATTCTTATTATTCATTGTTAATGCCTATATTAACCTTGACCAGGGCTTTAACTGCAGAGTCACATGTTGACACTGACTTAACAAAGCCTTACATTAAGTGGGATTTGTTAAAATATGACTTCACGGAAGAGAGGTTAAAACTCTTTGACCGTTATTTTAAATATTGGGATCAGACATACCACCCAAATTGTGTTAACTGTTTGGATGACAGATGCATTCTGCATTGTGCAAACTTTAATGTTTTATTCTCTACAGTGTTCCCACCTACAAGTTTTGGACCACTAGTGAGAAAAATATTTGTTGATGGTGTTCCATTTGTAGTTTCAACTGGATACCACTTCAGAGAGCTAGGTGTTGTACATAATCAGGATGTAAACTTACATAGCTCTAGACTTAGTTTTAAGGAATTACTTGTGTATGCTGCTGACCCTGCTATGCACGCTGCTTCTGGTAATCTATTACTAGATAAACGCACTACGTGCTTTTCAGTAGCTGCACTTACTAACAATGTTGCTTTTCAAACTGTCAAACCCGGTAATTTTAACAAAGACTTCTATGACTTTGCTGTGTCTAAGGGTTTCTTTAAGGAAGGAAGTTCTGTTGAATTAAAACACTTCTTCTTTGCTCAGGATGGTAATGCTGCTATCAGCGATTATGACTACTATCGTTATAATCTACCAACAATGTGTGATATCAGACAACTACTATTTGTAGTTGAAGTTGTTGATAAGTACTTTGATTGTTACGATGGTGGCTGTATTAATGCTAACCAAGTCATCGTCAACAACCTAGACAAATCAGCTGGTTTTCCATTTAATAAATGGGGTAAGGCTAGACTTTATTATGATTCAATGAGTTATGAGGATCAAGATGCACTTTTCGCATATACAAAACGTAATGTCATCCCTACTATAACTCAAATGAATCTTAAGTATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTGGCACAACATGTTAAAAACTGTTTATAGTGATGTAGAAAACCCTCACCTTATGGGTTGGGATTATCCTAAATGTGATAGAGCCATGCCTAACATGCTTAGAATTATGGCCTCACTTGTTCTTGCTCGCAAACATACAACGTGTTGTAGCTTGTCACACCGTTTCTATAGATTAGCTAATGAGTGTGCTCAAGTATTGAGTGAAATGGTCATGTGTGGCGGTTCACTATATGTTAAACCAGGTGGAACCTCATCAGGAGATGCCACAACTGCTTATGCTAATAGTGTTTTTAACATTTGTCAAGCTGTCACGGCCAATGTTAATGCACTTTTATCTACTGATGGTAACAAAATTGCCGATAAGTATGTCCGCAATTTACAACACAGACTTTATGAGTGTCTCTATAGAAATAGAGATGTTGACACAGACTTTGTGAATGAGTTTTACGCATATTTGCGTAAACATTTCTCAATGATGATACTCTCTGACGATGCTGTTGTGTGTTTCAATAGCACTTATGCATCTCAAGGTCTAGTGGCTAGCATAAAGAACTTTAAGTCAGTTCTTTATTATCAAAACAATGTTTTTATGTCTGAAGCAAAATGTTGGACTGAGACTGACCTTACTAAAGGACCTCATGAATTTTGCTCTCAACATACAATGCTAGTTAAACAGGGTGATGATTATGTGTACCTTCCTTACCCAGATCCATCAAGAATCCTAGGGGCCGGCTGTTTTGTAGATGATATCGTAAAAACAGATGGTACACTTATGATTGAACGGTTCGTGTCTTTAGCTATAGATGCTTACCCACTTACTAAACATCCTAATCAGGAGTATGCTGATGTCTTTCATTTGTACTTACAATACATAAGAAAGCTACATGATGAGTTAACAGGACACATGTTAGACATGTATTCTGTTATGCTTACTAATGATAACACTTCAAGGTATTGGGAACCTGAGTTTTATGAGGCTATGTACACACCGCATACAGTCTTACAGGCTGTTGGGGCTTGTGTTCTTTGCAATTCACAGACTTCATTAAGATGTGGTGCTTGCATACGTAGACCATTCTTATGTTGTAAATGCTGTTACGACCATGTCATATCAACATCACATAAATTAGTCTTGTCTGTTAATCCGTATGTTTGCAATGCTCCAGGTTGTGATGTCACAGATGTGACTCAACTTTACTTAGGAGGTATGAGCTATTATTGTAAATCACATAAACCACCCATTAGTTTTCCATTGTGTGCTAATGGACAAGTTTTTGGTTTATATAAAAATACATGTGTTGGTAGCGATAATGTTACTGACTTTAATGCAATTGCAACATGTGACTGGACAAATGCTGGTGATTACATTTTAGCTAACACCTGTACTGAAAGACTCAAGCTTTTTGCAGCAGAAACGCTCAAAGCTACTGAGGAGACATTTAAACTGTCTTATGGTATTGCTACTGTACGTGAAGTGCTGTCTGACAGAGAATTACATCTTTCATGGGAAGTTGGTAAACCTAGACCACCACTTAACCGAAATTATGTCTTTACTGGTTATCGTGTAACTAAAAACAGTAAAGTACAAATAGGAGAGTACACCTTTGAAAAAGGTGACTATGGTGATGCTGTTGTTTACCGAGGTACAACAACTTACAAATTAAATGTTGGTGATTATTTTGTGCTGACATCACATACAGTAATGCCATTAAGTGCACCTACACTAGTGCCACAAGAGCACTATGTTAGAATTACTGGCTTATACCCAACACTCAATATCTCAGATGAGTTTTCTAGCAATGTTGCAAATTATCAAAAGGTTGGTATGCAAAAGTATTCTACACTCCAGGGACCACCTGGTACTGGTAAGAGTCATTTTGCTATTGGCCTAGCTCTCTACTACCCTTCTGCTCGCATAGTGTATACAGCTTGCTCTCATGCCGCTGTTGATGCACTATGTGAGAAGGCATTAAAATATTTGCCTATAGATAAATGTAGTAGAATTATACCTGCACGTGCTCGTGTAGAGTGTTTTGATAAATTCAAAGTGAATTCAACATTAGAACAGTATGTCTTTTGTACTGTAAATGCATTGCCTGAGACGACAGCAGATATAGTTGTCTTTGATGAAATTTCAATGGCCACAAATTATGATTTGAGTGTTGTCAATGCCAGATTACGTGCTAAGCACTATGTGTACATTGGCGACCCTGCTCAATTACCTGCACCACGCACATTGCTAACTAAGGGCACACTAGAACCAGAATATTTCAATTCAGTGTGTAGACTTATGAAAACTATAGGTCCAGACATGTTCCTCGGAACTTGTCGGCGTTGTCCTGCTGAAATTGTTGACACTGTGAGTGCTTTGGTTTATGATAATAAGCTTAAAGCACATAAAGACAAATCAGCTCAATGCTTTAAAATGTTTTATAAGGGTGTTATCACGCATGATGTTTCATCTGCAATTAACAGGCCACAAATAGGCGTGGTAAGAGAATTCCTTACACGTAACCCTGCTTGGAGAAAAGCTGTCTTTATTTCACCTTATAATTCACAGAATGCTGTAGCCTCAAAGATTTTGGGACTACCAACTCAAACTGTTGATTCATCACAGGGCTCAGAATATGACTATGTCATATTCACTCAAACCACTGAAACAGCTCACTCTTGTAATGTAAACAGATTTAATGTTGCTATTACCAGAGCAAAAGTAGGCATACTTTGCATAATGTCTGATAGAGACCTTTATGACAAGTTGCAATTTACAAGTCTTGAAATTCCACGTAGGAATGTGGCAACTTTACAAGCTGAAAATGTAACAGGACTCTTTAAAGATTGTAGTAAGGTAATCACTGGGTTACATCCTACACAGGCACCTACACACCTCAGTGTTGACACTAAATTCAAAACTGAAGGTTTATGTGTTGACATACCTGGCATACCTAAGGACATGACCTATAGAAGACTCATCTCTATGATGGGTTTTAAAATGAATTATCAAGTTAATGGTTACCCTAACATGTTTATCACCCGCGAAGAAGCTATAAGACATGTACGTGCATGGATTGGCTTCGATGTCGAGGGGTGTCATGCTACTAGAGAAGCTGTTGGTACCAATTTACCTTTACAGCTAGGTTTTTCTACAGGTGTTAACCTAGTTGCTGTACCTACAGGTTATGTTGATACACCTAATAATACAGATTTTTCCAGAGTTAGTGCTAAACCACCGCCTGGAGATCAATTTAAACACCTCATACCACTTATGTACAAAGGACTTCCTTGGAATGTAGTGCGTATAAAGATTGTACAAATGTTAAGTGACACACTTAAAAATCTCTCTGACAGAGTCGTATTTGTCTTATGGGCACATGGCTTTGAGTTGACATCTATGAAGTATTTTGTGAAAATAGGACCTGAGCGCACCTGTTGTCTATGTGATAGACGTGCCACATGCTTTTCCACTGCTTCAGACACTTATGCCTGTTGGCATCATTCTATTGGATTTGATTACGTCTATAATCCGTTTATGATTGATGTTCAACAATGGGGTTTTACAGGTAACCTACAAAGCAACCATGATCTGTATTGTCAAGTCCATGGTAATGCACATGTAGCTAGTTGTGATGCAATCATGACTAGGTGTCTAGCTGTCCACGAGTGCTTTGTTAAGCGTGTTGACTGGACTATTGAATATCCTATAATTGGTGATGAACTGAAGATTAATGCGGCTTGTAGAAAGGTTCAACACATGGTTGTTAAAGCTGCATTATTAGCAGACAAATTCCCAGTTCTTCACGACATTGGTAACCCTAAAGCTATTAAGTGTGTACCTCAAGCTGATGTAGAATGGAAGTTCTATGATGCACAGCCTTGTAGTGACAAAGCTTATAAAATAGAAGAATTATTCTATTCTTATGCCACACATTCTGACAAATTCACAGATGGTGTATGCCTATTTTGGAATTGCAATGTCGATAGATATCCTGCTAATTCCATTGTTTGTAGATTTGACACTAGAGTGCTATCTAACCTTAACTTGCCTGGTTGTGATGGTGGCAGTTTGTATGTAAATAAACATGCATTCCACACACCAGCTTTTGATAAAAGTGCTTTTGTTAATTTAAAACAATTACCATTTTTCTATTACTCTGACAGTCCATGTGAGTCTCATGGAAAACAAGTAGTGTCAGATATAGATTATGTACCACTAAAGTCTGCTACGTGTATAACACGTTGCAATTTAGGTGGTGCTGTCTGTAGACATCATGCTAATGAGTACAGATTGTATCTCGATGCTTATAACATGATGATCTCAGCTGGCTTTAGCTTGTGGGTTTACAAACAATTTGATACTTATAACCTCTGGAACACTTTTACAAGACTTCAGAGTTTAGAAAATGTGGCTTTTAATGTTGTAAATAAGGGACACTTTGATGGACAACAGGGTGAAGTACCAGTTTCTATCATTAATAACACTGTTTACACAAAAGTTGATGGTGTTGATGTAGAATTGTTTGAAAATAAAACAACATTACCTGTTAATGTAGCATTTGAGCTTTGGGCTAAGCGCAACATTAAACCAGTACCAGAGGTGAAAATACTCAATAATTTGGGTGTGGACATTGCTGCTAATACTGTGATCTGGGACTACAAAAGAGATGCTCCAGCACATATATCTACTATTGGTGTTTGTTCTATGACTGACATAGCCAAGAAACCAACTGAAACGATTTGTGCACCACTCACTGTCTTTTTTGATGGTAGAGTTGATGGTCAAGTAGACTTATTTAGAAATGCCCGTAATGGTGTTCTTATTACAGAAGGTAGTGTTAAAGGTTTACAACCATCTGTAGGTCCCAAACAAGCTAGTCTTAATGGAGTCACATTAATTGGAGAAGCCGTAAAAACACAGTTCAATTATTATAAGAAAGTTGATGGTGTTGTCCAACAATTACCTGAAACTTACTTTACTCAGAGTAGAAATTTACAAGAATTTAAACCCAGGAGTCAAATGGAAATTGATTTCTTAGAATTAGCTATGGATGAATTCATTGAACGGTATAAATTAGAAGGCTATGCCTTCGAACATATCGTTTATGGAGATTTTAGTCATAGTCAGTTAGGTGGTTTACATCTACTGATTGGACTAGCTAAACGTTTTAAGGAATCACCTTTTGAATTAGAAGATTTTATTCCTATGGACAGTACAGTTAAAAACTATTTCATAACAGATGCGCAAACAGGTTCATCTAAGTGTGTGTGTTCTGTTATTGATTTATTACTTGATGATTTTGTTGAAATAATAAAATCCCAAGATTTATCTGTAGTTTCTAAGGTTGTCAAAGTGACTATTGACTATACAGAAATTTCATTTATGCTTTGGTGTAAAGATGGCCATGTAGAAACATTTTACCCAAAATTACAATCTAGTCAAGCGTGGCAACCGGGTGTTGCTATGCCTAATCTTTACAAAATGCAAAGAATGCTATTAGAAAAGTGTGACCTTCAAAATTATGGTGATAGTGCAACATTACCTAAAGGCATAATGATGAATGTCGCAAAATATACTCAACTGTGTCAATATTTAAACACATTAACATTAGCTGTACCCTATAATATGAGAGTTATACATTTTGGTGCTGGTTCTGATAAAGGAGTTGCACCAGGTACAGCTGTTTTAAGACAGTGGTTGCCTACGGGTACGCTGCTTGTCGATTCAGATCTTAATGACTTTGTCTCTGATGCAGATTCAACTTTGATTGGTGATTGTGCAACTGTACATACAGCTAATAAATGGGATCTCATTATTAGTGATATGTACGACCCTAAGACTAAAAATGTTACAAAAGAAAATGACTCTAAAGAGGGTTTTTTCACTTACATTTGTGGGTTTATACAACAAAAGCTAGCTCTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGGGACACTTCGCATGGTGGACAGCCTTTGTTACTAATGTGAATGCGTCATCATCTGAAGCATTTTTAATTGGATGTAATTATCTTGGCAAACCACGCGAACAAATAGATGGTTATGTCATGCATGCAAATTACATATTTTGGAGGAATACAAATCCAATTCAGTTGTCTTCCTATTCTTTATTTGACATGAGTAAATTTCCCCTTAAATTAAGGGGTACTGCTGTTATGTCTTTAAAAGAAGGTCAAATCAATGATATGATTTTATCTCTTCTTAGTAAAGGTAGACTTATAATTAGAGAAAACAACAGAGTTGTTATTTCTAGTGATGTTCTTGTTAACAACTAAACGAACAATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTATACATGTCTCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAACGCTACTAATGTTGTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTGGGTGTTTATTACCACAAAAACAACAAAAGTTGGATGGAAAGTGAGTTCAGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAATTTGTGTTTAAGAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATTAATTTAGTGCGTGATCTCCCTCAGGGTTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTACTTGCTTTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGTGCTGCAGCTTATTATGTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATCTATCAAACTTCTAACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGACATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGATGTTAACTGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCCACAAATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGAATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCTTTAACTGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAATTTACAAAACACCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGATCCATCAAAACCAAGCAAGAGGTCATTTATTGAAGATCTACTTTTCAACAAAGTGACACTTGCAGATGCTGGCTTCATCAAACAATATGGTGATTGCCTTGGTGATATTGCTGCTAGAGACCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGAAATGATTGCTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAGGTGCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAGTTACACAGAATGTTCTCTATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTATTGGCAAAATTCAAGACTCACTTTCTTCCACAGCAAGTGCACTTGGAAAACTTCAAGATGTGGTCAACCAAAATGCACAAGCTTTAAACACGCTTGTTAAACAACTTAGCTCCAATTTTGGTGCAATTTCAAGTGTTTTAAATGATATCCTTTCACGTCTTGACAAAGTTGAGGCTGAAGTGCAAATTGATAGGTTGATCACAGGCAGACTTCAAAGTTTGCAGACATATGTGACTCAACAATTAATTAGAGCTGCAGAAATCAGAGCTTCTGCTAATCTTGCTGCTACTAAAATGTCAGAGTGTGTACTTGGACAATCAAAAAGAGTTGATTTTTGTGGAAAGGGCTATCATCTTATGTCCTTCCCTCAGTCAGCACCTCATGGTGTAGTCTTCTTGCATGTGACTTATGTCCCTGCACAAGAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACACTTTCCTCGTGAAGGTGTCTTTGTTTCAAATGGCACACACTGGTTTGTAACACAAAGGAATTTTTATGAACCACAAATCATTACTACAGACAACACATTTGTGTCTGGTAACTGTGATGTTGTAATAGGAATTGTCAACAACACAGTTTATGATCCTTTGCAACCTGAATTAGACTCATTCAAGGAGGAGTTAGATAAATATTTTAAGAATCATACATCACCAGATGTTGATTTAGGTGACATCTCTGGCATTAATGCTTCAGTTGTAAACATTCAAAAAGAAATTGACCGCCTCAATGAGGTTGCCAAGAATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTATGAGCAGTATATAAAATGGCCATGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGTAGTTGTCTCAAGGGCTGTTGTTCTTGTGGATCCTGCTGCAAATTTGATGAAGACGACTCTGAGCCAGTGCTCAAAGGAGTCAAATTACATTACACATAAACGAACTTATGGATTTGTTTATGAGAATCTTCACAATTGGAACTGTAACTTTGAAGCAAGGTGAAATCAAGGATGCTACTCCTTCAGATTTTGTTCGCGCTACTGCAACGATACCGATACAAGCCTCACTCCCTTTCGGATGGCTTATTGTTGGCGTTGCACTTCTTGCTGTTTTTCAGAGCGCTTCCAAAATCATAACCCTCAAAAAGAGATGGCAACTAGCACTCTCCAAGGGTGTTCACTTTGTTTGCAACTTGCTGTTGTTGTTTGTAACAGTTTACTCACACCTTTTGCTCGTTGCTGCTGGCCTTGAAGCCCCTTTTCTCTATCTTTATGCTTTAGTCTACTTCTTGCAGAGTATAAACTTTGTAAGAATAATAATGAGGCTTTGGCTTTGCTGGAAATGCCGTTCCAAAAACCCATTACTTTATGATGCCAACTATTTTCTTTGCTGGCATACTAATTGTTACGACTATTGTATACCTTACAATAGTGTAACTTCTTCAATTGTCATTACTTCAGGTGATGGCACAACAAGTCCTATTTCTGAACATGACTACCAGATTGGTGGTTATACTGAAAAATGGGAATCTGGAGTAAAAGACTGTGTTGTATTACACAGTTACTTCACTTCAGACTATTACCAGCTGTACTCAACTCAATTGAGTACAGACACTGGTGTTGAACATGTTACCTTCTTCATCTACAATAAAATTGTTGATGAGCCTGAAGAACATGTCCAAATTCACACAATCGACGGTTCATCCGGAGTTGTTAATCCAGTAATGGAACCAATTTATGATGAACCGACGACGACTACTAGCGTGCCTTTGTAAGCACAAGCTGATGAGTACGAACTTATGTACTCATTCGTTTCGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGTGGTATTCTTGCTAGTTACACTAGCCATCCTTACTGCGCTTCGATTGTGTGCGTACTGCTGCAATATTGTTAACGTGAGTCTTGTAAAACCTTCTTTTTACGTTTACTCTCGTGTTAAAAATCTGAATTCTTCTAGAGTTCCTGATCTTCTGGTCTAAACGAACTAAATATTATATTAGTTTTTCTGTTTGGAACTTTAATTTTAGCCATGGCAGATTCCAACGGTACTATTACCGTTGAAGAGCTTAAAAAGCTCCTTGAACAATGGAACCTAGTAATAGGTTTCCTATTCCTTACATGGATTTGTCTTCTACAATTTGCCTATGCCAACAGGAATAGGTTTTTGTATATAATTAAGTTAATTTTCCTCTGGCTGTTATGGCCAGTAACTTTAGCTTGTTTTGTGCTTGCTGCTGTTTACAGAATAAATTGGATCACCGGTGGAATTGCTATCGCAATGGCTTGTCTTGTAGGCTTGATGTGGCTCAGCTACTTCATTGCTTCTTTCAGACTGTTTGCGCGTACGCGTTCCATGTGGTCATTCAATCCAGAAACTAACATTCTTCTCAACGTGCCACTCCATGGCACTATTCTGACCAGACCGCTTCTAGAAAGTGAACTCGTAATCGGAGCTGTGATCCTTCGTGGACATCTTCGTATTGCTGGACACCATCTAGGACGCTGTGACATCAAGGACCTGCCTAAAGAAATCACTGTTGCTACATCACGAACGCTTTCTTATTACAAATTGGGAGCTTCGCAGCGTGTAGCAGGTGACTCAGGTTTTGCTGCATACAGTCGCTACAGGATTGGCAACTATAAATTAAACACAGACCATTCCAGTAGCAGTGACAATATTGCTTTGCTTGTACAGTAAGTGACAACAGATGTTTCATCTCGTTGACTTTCAGGTTACTATAGCAGAGATATTACTAATTATTATGAGGACTTTTAAAGTTTCCATTTGGAATCTTGATTACATCATAAACCTCATAATTAAAAATTTATCTAAGTCACTAACTGAGAATAAATATTCTCAATTAGATGAAGAGCAACCAATGGAGATTGATTAAACGAACATGAAAATTATTCTTTTCTTGGCACTGATAACACTCGCTACTTGTGAGCTTTATCACTACCAAGAGTGTGTTAGAGGTACAACAGTACTTTTAAAAGAACCTTGCTCTTCTGGAACATACGAGGGCAATTCACCATTTCATCCTCTAGCTGATAACAAATTTGCACTGACTTGCTTTAGCACTCAATTTGCTTTTGCTTGTCCTGACGGCGTAAAACACGTCTATCAGTTACGTGCCAGATCAGTTTCACCTAAACTGTTCATCAGACAAGAGGAAGTTCAAGAACTTTACTCTCCAATTTTTCTTATTGTTGCGGCAATAGTGTTTATAACACTTTGCTTCACACTCAAAAGAAAGACAGAATGATTGAACTTTCATTAATTGACTTCTATTTGTGCTTTTTAGCCTTTCTGCTATTCCTTGTTTTAATTATGCTTATTATCTTTTGGTTCTCACTTGAACTGCAAGATCATAATGAAACTTGTCACGCCTAAACGAACATGAAATTTCTTGTTTTCTTAGGAATCATCACAACTGTAGCTGCATTTCACCAAGAATGTAGTTTACAGTCATGTACTCAACATCAACCATATGTAGTTGATGACCCGTGTCCTATTCACTTCTATTCTAAATGGTATATTAGAGTAGGAGCTAGAAAATCAGCACCTTTAATTGAATTGTGCGTGGATGAGGCTGGTTCTAAATCACCCATTCAGTACATCGATATCGGTAATTATACAGTTTCCTGTTTACCTTTTACAATTAATTGCCAGGAACCTAAATTGGGTAGTCTTGTAGTGCGTTGTTCGTTCTATGAAGACTTTTTAGAGTATCATGACGTTCGTGTTGTTTTAGATTTCATCTAAACGAACAAACTAAAATGTCTGATAATGGACCCCAAAATCAGCGAAATGCACCCCGCATTACGTTTGGTGGACCCTCAGATTCAACTGGCAGTAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAACGTCGGCCCCAAGGTTTACCCAATAATACTGCGTCTTGGTTCACCGCTCTCACTCAACATGGCAAGGAAGACCTTAAATTCCCTCGAGGACAAGGCGTTCCAATTAACACCAATAGCAGTCCAGATGACCAAATTGGCTACTACCGAAGAGCTACCAGACGAATTCGTGGTGGTGACGGTAAAATGAAAGATCTCAGTCCAAGATGGTATTTCTACTACCTAGGAACTGGGCCAGAAGCTGGACTTCCCTATGGTGCTAACAAAGACGGCATCATATGGGTTGCAACTGAGGGAGCCTTGAATACACCAAAAGATCACATTGGCACCCGCAATCCTGCTAACAATGCTGCAATCGTGCTACAACTTCCTCAAGGAACAACATTGCCAAAAGGCTTCTACGCAGAAGGGAGCAGAGGCGGCAGTCAAGCCTCTTCTCGTTCCTCATCACGTAGTCGCAACAGTTCAAGAAATTCAACTCCAGGCAGCAGTAGGGGAACTTCTCCTGCTAGAATGGCTGGCAATGGCGGTGATGCTGCTCTTGCTTTGCTGCTGCTTGACAGATTGAACCAGCTTGAGAGCAAAATGTCTGGTAAAGGCCAACAACAACAAGGCCAAACTGTCACTAAGAAATCTGCTGCTGAGGCTTCTAAGAAGCCTCGGCAAAAACGTACTGCCACTAAAGCATACAATGTAACACAAGCTTTCGGCAGACGTGGTCCAGAACAAACCCAAGGAAATTTTGGGGACCAGGAACTAATCAGACAAGGAACTGATTACAAACATTGGCCGCAAATTGCACAATTTGCCCCCAGCGCTTCAGCGTTCTTCGGAATGTCGCGCATTGGCATGGAAGTCACACCTTCGGGAACGTGGTTGACCTACACAGGTGCCATCAAATTGGATGACAAAGATCCAAATTTCAAAGATCAAGTCATTTTGCTGAATAAGCATATTGACGCATACAAAACATTCCCACCAACAGAGCCTAAAAAGGACAAAAAGAAGAAGGCTGATGAAACTCAAGCCTTACCGCAGAGACAGAAGAAACAGCAAACTGTGACTCTTCTTCCTGCTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAGCAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGGCAGATGGGCTATATAAACGTTTTCGCTTTTCCGTTTACGATATATAGTCTACTCTTGTGCAGAATGAATTCTCGTAACTACATAGCACAAGTAGATGTAGTTAACTTTAATCTCACATAGCAATCTTTAATCAGTGTGTAACATTAGGGAGGACTTGAAAGAGCCACCACATTTTCACCGAGGCCACGCGGAGTACGATCGAGTGTACAGTGAACAATGCTAGGGAGAGCTGCCTATATGGAAGAGCCCTAATGTGTAAAATTAATTTTAGTAGTGCTATCCCCATGTGATTTTAATAGCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	var resultSeq = wuhan_hu_1;
	_.each(mods, function(mod) {
		resultSeq = applyMod(resultSeq, mod);
	});
	//glue.logInfo("testSeq", resultSeq);
	return resultSeq;
}



function applyMod(resultSeq, mod) {
	if(mod.type == "replaceSingle") {
		resultSeq = resultSeq.slice(0, mod.loc-1) + mod.replacement + resultSeq.slice(mod.loc);
	} else if(mod.type == "replaceRegion") {
		resultSeq = resultSeq.slice(0, mod.locStart-1) + mod.replacement + resultSeq.slice(mod.locEnd);
	} else if(mod.type == "insertRegion") {
		resultSeq = resultSeq.slice(0, mod.afterLoc) + mod.insertion + resultSeq.slice(mod.afterLoc);
	} 
	return resultSeq;
}

function reportMultiFasta(fastaFilePath) {
	var fastaDoc;
	glue.inMode("module/covFastaUtility", function() {
		fastaDoc = glue.command(["load-nucleotide-fasta", fastaFilePath]);
	});
	var issuesFileLines = ["multi_fasta_path\tfasta_id\tassay_name\tprimer_name\tprimer_sequence\tref_start\tref_end\tissue\tassay_publication\tassay_url\n"]
	var summaryFileLines = ["multi_fasta_path\tfasta_id\tnum_issues\n"]
	_.each(fastaDoc.nucleotideFasta.sequences, function(seq) {
		var queryName = seq.id;
		var queryNucleotides = seq.sequence;
		var resultDoc = singleSequenceReportAux(fastaFilePath, queryName, queryNucleotides);
		var numUnknown = 0;
		_.each(resultDoc.covPPReport.assays, function(assayObj) {
			_.each(assayObj.ppObjs, function(ppObj) {
				numUnknown += (ppObj.numIssues - ppObj.numKnownIssues);
				_.each(ppObj.unknownDisplayIssues, function(udi) {
					issuesFileLines.push(fastaFilePath+"\t"+queryName+"\t"+
							assayObj.display_name+"\t"+ppObj.display_name+"\t"+ppObj.sequence_to_scan+"\t"+
							ppObj.ref_start+"\t"+ppObj.ref_end+"\t"+
							udi+"\t"+assayObj.organisation+"\t"+assayObj.url+"\n");
				});
			});
		});
		summaryFileLines.push(fastaFilePath+"\t"+queryName+"\t"+numUnknown+"\n");
	});
	var issuesFilePath = fastaFilePath.substring(0, fastaFilePath.lastIndexOf("."))+"_ppReport_issues.txt";
	glue.command(["file-util", "save-string", issuesFileLines.join(""), issuesFilePath]);
	var summaryFilePath = fastaFilePath.substring(0, fastaFilePath.lastIndexOf("."))+"_ppReport_summary.txt";
	glue.command(["file-util", "save-string", summaryFileLines.join(""), summaryFilePath]);
}

