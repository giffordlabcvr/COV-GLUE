function setSubtreeLineage(subtree) {
	var userData;
	if(subtree.internal != null) { // internal node
		userData = subtree.internal.userData;
		if(userData == null) {
			userData = {};
			subtree.internal.userData = userData;
		}
		var branches = subtree.internal.branch;
		var childLineages = [];
		_.each(branches, function(branch) {
			childLineages.push(setSubtreeLineage(branch));
		});
		subTreeLineage = findCommonLineage(childLineages);
		userData.lineage = subTreeLineage;
		return subTreeLineage;
	} else { // leaf node
		userData = subtree.leaf.userData;
		var seqPath = userData["name"].replace("alignment/AL_GISAID_UNCONSTRAINED/member/", "");
		var leafLineage;
		glue.inMode("sequence/"+seqPath, function() {
			leafLineage = glue.command(["show", "property", "pang_lineage"]).propertyValueResult.value;
		});
		userData.lineage = leafLineage;
		return leafLineage;
	}
}

function findPlacementLineages(subtree, lineages, queryName) {
	if(subtree.internal != null) { // internal node
		userData = subtree.internal.userData;
		var lineage = null;
		if(userData != null) {
			lineage = userData.lineage;
			if(lineage == "") {lineage = "SARS-CoV-2";}
			if(lineage != null) {
				lineages.push(lineage);
			}
		}
		var branches = subtree.internal.branch;
		var lineagesResult = null;
		var resultBranch = null;
		var nonResultBranches = [];
		for(var i = 0; i < branches.length; i++) {
			var branch = branches[i];
			var branchResult = findPlacementLineages(branch, lineages, queryName);
			if(branchResult != null) {
				lineagesResult = branchResult;
				resultBranch = branch;
			} else {
				nonResultBranches.push(branch);
			}
		}
		// add the common lineage of any sister subtrees of the query taxon, that are on a near-zero length branch.
		if(resultBranch != null && resultBranch.leaf != null) {
			var additionalLineage = null;
			_.each(nonResultBranches, function(nonResultBranch) {
				var length = null;
				if(nonResultBranch.userData.length < 0.000006) { // this is approx 1/5 of a nucleotide in SARS-COV-2
					var nonResultLineage = null;
					if(nonResultBranch.internal != null) {
						nonResultLineage = nonResultBranch.internal.userData.lineage;
					} else {
						nonResultLineage = nonResultBranch.leaf.userData.lineage;
					}
					if(additionalLineage == null) {
						additionalLineage = nonResultLineage;
					} else {
						additionalLineage = findCommonLineage([additionalLineage, nonResultLineage]);
					}
				} 
			});
			if(additionalLineage != null) {
				if(additionalLineage == "") {additionalLineage = "SARS-CoV-2";}
				lineagesResult.push(additionalLineage);
			}
		}
		if(lineagesResult != null) {
			return lineagesResult;
		}
		if(lineage != null) {
			lineages.pop();
		}
		return null;
	} else { // leaf node
		userData = subtree.leaf.userData;
		if(userData != null && userData["name"] == queryName) {
			return lineages;
		}
		return null;
	}
}


function checkMonophyletic(current, completed, subtree) {
	var userData;
	var subtreeType;
	var newlyEntered;
	if(subtree.internal != null) { // internal node
		newlyEntered = checkLineageStatus(current, completed, subtree.internal.userData.lineage);
		_.each(subtree.internal.branch, function(branch) {
			checkMonophyletic(current, completed, branch);
		});
		//glue.log("FINEST", "Exiting '"+subtree.internal.userData.lineage+"'");
	} else {
		newlyEntered = checkLineageStatus(current, completed, subtree.leaf.userData.lineage);
	}
	_.each(newlyEntered, function(lineage) {
		delete current[lineage];
		glue.log("FINEST", "Exited lineage '"+lineage+"' ");
		completed[lineage] = "yes";
	});
}

function checkLineageStatus(current, completed, lineage) {
	var newlyEntered = [];
	var lineageBits = lineage.split(".");
	var ancestors = {};
	for(var i = 0; i <= lineageBits.length; i++) {
		var ancestorLineage = lineageBits.slice(0, i).join(".");
		ancestors[ancestorLineage] = "yes";
		if(completed[ancestorLineage]) {
			glue.log("INFO", "Warning: Lineage '"+ancestorLineage+"' is not monophyletic");
		}
	}
	_.each(_.keys(ancestors), function(lineage) {
		if(current[lineage] == null) {
			current[lineage] = "yes";
			glue.log("FINEST", "Entered lineage '"+lineage+"' ");
			newlyEntered.push(lineage);
		}
	});
	return newlyEntered;
}
	
function setLineagesInTree() {
	// set the lineage field on internal nodes and leaves in the tree.
	var glueTree;
	glue.inMode("module/covPhyloUtility", function() {
		glueTree = glue.command(["read-alignment-phylogeny", "AL_GISAID_UNCONSTRAINED", "phylogeny"]);
	});
	setSubtreeLineage(glueTree.phyloTree.root);
	// take all the internal nodes assigned to a given lineage X. 
	// there should be a single internal node such that all its internal node descendents
	// are assigned to lineage X or its sublineages
	checkMonophyletic({}, {}, glueTree.phyloTree.root);

	var glueTreeJson = JSON.stringify(glueTree);
	//glue.log("FINEST", "glueTreeJson", glueTreeJson);
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
		glue.command(["set", "field", "phylogeny", glueTreeJson]);
	});
	
}

function findCommonLineage(childLineages) {
	if(childLineages.length == 2) {
		if(isAncestorLineage(childLineages[0], childLineages[1])) {
			return childLineages[0];
		}
		if(isAncestorLineage(childLineages[1], childLineages[0])) {
			return childLineages[1];
		}
	}
	
	var commonLineageBits = null;
	
	_.each(childLineages, function(childLineage) {
		var childLineageBits = childLineage.split(".");
		if(commonLineageBits == null) {
			commonLineageBits = childLineageBits;
		} else {
			var i = 0;
			while(i < commonLineageBits.length && i < childLineageBits.length) {
				if(childLineageBits[i] != commonLineageBits[i]) {
					break;
				}
				i++;
			}
			commonLineageBits = commonLineageBits.slice(0, i);
		}
	});
	var commonLineage = commonLineageBits.join(".");
	return commonLineage;
}

function placeSequenceToFile(sequenceID, filePath) {
	var seqNts;
	glue.inMode("sequence/cov-gisaid/"+sequenceID, function() {
		seqNts = glue.command(["show", "nucleotides"]).nucleotidesResult.nucleotides;
	});
	var fastaDocument = { 
		nucleotideFasta: {
			sequences: [
				{
					id: sequenceID,
					sequence: seqNts
				}
			]
		}
	};
	var placerResult;
	glue.inMode("module/covMaxLikelihoodPlacer", function() {
		placerResult = glue.command({
			place: {
				"fasta-document": {
					fastaCommandDocument: fastaDocument
				}
			}
		});
	});
	var placerResultString = JSON.stringify(placerResult);
	glue.command(["file-util", "save-string", placerResultString, filePath]);

}

// note this is for testing, the input is a json, not XML.
function assignLineagesFromPlacerFile(filePath) {
	var fileString = glue.command(["file-util", "load-string", filePath]).fileUtilLoadStringResult.loadedString;
	var placerResult = JSON.parse(fileString);
	return assignLineagesFromPlacerDocument(placerResult);
}

// return true if lineage1 is an ancestor of lineage2 (lineages are their own ancestors)
function isAncestorLineage(lineage1, lineage2) {
	if(lineage1 == "SARS-CoV-2") {
		return true;
	}
	var lineage1Bits = lineage1.split(".")
	var lineage2Bits = lineage2.split(".")
	if(lineage1Bits.length <= lineage2Bits.length) {
		for(var i = 0; i < lineage1Bits.length; i++) {
			if(lineage1Bits[i] != lineage2Bits[i]) {
				return false;
			}
		}
		return true;
	}
	return false;
}

function assignLineagesDocumentForSequenceBatch(whereClause, offset, batchSize) {
	var fastaDocument;
	glue.inMode("module/covAssignLineagesFastaExporter", function() {
		fastaDocument = glue.command(["export", "-w", whereClause, "--offset", 
			offset, "--batchSize", batchSize, "--preview"]);
	});
	return assignLineagesDocumentForFastaDocument(fastaDocument);
}

function assignLineagesForSequenceBatch(whereClause, offset, batchSize) {
	return assignLineagesDocumentToObjectList(assignLineagesDocumentForSequenceBatch(whereClause, offset, batchSize));
}

function assignLineagesDocumentForSequences(whereClause) {
	var fastaDocument;
	glue.inMode("module/covAssignLineagesFastaExporter", function() {
		fastaDocument = glue.command(["export", "-w", whereClause, "--preview"]);
	});
	return assignLineagesDocumentForFastaDocument(fastaDocument);
}

function assignLineagesForSequences(whereClause) {
	return assignLineagesDocumentToObjectList(assignLineagesDocumentForSequences(whereClause));
}

function assignLineagesDocumentForFile(fastaFilePath) {
	var fastaDocument;
	glue.inMode("module/covFastaUtility", function() {
		fastaDocument = glue.command(["load-nucleotide-fasta", fastaFilePath]);
	});
	return assignLineagesDocumentForFastaDocument(fastaDocument);
}

function assignLineagesForFile(fastaFilePath) {
	return assignLineagesDocumentToObjectList(assignLineagesDocumentForFile(fastaFilePath));
}


function assignLineagesDocumentForFastaDocument(fastaDocument) {
	var placerDocument;
	glue.inMode("module/covMaxLikelihoodPlacer", function() {
		placerDocument = glue.command({
			"place": {
				"fasta-document": {
					"fastaCommandDocument": fastaDocument
				}
			}
		});
	});
	return assignLineagesFromPlacerDocument(placerDocument);
}

function assignLineagesForFastaDocument(fastaDocument) {
	return assignLineagesDocumentToObjectList(assignLineagesDocumentForFastaDocument(fastaDocument));
}


function assignLineagesDocumentToObjectList(documentResult) {
	glue.log("FINEST", "documentResult", documentResult);
	return _.map(documentResult.covAssignLineagesResult.queryLineageResults, 
			function(qlr) { return { 
				queryName: qlr.queryName, 
				lineage: qlr.bestLineage, 
				likelihoodWeightRatio: qlr.bestLikelihoodWeightRatio
			} });
}

function assignLineagesFromPlacerDocument(placerResult) {
	
	var queryObjs;

	glue.inMode("module/covMaxLikelihoodPlacer", function() {
		queryObjs = glue.tableToObjects(glue.command({
			"list": {
				"query-from-document": {
					"placerResultDocument": placerResult
				}
			}
		}));
	});

	var queryLineageResults = [];
	
	_.each(queryObjs, function(queryObj) {
		var queryName = queryObj.queryName;
		
		var placementObjs;
		glue.inMode("module/covMaxLikelihoodPlacer", function() {
			placementObjs = glue.tableToObjects(glue.command({
				"list": {
					"placement-from-document": {
						"placerResultDocument": placerResult,
						"queryName": queryName
					}
				}
			}));
		});
		
		var placementLineageResults = [];
		_.each(placementObjs, function(placementObj) {
			// generate a tree for the placement, as a command document.
			var glueTree;
			glue.inMode("module/covMaxLikelihoodPlacer", function() {
				glueTree = glue.command({
						"export": {
							"placement-from-document": {
								"phylogeny": {
									"placerResultDocument": placerResult,
									"placementIndex": placementObj.placementIndex,
									"queryName": queryName, 
									"leafName": queryName
								}
							}
						}
				});
			});
			// write out the tree for a specific placement, useful for debugging.
			// good online JSON browser here: https://codebeautify.org/jsonviewer
			//if(placementObj.placementIndex == 1) {
			//	glue.command(["file-util", "save-string", JSON.stringify(glueTree, null, 2), "tree.json"]);
			//}
			// generate a lineage result for this placement.
			var placementLineages = findPlacementLineages(glueTree.phyloTree.root, [], queryName);
			var placementLineageResult = {
					placementIndex: placementObj.placementIndex,
					likeWeightRatio: placementObj.likeWeightRatio,
					lineages: _.uniq(placementLineages)
			};
			placementLineageResults.push(placementLineageResult);
		});
		// somehow combine the lineage results into a result for the query
		var lineageFractionsMap = {};
		_.each(placementLineageResults, function(placementLineageResult) {
			_.each(placementLineageResult.lineages, function(lineage) {
				var oldFraction = lineageFractionsMap[lineage];
				if(oldFraction == null) {
					oldFraction = 0.0;
				}
				var newFraction = oldFraction + placementLineageResult.likeWeightRatio;
				lineageFractionsMap[lineage] = newFraction;
			});
		});
		var lineageFractions = [];
		_.each(_.pairs(lineageFractionsMap), function(pair) {
			var lineage = pair[0];
			if(lineage == "") {
				lineage = "SARS-CoV-2";
			}
			lineageFractions.push({
				lineage: lineage,
				totalLikelihoodWeightRatio: pair[1]
			});
		});
		var bestLineage = null;
		var bestLikelihoodWeightRatio;
		// best lineage is defined as follows:
		// take the set of lineages with totalLikelihoodWeightRatio > 50%
		// which of these is not an ancestor of any other in the set
		var candidateLineageFractions = _.filter(lineageFractions, function(lf) {
			return lf.totalLikelihoodWeightRatio > 0.5;
		});
		_.each(candidateLineageFractions, function(clf) {
			if(bestLineage == null || isAncestorLineage(bestLineage, clf.lineage)) {
				bestLineage = clf.lineage;
				bestLikelihoodWeightRatio = clf.totalLikelihoodWeightRatio;
			}
		});
		
		queryLineageResults.push({
				queryName: queryName,
				placementLineageResults: placementLineageResults,
				lineageFractions: lineageFractions,
				bestLineage: bestLineage,
				bestLikelihoodWeightRatio: bestLikelihoodWeightRatio
		});
		
	});
	
	return {
		covAssignLineagesResult: {
			queryLineageResults: queryLineageResults
		}
	}
}

