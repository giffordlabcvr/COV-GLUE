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
			if(lineage != null) {
				lineages.push(lineage);
			}
		}
		var branches = subtree.internal.branch;
		for(var i = 0; i < branches.length; i++) {
			var branch = branches[i];
			var branchResult = findPlacementLineages(branch, lineages, queryName);
			if(branchResult != null) {
				return branchResult;
			}
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
		//glue.logInfo("Entering '"+subtree.internal.userData.lineage+"'");
		_.each(subtree.internal.branch, function(branch) {
			checkMonophyletic(current, completed, branch);
		});
		//glue.logInfo("Exiting '"+subtree.internal.userData.lineage+"'");
	} else {
		newlyEntered = checkLineageStatus(current, completed, subtree.leaf.userData.lineage);
	}
	_.each(newlyEntered, function(lineage) {
		delete current[lineage];
		glue.logInfo("Exited lineage '"+lineage+"' ");
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
			glue.logInfo("Entered lineage '"+lineage+"' ");
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
	//glue.logInfo("glueTreeJson", glueTreeJson);
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

function assignLineagesFromPlacerFile(filePath) {
	var fileString = glue.command(["file-util", "load-string", filePath]).fileUtilLoadStringResult.loadedString;
	var placerResult = JSON.parse(fileString);
	var document = {
			inputDocument: {
				placerResult: placerResult
			}
	};
	return assignLineagesFromPlacerDocument(document);
}

// return true if lineage1 is an ancestor of lineage2 (lineages are their own ancestors)
function isAncestorLineage(lineage1, lineage2) {
	if(lineage2.indexOf(lineage1) == 0) {
		return true;
	}
	if(lineage1 == "B.7" && lineage2 == "B.10") {
		return true;
	}
	return false;
}

function assignLineagesFromPlacerDocument(document) {
	
	var placerResult = document.inputDocument.placerResult;
	
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
			// generate a lineage result for this placement.
			var placementLineages = findPlacementLineages(glueTree.phyloTree.root, [], queryName);
			var placementLineageResult = {
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
			lineageFractions.push({
				lineage: pair[0],
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

