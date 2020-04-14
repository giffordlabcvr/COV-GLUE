// recursive function to set the sequenceIDs throughout the tree.
function setLeafSeqIDs(subtree) {
	var userData;
	var subtreeType;
	if(subtree.internal != null) { // internal node
		userData = subtree.internal.userData;
		subtreeType = "internal";
	} else { // leaf node
		userData = subtree.leaf.userData;
		subtreeType = "leaf";
	}
	if(subtreeType == "leaf") {
		if(userData["name"].indexOf("cov-gisaid/") > 0) {
			userData["treevisualiser-leafSourceName"] = "cov-gisaid";
			userData["treevisualiser-leafSequenceID"] = userData["name"].substring(userData["name"].lastIndexOf("/")+1);
		} else {
			userData["treevisualiser-leafSourceName"] = "submitted";
			userData["treevisualiser-leafSequenceID"] = userData["name"];
		}
	} else { // internal node
		var branches = subtree.internal.branch;
		_.each(branches, function(branch) {
			setLeafSeqIDs(branch);
		});
	}
}

function setInternalLineage(lineageToParent, subtree) {
	var userData;
	var subtreeType;
	if(subtree.internal != null) { // internal node
		userData = subtree.internal.userData;
		if(userData == null) {
			userData = {};
			subtree.internal.userData = userData;
		}
		subtreeType = "internal";
		var branches = subtree.internal.branch;
		var childLineages = [];
		_.each(branches, function(branch) {
			childLineages.push(setInternalLineage(lineageToParent, branch));
		});
		var commonLineage = findCommonLineage(childLineages);
		_.each(childLineages, function(childLineage) {
			if(commonLineage != childLineage) {
				var oldParent = lineageToParent[childLineage];
				if(oldParent != null) {
					throw new Error("Lineage "+childLineage+" is parented in two locations "+oldParent+" and "+commonLineage);
				}
				lineageToParent[childLineage] = commonLineage;
			}
		});
		userData.lineage = commonLineage;
		return commonLineage;
	} else { // leaf node
		userData = subtree.leaf.userData;
		subtreeType = "leaf";
		var seqPath = userData["name"].replace("alignment/AL_GISAID_UNCONSTRAINED/member/", "");
		var leafLineage;
		glue.inMode("sequence/"+seqPath, function() {
			leafLineage = glue.command(["show", "property", "pang_lineage"]).propertyValueResult.value;
		});
		userData.lineage = leafLineage;
		return leafLineage;
	}
}

function checkLineages() {
	var glueTree;
	glue.inMode("module/covPhyloUtility", function() {
		glueTree = glue.command(["read-alignment-phylogeny", "AL_GISAID_UNCONSTRAINED", "phylogeny"]);
	});
	var lineageToParent = {};
	setInternalLineage(lineageToParent, glueTree.phyloTree.root);
	glue.logInfo("lineageToParent", lineageToParent);
	
}

function findCommonLineage(childLineages) {
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
	return commonLineageBits.join(".");
}

function assignLineages(document) {
	
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

	var queryNameToLineageResult = {};
	
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
			// somehow generate a lineage result for this placement.
		});
		
		// somehow combine the lineage results into a result for the query
	});
	
	return {
		covAssignLineagesResult: {
		}
	}
}