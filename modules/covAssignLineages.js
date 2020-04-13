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