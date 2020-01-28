
function visualisePhyloAsSvg(document) {
	var glueTree;
	
	var queryName = document.inputDocument.queryName;
	var placementIndex = document.inputDocument.placementIndex;
	var placerResult = document.inputDocument.placerResult;
	
	// generate a tree for the placement, as a command document.
	glue.inMode("module/covMaxLikelihoodPlacer", function() {
		glueTree = glue.command({
				"export": {
					"placement-from-document": {
						"phylogeny": {
							"placerResultDocument": placerResult,
							"placementIndex": placementIndex,
							"queryName": queryName, 
							"leafName": queryName
						}
					}
				}
		});
	});

	glue.inMode("module/covPhyloUtility", function() {
		// set leaf node to highlighted
		glueTree = glue.command({
			"update-leaves" : {
				propertyName: "treevisualiser-highlighted",
				propertyValue: "true",
				leafNodeNames : [queryName], 
				inputPhyloTree: glueTree
			}
		});
		// set leaf node to non-member
		glueTree = glue.command({
			"update-leaves" : {
				propertyName: "treevisualiser-nonmember",
				propertyValue: "true",
				leafNodeNames : [queryName], 
				inputPhyloTree: glueTree
			}
		});
		// set ancestor branches of leaf node to highlighted
		glueTree = glue.command({
			"update-ancestor-branches" : {
				propertyName: "treevisualiser-highlighted",
				propertyValue: "true",
				leafNodeNames : [queryName], 
				inputPhyloTree: glueTree
			}
		});
	});

	
	// generate a visualisation document for the tree, 
	// with the visualisation maths etc. done
	var visualiseTreeResult;

	glue.inMode("module/covTreeVisualiser", function() {
		visualiseTreeResult = glue.command({
			"visualise" : {
				"tree-document": {
					"treeDocument" : glueTree, 
					"pxWidth" : document.inputDocument.pxWidth, 
					"pxHeight" : document.inputDocument.pxHeight,
					"legendPxWidth" : document.inputDocument.legendPxWidth, 
					"legendPxHeight" : document.inputDocument.legendPxHeight,
					"leafTextAnnotationName": "isolate"
				}
			}
		});
	});

	// from the visualisation documents, generate SVGs as GLUE web files.
	var treeTransformResult;
	glue.inMode("module/covTreeVisualisationTransformer", function() {
		treeTransformResult = glue.command({ "transform-to-web-file": 
			{
				"webFileType": "WEB_PAGE",
				"commandDocument":{
					transformerInput: {
						treeVisualisation: visualiseTreeResult.visDocument.treeVisualisation
					}
				},
				"outputFile": document.inputDocument.fileName
			}
		});
	});

	var legendTransformResult;
	glue.inMode("module/covTreeVisualisationLegendTransformer", function() {
		legendTransformResult = glue.command({ "transform-to-web-file": 
			{
				"webFileType": "WEB_PAGE",
				"commandDocument":{
					transformerInput: {
						treeVisualisationLegend: visualiseTreeResult.visDocument.treeVisualisationLegend
					}
				},
				"outputFile": document.inputDocument.legendFileName
			}
		});
	});

	glue.logInfo("legendTransformResult", legendTransformResult);
	
	return {
		visualisePhyloAsSvgResult: {
			treeTransformResult: treeTransformResult,
			legendTransformResult: legendTransformResult
		}
	}
}