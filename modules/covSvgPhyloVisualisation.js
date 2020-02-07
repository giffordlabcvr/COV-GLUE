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

function visualisePhyloAsSvg(document) {
	var glueTree;
	
	var queryName = document.inputDocument.queryName;
	var placementIndex = document.inputDocument.placementIndex;
	var placerResult = document.inputDocument.placerResult;
	var queryNucleotides = document.inputDocument.queryNucleotides;
	var targetRefName = document.inputDocument.targetReferenceName;
	var queryToTargetRefSegs = document.inputDocument.queryToTargetRefSegments;

	var includeQuerySequence = false;
	if(queryName != null && placementIndex != null && placerResult != null && 
			queryNucleotides != null && targetRefName != null && queryToTargetRefSegs != null) {
		includeQuerySequence = true;
	}
	
	var aaVisFeatureName = document.inputDocument.aaVisFeature;
	var aaVisCodonLabel = document.inputDocument.aaVisCodonLabel;
	var aaVisDeletionStart = document.inputDocument.aaVisDeletionStart;
	var aaVisDeletionEnd = document.inputDocument.aaVisDeletionEnd;

	// map seqID to the AA residue at the requested loc
	var gisaidSeqIdToAA = {};
	// map seqID to string of multiple AA residue at the requested loc, used when the nucleotides are ambiguous
	var gisaidSeqIdToMultipleResidues = {};

	// map seqID to boolean, whether the sequence contains the requested deletion (or a deletion which contains the requested deletion)
	var gisaidSeqIdToDeletion = {};

	glue.logInfo("document.inputDocument", document.inputDocument);

	var almtMemberObjs;
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
		almtMemberObjs = glue.tableToObjects(glue.command(["list", "member"]));
		// generate a map of sequenceID to AA value for the GISAID sequences.
		_.each(almtMemberObjs, function(almtMemberObj) {
			var sequenceID = almtMemberObj["sequence.sequenceID"];
			glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+sequenceID, function() {
				
				if(aaVisCodonLabel != null) {
					var memberAa;
					var multipleResidues;
					var aaRows = glue.tableToObjects(glue.command(["amino-acid",
						"-r", "REF_MASTER_WUHAN_HU_1", "-f", aaVisFeatureName, 
						"-c", aaVisCodonLabel, aaVisCodonLabel]));
					if(aaRows.length == 0) {
						memberAa = "-";
						multipleResidues = "";
					} else {
						var aaObj = aaRows[0];
						memberAa = aaObj.aminoAcid;
						multipleResidues = "";

						if(memberAa == "X" && aaObj.definiteAas != null && aaObj.definiteAas != "" &&
								aaObj.definiteAas.length > 1 && aaObj.codonNts.indexOf('N') < 0) {
							memberAa = "?";
							multipleResidues = aaObj.definiteAas;
						}

					}
					gisaidSeqIdToAA[sequenceID] = memberAa;
					gisaidSeqIdToMultipleResidues[sequenceID] = multipleResidues;
				}
				if(aaVisDeletionStart && aaVisDeletionEnd != null) {
					var memberDelObjs = glue.tableToObjects(glue.command(["variation", "scan", 
						"-r", "REF_MASTER_WUHAN_HU_1", "-f", aaVisFeatureName,
						"--whereClause", "name = 'cov_aa_del_detect:"+aaVisFeatureName+"'", 
						"--excludeAbsent", "--showMatchesAsTable"]));
					gisaidSeqIdToDeletion[sequenceID] = false;
					_.each(memberDelObjs, function(memberDelObj) {
						if(parseInt(memberDelObj.refFirstCodonDeleted) <= parseInt(aaVisDeletionStart) &&
								parseInt(memberDelObj.refLastCodonDeleted) >= parseInt(aaVisDeletionEnd)) {
							gisaidSeqIdToDeletion[sequenceID] = true;
						}
					});
					
				}
			});
		});
	});
	
	
	if(includeQuerySequence) {
		// translate the aaVisFeature/CodonLabel for the submitted sequence
		var queryAa;
		var multipleResidues;
		glue.inMode("module/covFastaSequenceReporter", function() {
			var aaRows = glue.tableToObjects(glue.command({
				"string-plus-alignment": { 
					"amino-acid": {
						"fastaString": queryNucleotides,
						"queryToTargetSegs": {
							queryToTargetSegs: {
								alignedSegment: queryToTargetRefSegs
							}
						},
						"labelledCodon": true,
						"lcStart": aaVisCodonLabel,
						"lcEnd": aaVisCodonLabel,
						"targetRefName":targetRefName,
						"relRefName":"REF_MASTER_WUHAN_HU_1",
						"linkingAlmtName":"AL_GISAID_UNCONSTRAINED",
						"featureName":aaVisFeatureName
					}
				}
			}));
			if(aaRows.length == 0) {
				queryAa = "-";
				multipleResidues = "";
			} else {
				var aaObj = aaRows[0];
				queryAa = aaObj.aminoAcid;
				multipleResidues = "";

				if(queryAa == "X" && aaObj.definiteAas != null && aaObj.definiteAas != "" &&
						aaObj.definiteAas.length > 1 && aaObj.codonNts.indexOf('N') < 0) {
					queryAa = "?";
					multipleResidues = aaObj.definiteAas;
				}
				
			}
		});


		gisaidSeqIdToAA[queryName] = queryAa;
		gisaidSeqIdToMultipleResidues[queryName] = multipleResidues;
	}	

	if(includeQuerySequence) {
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
			// set query leaf node to highlighted
			glueTree = glue.command({
				"update-leaves" : {
					propertyName: "treevisualiser-highlighted",
					propertyValue: "true",
					leafNodeNames : [queryName], 
					inputPhyloTree: glueTree
				}
			});
			// set query leaf node to non-member
			glueTree = glue.command({
				"update-leaves" : {
					propertyName: "treevisualiser-nonmember",
					propertyValue: "true",
					leafNodeNames : [queryName], 
					inputPhyloTree: glueTree
				}
			});
			// set ancestor branches of query leaf node to highlighted
			glueTree = glue.command({
				"update-ancestor-branches" : {
					propertyName: "treevisualiser-highlighted",
					propertyValue: "true",
					leafNodeNames : [queryName], 
					inputPhyloTree: glueTree
				}
			});
		});
	} else {
		glue.inMode("module/covPhyloUtility", function() {
			glueTree = glue.command(["read-alignment-phylogeny", "AL_GISAID_UNCONSTRAINED", "phylogeny"]);
		});
	}

	// set treevisualiser-leafSourceName and treevisualiser-leafSequenceID throughout the tree.
	// this will cause TreeVisualiser to set leafSourceName and leafSequenceID properties on the leaf objects within visualiseTreeResult
	setLeafSeqIDs(glueTree.phyloTree.root);
	
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
					"leafTextAnnotationName": document.inputDocument.tipAnnotation
				}
			}
		});
	});

	// use the leafSequenceID to look up the AA value in the map, then set this as an additional property.
	_.each(visualiseTreeResult.visDocument.treeVisualisation.leafNodes, function(leafNode) {
		var seqID = leafNode.properties.leafSequenceID;
		if(aaVisCodonLabel != null) {
			leafNode.properties.aaValue = gisaidSeqIdToAA[seqID];
			var multipleResidues = gisaidSeqIdToMultipleResidues[seqID];
			if(multipleResidues.length > 0) {
				leafNode.properties.multipleResidues = multipleResidues;
			}
		}
		if(aaVisDeletionStart && aaVisDeletionEnd != null) {
			leafNode.properties.deletionPresent = gisaidSeqIdToDeletion[seqID];
		}
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