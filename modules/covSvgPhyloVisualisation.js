// recursive function to set the treevisualiser hints throughout the tree.
function setTreeVisualiserHints(subtree, gisaidSeqIdToPopupLines) {
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
			var sequenceID = userData["name"].substring(userData["name"].lastIndexOf("/")+1);
			userData["treevisualiser-leafSequenceID"] = sequenceID;
			var popupLines = gisaidSeqIdToPopupLines[sequenceID];
			if(popupLines != null && popupLines.length > 0) {
				userData["treevisualiser-leafPopupLines"] = popupLines;
			}
		} else if(userData["name"].indexOf("cov-coguk-refs/") > 0) {
			userData["treevisualiser-leafSourceName"] = "cov-coguk-refs";
			var sequenceID = userData["name"].substring(userData["name"].lastIndexOf("/")+1);
			userData["treevisualiser-leafSequenceID"] = sequenceID;
			var popupLines = gisaidSeqIdToPopupLines[sequenceID];
			if(popupLines != null && popupLines.length > 0) {
				userData["treevisualiser-leafPopupLines"] = popupLines;
			}
		} else {
			userData["treevisualiser-leafSourceName"] = "submitted";
			userData["treevisualiser-leafSequenceID"] = userData["name"];
		}
	} else { // internal node
		var branches = subtree.internal.branch;
		_.each(branches, function(branch) {
			setTreeVisualiserHints(branch, gisaidSeqIdToPopupLines);
		});
	}
}

function wordWrap(longString, maxChars) {
	var lines = [];
	var remainingString = longString.trim();
	while(remainingString.length > 0) {
		if(remainingString.length <= maxChars) {
			lines.push(remainingString);
			break;
		} else {
			var lastSpaceIndex = null;
			for(var i = maxChars-1; i >= 0; i--) {
				if(remainingString[i] == ' ') {
					lastSpaceIndex = i;
					break;
				}
			}
			if(lastSpaceIndex != null) {
				lines.push(remainingString.substring(0, lastSpaceIndex));
				remainingString = remainingString.substring(lastSpaceIndex).trim();
			} else {
				lines.push(remainingString.substring(0, maxChars-2)+"-");
				remainingString = remainingString.substring(maxChars-2).trim();
			}
		}
	}
	return lines;
}

function visualisePhyloAsSvg(document) {
	var glueTree;
	
	var queryName = document.inputDocument.queryName;
	var placementIndex = document.inputDocument.placementIndex;
	var placerResult = document.inputDocument.placerResult;
	
	var webFileType = document.inputDocument.webFileType;
	
	var queryNucleotides = document.inputDocument.queryNucleotides;
	var targetRefName = document.inputDocument.targetReferenceName;
	var queryToTargetRefSegs = document.inputDocument.queryToTargetRefSegments;

	var includeQuerySequence = false;
	if(queryName != null && placementIndex != null && placerResult != null && 
			queryNucleotides != null && targetRefName != null && queryToTargetRefSegs != null) {
		includeQuerySequence = true;
	}
	
	var visFeatureName = document.inputDocument.visFeature;
	var aaVisCodonLabel = document.inputDocument.aaVisCodonLabel;
	var visDeletionStart = document.inputDocument.visDeletionStart;
	var visDeletionEnd = document.inputDocument.visDeletionEnd;
	var visInsertionLastNtBeforeStart = document.inputDocument.visInsertionLastNtBeforeStart;
	var visInsertionFirstNtAfterEnd = document.inputDocument.visInsertionFirstNtAfterEnd;
	var visInsertedNts = document.inputDocument.visInsertedNts;

	// map seqID to the AA residue at the requested loc
	var gisaidSeqIdToAA = {};
	// map seqID to string of multiple AA residue at the requested loc, used when the nucleotides are ambiguous
	var gisaidSeqIdToMultipleResidues = {};

	// map seqID to string code
	// "deletion" the sequence contains the requested deletion (or a deletion which contains the requested deletion)
	// "no_deletion" the sequence does not contain the requested deletion
	// "insufficient_coverage" insufficient coverage to detect whether or not.
	var gisaidSeqIdToDeletion = {};

	// map seqID to string code
	// "insertion" the sequence contains the requested insertion
	// "insertion_different_nts" the sequence contains the requested insertion, but different NTs
	// "no_insertion" the sequence does not contain the requested insertion
	// "insufficient_coverage" insufficient coverage to detect whether or not.
	var gisaidSeqIdToInsertion = {};
	// map seqID to inserted AAs, populated when code is "insertion_different_nts"
	var gisaidSeqIdToInsertedNTs = {};

	var gisaidSeqIdToPopupLines = {};
	
	glue.logInfo("document.inputDocument", document.inputDocument);

	var almtMemberObjs;
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
		almtMemberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", "sequence.include_in_ref_tree = true"]));
		// generate a map of sequenceID to AA value for the sequences.
		_.each(almtMemberObjs, function(almtMemberObj) {
			var sequenceID = almtMemberObj["sequence.sequenceID"];
			glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+sequenceID, function() {
				
				if(aaVisCodonLabel != null) {
					var memberAa;
					var multipleResidues;
					var aaRows = glue.tableToObjects(glue.command(["amino-acid",
						"-r", "REF_MASTER_WUHAN_HU_1", "-f", visFeatureName, 
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
				// "insufficient_coverage" not generated yet
				if(visDeletionStart != null && visDeletionEnd != null) {
					var memberDelObjs = glue.tableToObjects(glue.command(["variation", "scan", 
						"-r", "REF_MASTER_WUHAN_HU_1", "-f", visFeatureName,
						"--whereClause", "name = 'cov_aa_del_detect:"+visFeatureName+"'", 
						"--excludeAbsent", "--showMatchesAsTable"]));
					gisaidSeqIdToDeletion[sequenceID] = "no_deletion";
					for(var i = 0; i < memberDelObjs.length; i++) {
						var memberDelObj = memberDelObjs[i];
						if(parseInt(memberDelObj.refFirstNtDeleted) <= parseInt(visDeletionStart) &&
								parseInt(memberDelObj.refLastNtDeleted) >= parseInt(visDeletionStart)) {
							gisaidSeqIdToDeletion[sequenceID] = "deletion";
							break;
						}
					}
					
				}
				// "insufficient_coverage" not generated yet
				if(visInsertionLastNtBeforeStart != null && visInsertionFirstNtAfterEnd != null && visInsertedNts != null) {
					var memberInsObjs = glue.tableToObjects(glue.command(["variation", "scan", 
						"-r", "REF_MASTER_WUHAN_HU_1", "-f", visFeatureName,
						"--whereClause", "name = 'cov_aa_ins_detect:"+visFeatureName+"'", 
						"--excludeAbsent", "--showMatchesAsTable"]));
					gisaidSeqIdToInsertion[sequenceID] = "no_insertion";
					for(var i = 0; i < memberInsObjs.length; i++) {
						var memberInsObj = memberInsObjs[i];
						if(parseInt(memberInsObj.refLastNtBeforeIns) == parseInt(visInsertionLastNtBeforeStart) &&
								parseInt(memberInsObj.refFirstNtAfterIns) == parseInt(visInsertionFirstNtAfterEnd)) {
							if(memberInsObj.insertedQryNts == visInsertedNts) {
								gisaidSeqIdToInsertion[sequenceID] = "insertion";
							} else {
								gisaidSeqIdToInsertion[sequenceID] = "insertion_different_nts";
								gisaidSeqIdToInsertedNTs[sequenceID] = memberInsObj.insertedQryNts;
							}
							break;
						}
					}
					
				}
			});
		});
	});

	
	var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", 
		"include_in_ref_tree = true", "source.name", "sequenceID", "collection_date", "m49_country.display_name", 
		"gisaid_authors_short", "gisaid_originating_lab", "gisaid_submitting_lab"]));

	var maxChars = 50; // word wrapping param.
	
	_.each(seqObjs, function(seqObj) {
		var sequenceID = seqObj["sequenceID"];
		var sourceName = seqObj["source.name"];
		var date = seqObj["collection_date"];
		var country = seqObj["m49_country.display_name"];
		var authors = seqObj["gisaid_authors_short"];
		var originatingLab = seqObj["gisaid_originating_lab"];
		var submittingLab = seqObj["gisaid_submitting_lab"];
		var popupLines = [];
		if(sourceName == "cov-gisaid") {
			popupLines.push("GISAID: "+sequenceID);
		} else if(sourceName == "cov-coguk-refs") {
			popupLines.push("COG-UK: "+sequenceID);
		}
		if(date != null) {
			popupLines.push("Collected: "+date);
		}
		if(country != null) {
			popupLines.push("Country: "+country);
		}
		if(originatingLab != null) {
			popupLines = popupLines.concat(wordWrap("Originating lab: "+originatingLab, maxChars))
		}
		if(submittingLab != null) {
			popupLines = popupLines.concat(wordWrap("Submitting lab: "+submittingLab, maxChars))
		}
		if(authors != null) {
			popupLines.push("Authors: "+authors);
		}
		gisaidSeqIdToPopupLines[sequenceID] = popupLines;
	})
	
	
	if(includeQuerySequence) {
		// translate the visFeature/CodonLabel for the submitted sequence
		var queryAa;
		var multipleResidues;
		glue.inMode("module/covFastaSequenceReporter", function() {
			if(aaVisCodonLabel != null) {
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
							"linkingAlmtName":"AL_UNCONSTRAINED_DUMMY",
							"featureName":visFeatureName
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
			}
			if(visDeletionStart != null && visDeletionEnd != null) {
				var queryDelObjs = glue.tableToObjects(glue.command({
					"string-plus-alignment": { 
						"variation": {
							"scan" :{
								"fastaString": queryNucleotides,
								"queryToTargetSegs": {
									queryToTargetSegs: {
										alignedSegment: queryToTargetRefSegs
									}
								},
								"whereClause": "name = 'cov_aa_del_detect:"+visFeatureName+"'",
								"excludeAbsent": true,
								"showMatchesAsTable": true,
								"showMatchesAsDocument": false,
								"targetRefName":targetRefName,
								"relRefName":"REF_MASTER_WUHAN_HU_1",
								"linkingAlmtName":"AL_UNCONSTRAINED_DUMMY",
								"featureName":visFeatureName
							}
						}
					}
				}));
				gisaidSeqIdToDeletion[queryName] = "no_deletion";
				_.each(queryDelObjs, function(queryDelObj) {
					if(parseInt(queryDelObj.refFirstNtDeleted) <= parseInt(visDeletionStart) &&
							parseInt(queryDelObj.refLastNtDeleted) >= parseInt(visDeletionEnd)) {
						gisaidSeqIdToDeletion[queryName] = "deletion";
					}
				});
			}
			if(visInsertionLastNtBeforeStart != null && visInsertionFirstNtAfterEnd != null && visInsertedNts != null) {
				var queryInsObjs = glue.tableToObjects(glue.command({
					"string-plus-alignment": { 
						"variation": {
							"scan" :{
								"fastaString": queryNucleotides,
								"queryToTargetSegs": {
									queryToTargetSegs: {
										alignedSegment: queryToTargetRefSegs
									}
								},
								"whereClause": "name = 'cov_aa_ins_detect:"+visFeatureName+"'",
								"excludeAbsent": true,
								"showMatchesAsTable": true,
								"showMatchesAsDocument": false,
								"targetRefName":targetRefName,
								"relRefName":"REF_MASTER_WUHAN_HU_1",
								"linkingAlmtName":"AL_UNCONSTRAINED_DUMMY",
								"featureName":visFeatureName
							}
						}
					}
				}));
				gisaidSeqIdToInsertion[queryName] = "no_insertion";
				for(var i = 0; i < queryInsObjs.length; i++) {
					var queryInsObj = queryInsObjs[i];
					if(parseInt(queryInsObj.refLastNtBeforeIns) == parseInt(visInsertionLastNtBeforeStart) &&
							parseInt(queryInsObj.refFirstNtAfterIns) == parseInt(visInsertionFirstNtAfterEnd)) {
						if(queryInsObj.insertedQryNts == visInsertedNts) {
							gisaidSeqIdToInsertion[queryName] = "insertion";
						} else {
							gisaidSeqIdToInsertion[queryName] = "insertion_different_nts";
							gisaidSeqIdToInsertedNTs[queryName] = queryInsObj.insertedQryNts;
						}
						break;
					}
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

	// set treevisualiser-leafSourceName, treevisualiser-leafSequenceID and treevisualiser-leafPopupLines
	// throughout the tree.
	// this will cause TreeVisualiser to set various properties on the leaf objects within visualiseTreeResult
	setTreeVisualiserHints(glueTree.phyloTree.root, gisaidSeqIdToPopupLines);
	
	// generate a visualisation document for the tree, 
	// with the visualisation maths etc. done
	var visualiseTreeResult;

	/* hacky way to dynamically change the pxHeight based on the number of taxa, so the text does not go too small. */
	var pxHeight = document.inputDocument.pxHeight;
	if(pxHeight == "auto") {
		var numTaxa = glue.command(["count", "sequence", "-w", "include_in_ref_tree = true"]).countResult.count;
		var leafHeightPx = 15;
		pxHeight = Math.ceil((numTaxa * leafHeightPx) / 0.96); // the division here is to account for the 2% top/bottom margin which the tree visualiser will apply.
	}
	
	glue.inMode("module/covTreeVisualiser", function() {
		visualiseTreeResult = glue.command({
			"visualise" : {
				"tree-document": {
					"treeDocument" : glueTree, 
					"pxWidth" : document.inputDocument.pxWidth, 
					"pxHeight" : pxHeight,
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
		if(visDeletionStart != null && visDeletionEnd != null) {
			leafNode.properties.deletionCode = gisaidSeqIdToDeletion[seqID];
		}
		if(visInsertionLastNtBeforeStart != null && visInsertionFirstNtAfterEnd != null && visInsertedNts != null) {
			leafNode.properties.insertionCode = gisaidSeqIdToInsertion[seqID];
			var insertedNTs = gisaidSeqIdToInsertedNTs[seqID];
			if(insertedNTs != null) {
				leafNode.properties.insertedNTs = insertedNTs;
			}
		}
		var popupLines = gisaidSeqIdToPopupLines[seqID];
		if(popupLines != null) {
			leafNode.properties.popupLines = popupLines;
		}
	});
	
	var includeLegendInTree = false;
	if(document.inputDocument.legendFileName == null) {
		includeLegendInTree = true;
	}
	
	// from the visualisation documents, generate SVGs as GLUE web files.
	var treeTransformResult;
	glue.inMode("module/covTreeVisualisationTransformer", function() {
		var transformerInput = {
			treeVisualisation: visualiseTreeResult.visDocument.treeVisualisation
		};
		if(includeLegendInTree) {
			transformerInput.treeVisualisationLegend = visualiseTreeResult.visDocument.treeVisualisationLegend
		}
		treeTransformResult = glue.command({ "transform-to-web-file": 
			{
				"webFileType": webFileType,
				"commandDocument":{
					transformerInput: transformerInput
				},
				"outputFile": document.inputDocument.fileName
			}
		});
	});

	var visualisePhyloAsSvgResult = {
		treeTransformResult: treeTransformResult
	};
	
	if(!includeLegendInTree) {
		var legendTransformResult;
		glue.inMode("module/covTreeVisualisationLegendTransformer", function() {
			legendTransformResult = glue.command({ "transform-to-web-file": 
				{
					"webFileType": webFileType,
					"commandDocument":{
						transformerInput: {
							treeVisualisationLegend: visualiseTreeResult.visDocument.treeVisualisationLegend
						}
					},
					"outputFile": document.inputDocument.legendFileName
				}
			});
		});
	}
	visualisePhyloAsSvgResult.legendTransformResult = legendTransformResult;
	
	return {
		visualisePhyloAsSvgResult: visualisePhyloAsSvgResult
	}
}