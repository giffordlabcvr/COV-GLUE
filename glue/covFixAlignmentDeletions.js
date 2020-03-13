/*
if non-codon-aligned deletions are detected at this step, we try to fix them if the deletion size is a multiple of 3.
This is because they are probably due to deletions which MAFFT did not properly codon-align.
Some examples:
nsp15:227:del in EPI_ISL_408670. 
nsp1:82-86:del in EPI_ISL_410044, EPI_ISL_413623
nsp2:268del in EPI_ISL_413564, EPI_ISL_413568, EPI_ISL_413573, EPI_ISL_413577, 
               EPI_ISL_413580, EPI_ISL_413581, EPI_ISL_413582, EPI_ISL_413583
S:145del in EPI_ISL_413522
the fix involves 
(a) removing the segments of the sequence from AL_GISAID_UNCONSTRAINED 
(b) adding the sequence to AL_GISAID_CONSTRAINED 
(c) recomputing this row of AL_GISAID_CONSTRAINED using the codon-aware compound aligner.
(d) deriving the segments back to AL_GISAID_UNCONSTRAINED.
*/


glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del_detect%'"]);

var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));

// create "detection" variation objects to detect deletions in each of the features respectively.
var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
_.each(featuresList, function(featureObj) {
	glue.inMode("reference/"+comparisonRefName+"/feature-location/"+featureObj.name, function() {
		var codonLabels= glue.getTableColumn(glue.command(["list", "labeled-codon"]), "codonLabel");
		var codonLabel2 = codonLabels[1];
		var codonLabelPenultimate = codonLabels[codonLabels.length - 2];
		glue.command(["create", "variation", "cov_aa_del_detect:"+featureObj.name, 
			"-t", "aminoAcidDeletion", "-c", codonLabel2, codonLabelPenultimate]);
	});
});

var realignedMemberObjs = [];

glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
	var almtMemberObjs = glue.tableToObjects(glue.command(["list", "member"]));
	_.each(almtMemberObjs, function(almtMemberObj) {
		glue.inMode("member/"+almtMemberObj["sequence.source.name"]+"/"+almtMemberObj["sequence.sequenceID"], function() {
			
			var allFeatureDelObjs = [];
			
			_.each(featuresList, function(featureObj) {
				var memberDelObjs = glue.tableToObjects(glue.command(["variation", "scan", 
					"-r", comparisonRefName, "-f", featureObj.name, 
					"--whereClause", "name = 'cov_aa_del_detect:"+featureObj.name+"'", 
					"--excludeAbsent", "--showMatchesAsTable"]));
				
				allFeatureDelObjs = allFeatureDelObjs.concat(memberDelObjs);
			});
			var realign = false;
			_.each(allFeatureDelObjs, function(memberDelObj) {
				var codonAligned = memberDelObj.deletionIsCodonAligned;
				if((!codonAligned) && memberDelObj.deletedRefNts.length % 3 == 0) { // deletion is multiple of 3, so attempt realignment.
					realign = true;
				}
			});
			if(realign) {
				realignedMemberObjs.push(almtMemberObj);
			}
		});
	}); 
});

_.each(realignedMemberObjs, function(realignedMemberObj) {
	var sourceName = realignedMemberObj["sequence.source.name"];
	var sequenceID = realignedMemberObj["sequence.sequenceID"];
	
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED/member/"+sourceName+"/"+sequenceID, function() {
		glue.command(["remove", "segment", "-a"]);
	});

	glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
		glue.command(["add", "member", sourceName, sequenceID]);
	});

	glue.command(["compute", "alignment", "AL_GISAID_CONSTRAINED", "covCompoundAligner", 
		"-w", "sequence.source.name = '"+sourceName+"' and sequence.sequenceID = '"+sequenceID+"'"]);

	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED", function() {
		glue.command(["derive", "segments", "AL_GISAID_CONSTRAINED", 
			"-w", "sequence.source.name = '"+sourceName+"' and sequence.sequenceID = '"+sequenceID+"'"]);
	});
});

