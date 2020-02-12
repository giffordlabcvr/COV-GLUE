var refSeqObjs = glue.tableToObjects(glue.command(["list", "reference", "name"]));

var codingFeaturesToCheck = glue.getTableColumn(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true"]), "name");

var problematicRefs = {};

var featureNameToMaxCodonNumber = {};

_.each(refSeqObjs, function(refSeqObj) {
	glue.inMode("reference/"+refSeqObj.name, function() {
		_.each(codingFeaturesToCheck, function(featureName) {
			glue.inMode("feature-location/"+featureName, function() {
				var aaRows = glue.tableToObjects(glue.command(["amino-acid"]));
				featureNameToMaxCodonNumber[featureName] = aaRows[aaRows.length - 1].codonLabel;
				
				glue.logInfo("Checking reference "+refSeqObj.name+", feature "+featureName+", "+aaRows.length+" amino acids.");
				for(var i = 0; i < aaRows.length; i++) {
					var aa = aaRows[i].aminoAcid;
					var nsp = featureName.indexOf("NSP") == 0;
					if(i == 0 && aa != "M" && !nsp) {
						glue.log("WARNING", "Residue "+aaRows[i].codonLabel+" of feature "+featureName+" on reference "+refSeqObj.name+" should be M");
						problematicRefs[refSeqObj.name] = "yes";
					}
					if(i < aaRows.length-1 && aa == "*") {
						glue.log("WARNING", "Residue "+aaRows[i].codonLabel+" of feature "+featureName+" on reference "+refSeqObj.name+" should not be *");
						problematicRefs[refSeqObj.name] = "yes";
					}
					if(i < aaRows.length-1 && aa == "X") {
						glue.log("WARNING", "Residue "+aaRows[i].codonLabel+" of feature "+featureName+" on reference "+refSeqObj.name+" should not be X");
						problematicRefs[refSeqObj.name] = "yes";
					}
					if(i == aaRows.length-1 && aa != "*" && !nsp) {
						glue.log("WARNING", "Residue "+aaRows[i].codonLabel+" of feature "+featureName+" on reference "+refSeqObj.name+" should be *");
						problematicRefs[refSeqObj.name] = "yes";
					}
				}
			});
			
		});
	});
});
_.each(_.pairs(featureNameToMaxCodonNumber), function(pair) {
	var featureName = pair[0];
	var maxCodonNumber = pair[1];
	glue.inMode("feature/"+featureName, function() {
		glue.command(["set", "field", "max_codon_number", maxCodonNumber]);
	})
});

var problematicReferenceSequences = _.keys(problematicRefs);
if(problematicReferenceSequences.length > 0) {
	glue.log("SEVERE", "problematicReferenceSequences", problematicReferenceSequences);
	throw new Error("Issues found in coding feature locations, see log for details");
}
