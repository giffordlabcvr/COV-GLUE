var refSeqObjs = glue.tableToObjects(glue.command(["list", "reference", "name"]));

var codingFeaturesToCheck = ["S", "ORF_8", "ORF_7a", "ORF_6", "ORF_3a", "ORF_1ab", "ORF_10", "N", "M", "E"];

var problematicRefs = {};

_.each(refSeqObjs, function(refSeqObj) {
	glue.inMode("reference/"+refSeqObj.name, function() {
		_.each(codingFeaturesToCheck, function(featureName) {
			glue.inMode("feature-location/"+featureName, function() {
				var aaRows = glue.tableToObjects(glue.command(["amino-acid"]));
				glue.logInfo("Checking reference "+refSeqObj.name+", feature "+featureName+", "+aaRows.length+" amino acids.");
				for(var i = 0; i < aaRows.length; i++) {
					var aa = aaRows[i].aminoAcid;
					if(i == 0 && aa != "M") {
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
					if(i == aaRows.length-1 && aa != "*") {
						glue.log("WARNING", "Residue "+aaRows[i].codonLabel+" of feature "+featureName+" on reference "+refSeqObj.name+" should be *");
						problematicRefs[refSeqObj.name] = "yes";
					}
				}
			});
			
		});
	});
});
var problematicReferenceSequences = _.keys(problematicRefs);
if(problematicReferenceSequences.length > 0) {
	glue.log("SEVERE", "problematicReferenceSequences", problematicReferenceSequences);
	throw new Error("Issues found in coding feature locations, see log for details");
}
