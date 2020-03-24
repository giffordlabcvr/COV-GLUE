glue.command(["multi-delete", "variation", "-w", "name like 'cov:clvg%'"]);

glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/ORF_1ab", function() {
	var variationName = "cov:clvg_papain";
	glue.command(["create", "variation", variationName, 
		"-t", "aminoAcidRegexPolymorphism", 
		"--labeledCodon", 1, 7097]);
	glue.inMode("variation/"+variationName, function() {
		glue.command(["set", "metatag", "REGEX_AA_PATTERN", "L.GG"]);
	});
});


var dtu3clPredictions = [
	{ motif: "SAVLQSGFRK"},
	{ motif: "VATVQSKMSD"},
	{ motif: "RATLQAIASE"},
	{ motif: "AVKLQNNELS"},
	{ motif: "TVRLQAGNAT"},
	{ motif: "EPMLQSADAQ"},
	{ motif: "HTVLQAVGAC"},
	{ motif: "VATLQAENVT"},
	{ motif: "FTRLQSLENV"},
	{ motif: "YPKLQSSQAW"}
];

var featureNames = [
	"NSP1",
	"NSP2",
	"NSP3",
	"NSP4",
	"NSP5",
	"NSP6",
	"NSP7",
	"NSP8",
	"NSP9",
	"NSP10",
	"NSP11",
	"NSP12",
	"NSP13",
	"NSP14",
];

_.each(dtu3clPredictions, function(dtu3clPrediction) {
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/ORF_1ab", function() {
		var variationName = "cov:clvg_dtu3cl_"+dtu3clPrediction.motif;
		glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidSimplePolymorphism", 
			"--labeledCodon", 1, 7097]);
		glue.inMode("variation/"+variationName, function() {
			glue.command(["set", "metatag", "SIMPLE_AA_PATTERN", dtu3clPrediction.motif]);
		});
	});
});

var dtu3clMatches;
var papainMatches;

glue.inMode("alignment/AL_GISAID_CONSTRAINED/member/cov-gisaid/EPI_ISL_402125", function() {
	dtu3clMatches = glue.tableToObjects(
			glue.command(["variation", "scan", "-r", "REF_MASTER_WUHAN_HU_1", "-f", "ORF_1ab", 
				"--whereClause", "name like 'cov:clvg_dtu3cl%'", "--showMatchesAsTable"]));
	papainMatches = glue.tableToObjects(
			glue.command(["variation", "scan", "-r", "REF_MASTER_WUHAN_HU_1", "-f", "ORF_1ab", 
				"--whereClause", "name = 'cov:clvg_papain'", "--showMatchesAsTable"]));
});

var allMatches = papainMatches.concat(dtu3clMatches);

var featureStartCodon = 1;
var featureEndCodon = 1;

var featureNameIndex = 0;

_.each(allMatches, function(match) {
	if(match.variationName.indexOf("papain") >= 0) {
		// for papain-like protease, cleavage happens after the diglycine residues (GG)
		// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1316023
		featureEndCodon = parseInt(match.lastRefCodon);
	} else if(match.variationName.indexOf("dtu3cl") >= 0) {
		// for 3CL protease, cleavage happens just after the Glutamine (Q) at position 5 in the motif.
		featureEndCodon = parseInt(match.firstRefCodon)+4;
	}
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/ORF_1ab", function() {
		
	});	
	glue.logInfo("featureLoc", {
		featureName: featureNames[featureNameIndex],
		startCodon:featureStartCodon,
		endCodon:featureEndCodon,
		proteinSize:(featureEndCodon - featureStartCodon) + 1
	});
	featureNameIndex++;
	featureStartCodon = featureEndCodon+1;
});
