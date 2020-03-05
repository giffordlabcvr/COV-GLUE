glue.command(["delete", "alignment", "AL_CVR_CONSTRAINED"]);

glue.command(["create", "alignment", "AL_CVR_CONSTRAINED", "-r", "REF_MASTER_WUHAN_HU_1"]);

glue.inMode("alignment/AL_CVR_CONSTRAINED", function() {
	// Wuhan-Hu-1
	glue.command(["add", "member", "cov-gisaid", "EPI_ISL_402125"]);
	
	glue.inMode("member/cov-gisaid/EPI_ISL_402125", function() {
		glue.command(["add", "segment", 1, 29903, 1, 29903]);
	});
	// Germany/BavPat1/2020
	glue.command(["add", "member", "cov-gisaid", "EPI_ISL_406862"]);
	//Mexico/CDMX/InDRE_01/2020
	glue.command(["add", "member", "cov-gisaid", "EPI_ISL_412972"]);
	// Germany/Baden-Wuerttemberg-1/2020
	glue.command(["add", "member", "cov-gisaid", "EPI_ISL_412912"]);
	// Italy/CDG1/2020
	glue.command(["add", "member", "cov-gisaid", "EPI_ISL_412973"]);
	// CVR01
	glue.command(["add", "member", "cvr", "CVR01"]);

	glue.command(["derive", "segments", "AL_GISAID_UNCONSTRAINED", "--existingMembersOnly", "-w", "sequence.sequenceID != 'EPI_ISL_402125'"]);
});

var almtRows;

glue.inMode("module/covFastaAlignmentExporterSeqIdOnly", function() {
	almtRows = glue.command(["export", "AL_CVR_CONSTRAINED", "--allMembers", "--preview"]).nucleotideFasta.sequences;
});

var comparisons = [
	{name: "Wuhan-Hu-1", id: "EPI_ISL_402125"},
	{name: "Germany/BavPat1/2020", id: "EPI_ISL_406862"},
	{name: "Mexico/CDMX/InDRE_01/2020", id: "EPI_ISL_412972"},
	{name: "Germany/Baden-Wuerttemberg-1/2020", id: "EPI_ISL_412912"},
	{name: "Italy/CDG1/2020", id: "EPI_ISL_412973"}
]

var CVR01Seq = _.find(almtRows, function(almtRow) {
	if(almtRow.id == "CVR01") {
		return true;
	}
}).sequence;

var table = "VirusName\tGISAID_ID\tReferenceNtCoord\tComparisonNt\tCVR01Nt\n";

_.each(comparisons, function(comparison) {
	var comparisonSeq = _.find(almtRows, function(almtRow) {
		if(almtRow.id == comparison.id) {
			return true;
		}
	}).sequence;
	for(var nt = 1; nt <= 29903; nt++) {
		var charCVR01 = CVR01Seq[nt-1];
		var charComp = comparisonSeq[nt-1];
		if(charCVR01 != '-' && charComp != '-' && charCVR01 != 'N' && charCVR01 != charComp) {
			glue.logInfo("COV01 difference", {
				name: comparison.name,
				id: comparison.id,
				refNt: nt,
				char: charComp,
				charCVR01: charCVR01,
			});
			table += comparison.name+"\t"
			table += comparison.id+"\t"
			table += nt+"\t"
			table += charComp+"\t"
			table += charCVR01+"\n"
		}
	}
});

glue.command(["file-util", "save-string", table, "CVR01_nt_diffs.txt"]);