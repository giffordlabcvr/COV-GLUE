

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", 
	"source.name = 'cov-gisaid' and pang_representative = false and pang_lineage != null and cov_glue_lineage != null",
	"sequenceID",
	"isolate",
	"pang_lineage",
	"cov_glue_lineage",
	"cov_glue_lw_ratio"]));

glue.logInfo("seqObjs", seqObjs);

var totalSeqs = seqObjs.length;
var recallFailures = 0;
var recallFailures_CGL_parent_PL = 0;
var recallFailures_PL_parent_CGL = 0;
var recallFailures_CGL_ancestor_PL = 0;
var recallFailures_PL_ancestor_CGL = 0;
var recallFailures_CGL_or_PL_putative = 0;
var recallFailures_CGL_or_PL_B_1_p16 = 0;

_.each(seqObjs, function(seqObj) {
	var pang_lineage = seqObj["pang_lineage"];
	var cov_glue_lineage = seqObj["cov_glue_lineage"];
	if(pang_lineage != cov_glue_lineage) {
		recallFailures++;
		if(isAncestorLineage(pang_lineage, cov_glue_lineage)) {
			recallFailures_PL_ancestor_CGL++;
		}
		if(isAncestorLineage(cov_glue_lineage, pang_lineage)) {
			recallFailures_CGL_ancestor_PL++;
		}
		if(isParentLineage(pang_lineage, cov_glue_lineage)) {
			recallFailures_PL_parent_CGL++;
		}
		if(isParentLineage(cov_glue_lineage, pang_lineage)) {
			recallFailures_CGL_parent_PL++;
		}
		if(pang_lineage.indexOf("p") >= 0 || cov_glue_lineage.indexOf("p") >= 0) {
			recallFailures_CGL_or_PL_putative++;
		}
		if(pang_lineage == "B.1.p16" || cov_glue_lineage == "B.1.p16") {
			recallFailures_CGL_or_PL_B_1_p16++;
		}
	}
});

glue.logInfo("results", {
	totalSeqs : totalSeqs,
	recallFailures : recallFailures,
	recallFailures_CGL_parent_PL : recallFailures_CGL_parent_PL,
	recallFailures_PL_parent_CGL : recallFailures_PL_parent_CGL,
	recallFailures_CGL_ancestor_PL : recallFailures_CGL_ancestor_PL,
	recallFailures_PL_ancestor_CGL : recallFailures_PL_ancestor_CGL,
	recallFailures_CGL_or_PL_putative : recallFailures_CGL_or_PL_putative,
	recallFailures_CGL_or_PL_B_1_p16 : recallFailures_CGL_or_PL_B_1_p16
});

//return true if lineage1 is an ancestor of lineage2 (lineages are their own ancestors)
function isAncestorLineage(lineage1, lineage2) {
	if(lineage1 == "SARS-CoV-2") {
		return true;
	}
	var lineage1Bits = lineage1.split(".")
	var lineage2Bits = lineage2.split(".")
	if(lineage1Bits.length <= lineage2Bits.length) {
		for(var i = 0; i < lineage1Bits.length; i++) {
			if(lineage1Bits[i] != lineage2Bits[i]) {
				return false;
			}
		}
		return true;
	}
	return false;
}

function isParentLineage(lineage1, lineage2) {
	if(lineage1 == "SARS-CoV-2" && (lineage2 == "A" || lineage2 == "B")) {
		return true;
	}
	var lineage1Bits = lineage1.split(".")
	var lineage2Bits = lineage2.split(".")
	if(lineage1Bits.length == (lineage2Bits.length-1)) {
		for(var i = 0; i < lineage1Bits.length; i++) {
			if(lineage1Bits[i] != lineage2Bits[i]) {
				return false;
			}
		}
		return true;
	}
	return false;
}