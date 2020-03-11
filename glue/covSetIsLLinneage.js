glue.command(["multi-unset", "field", "sequence", "-a", "is_l_lineage"]);

var seqIdToIsLLineage = {};
var proteinAlmt;

glue.inMode("module/covFastaProteinAlignmentExporterSeqIdOnly", function() {
	proteinAlmt = glue.command(["export", "AL_GISAID_UNCONSTRAINED", "-r", "REF_MASTER_WUHAN_HU_1", "-f", "ORF_8", "-l", "84", "84", "-a", "-p"]);
});

_.each(proteinAlmt.aminoAcidFasta.sequences, function(membObj) {
	if(membObj.sequence == "L") {
		seqIdToIsLLineage[membObj.id] = true;
	} else if(membObj.sequence == "S") {
		seqIdToIsLLineage[membObj.id] = false;
	} else if(membObj.id = "EPI_ISL_404253") {
		// actually a mix of S and L.
		seqIdToIsLLineage[membObj.id] = true;
	} else {
		glue.log("SEVERE", "Unexpected residue at ORF 8 position 84", membObj);
		throw new Error("Unexpected residue at ORF 8 position 84: sequenceID "+membObj.id+", residue: "+membObj.sequence);
	}
});

_.each(_.pairs(seqIdToIsLLineage), function(pair) {
	glue.inMode("sequence/cov-gisaid/"+pair[0], function() {
		glue.command(["set", "field", "is_l_lineage", pair[1]]);
	});
});