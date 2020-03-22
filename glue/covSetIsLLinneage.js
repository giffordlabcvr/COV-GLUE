glue.command(["multi-unset", "field", "sequence", "-a", "is_l_lineage"]);

var proteinAlmt;

glue.inMode("module/covFastaProteinAlignmentExporterSeqIdOnly", function() {
	proteinAlmt = glue.command(["export", "AL_GISAID_UNCONSTRAINED", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", "ORF_8", "-l", "84", "84", 
			"-w", "sequence.include_in_ref_tree = true", 
			"-p"]);
});

var seqIdToOrf8Codon84Aa = {};

_.each(proteinAlmt.aminoAcidFasta.sequences, function(membObj) {
	seqIdToOrf8Codon84Aa[membObj.id] = membObj.sequence;
});

// look at the sequence plus its 10 nearest neighbours in the tree, what is happening in ORF 8 
// take a majority vote. If equivocal then throw an error.
var seqIdToNeighbours;

glue.inMode("module/covPhyloUtility", function() {
	var resultRows = glue.tableToObjects(glue.command(["alignment-phylogeny", "list", "neighbours", "AL_GISAID_UNCONSTRAINED", "phylogeny", 
		"-n", "10"]));
	seqIdToNeighbours = _.groupBy(resultRows, function(resultRow) { return resultRow["startSequenceID"];});
});

_.each(_.pairs(seqIdToNeighbours), function(pair) {
	var isLLineage;
	var numL = 0;
	var numS = 0;
	var seqID = pair[0];
	if(seqIdToOrf8Codon84Aa[seqID] == "L") {
		numL++;
	} else if(seqIdToOrf8Codon84Aa[seqID] == "S") {
		numS++;
	}
	_.each(pair[1], function(resultRow) {
		if(seqIdToOrf8Codon84Aa[resultRow.neighbourSequenceID] == "L") {
			numL++;
		} else if(seqIdToOrf8Codon84Aa[resultRow.neighbourSequenceID] == "S") {
			numS++;
		}
	});
	if(numL == numS) {
		throw new Error("Equal number of L ("+numL+") and S ("+numS+") in phylogenetic neighbourhood.");
	}
	isLLineage = numL > numS;
	
	glue.inMode("sequence/cov-gisaid/"+seqID, function() {
		glue.command(["set", "field", "is_l_lineage", isLLineage]);
	});
});