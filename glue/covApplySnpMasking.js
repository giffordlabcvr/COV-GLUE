// Apply SNPs to alignment from CSV file imported from the 
// Rambaut et al. lineages repo.

// reference coding spanning region
var refNucleotides;

// list of objects containing each sequence ID and the coding spanning region nucleotides
var memberSeqObjs;

glue.inMode("module/covFastaAlignmentExporterSourcePlusSeqId", function() {
	refNucleotides = glue.command(["export", "AL_GISAID_UNCONSTRAINED", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", "coding_spanning_region", 
			"-w", "sequence.sequenceID = 'EPI_ISL_402125'", 
			"-p"]).nucleotideFasta.sequences[0].sequence;
	memberSeqObjs = glue.command(["export", "AL_GISAID_UNCONSTRAINED", 
		"-r", "REF_MASTER_WUHAN_HU_1", "-f", "coding_spanning_region", 
		"-w", "sequence.pang_representative = true", 
		"-p"]).nucleotideFasta.sequences;
});

_.each(memberSeqObjs, function(membSeqObj) {
	var bits = membSeqObj.id.split("/");
	var sourceName = bits[0];
	var sequenceID = bits[1];
	var snpsString;
	// list the SNPs that have been recorded against the sequence.
	glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
		snpsString = glue.command(["show", "property", "pang_masked_snps"]).propertyValueResult.value;
	});
	var snps;
	if(snpsString == null) {
		snps = [];
	} else {
		snps = snpsString.split(",");
	}
	// apply the SNP masking within the unconstrained alignment that is used for
	// the reference tree generation.
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED/member/"+sourceName+"/"+sequenceID, function() {
		_.each(snps, function(snp) {
			var snpRefNtChar = snp[0];
			var snpMembNtChar = snp[snp.length-1];
			var refCoord = parseInt(snp.substring(1, snp.length-1));
			// need to apply an offset here because of the coding spanning region.
			var refNtChar = refNucleotides[refCoord-266];
			var membNtChar = membSeqObj.sequence[refCoord-266];
			// check that the SNP proposed to be masked is correct.
			if(refNtChar != snpRefNtChar) {
				glue.log("WARNING", "Sequence "+sourceName+"/"+sequenceID+": Reference character for SNP "+snp+" did not match. Expected "+snpRefNtChar+", found "+refNtChar)
			}
			if(membNtChar != snpMembNtChar) {
				glue.log("WARNING", "Sequence "+sourceName+"/"+sequenceID+": Member character for SNP "+snp+" did not match. Expected "+snpMembNtChar+", found "+membNtChar)
			}
			glue.command(["mask", "positions", "-r", "REF_MASTER_WUHAN_HU_1", refCoord, refCoord]);
			glue.log("FINEST", "Applied masking for SNP "+snp+" in sequence "+sourceName+"/"+sequenceID);
		});
	});
});
