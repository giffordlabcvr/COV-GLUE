glue.command(["multi-unset", "link-target", "cov_primer_probe", "-a", "cov_primer_probe_assay"]);
glue.command(["multi-unset", "link-target", "cov_primer_probe", "-a", "seq_match"]);
glue.command(["multi-delete", "cov_primer_probe", "-a"]);
glue.command(["multi-delete", "cov_primer_probe_assay", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like '%cov_pp_seq_match_anywhere%'"]);

var ppObjs;

var ntCharToSubChars = {
    "A": ["A"],
	"C": ["C"],
	"G": ["G"],
	"T": ["T"],
	"R": ["A", "G"],
	"Y": ["C", "T"],
	"K": ["G", "T"],
	"M": ["A", "C"],
	"S": ["C", "G"],
	"W": ["A", "T"],
	"B": ["C", "G", "T"],
	"D": ["A", "G", "T"],
	"H": ["A", "C", "T"],
	"V": ["A", "C", "G"],
	"N": ["A", "C", "G", "T"]
}

var ntCharToRegex = {};

_.each(_.pairs(ntCharToSubChars), function(pair) {
	var ntChar = pair[0];
	var subChars = pair[1];
	var regex;
	if(subChars.length == 1) {
		regex = ntChar;
	} else {
		var fullSubChars = subChars.slice(); // copy array
		_.each(_.pairs(ntCharToSubChars), function(pair2) {
			if(_.difference(pair2[1], subChars).length == 0 && fullSubChars.indexOf(pair2[0]) < 0) {
				fullSubChars.push(pair2[0]);
			}
		});
		regex = "["+fullSubChars.join("")+"]";
	}
	ntCharToRegex[ntChar] = regex;
});

glue.inMode("module/tabularUtilityTab", function() {
	ppObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/faria_probes_primers/Primers_probes.txt"]));
});


_.each(ppObjs, function(ppObj) {
	var assayID = ppObj["Assay"].trim().replace(" ", "_");
	if(assayID == "See_note") {return;}
 	var createResult = glue.command(["create", "custom-table-row", 
 		"--allowExisting", "cov_primer_probe_assay", assayID]);
 	if(createResult.okResult.number != 0) {
 		glue.inMode("custom-table-row/cov_primer_probe_assay/"+assayID, function() {
 			glue.command(["set", "field", "display_name", ppObj["Assay"].trim()]);
 			var organisation = ppObj["Organisation"];
 			if(assayID.indexOf("HKU") == 0) {
 				organisation = organisation + " (HKU)";
 			}
 			glue.command(["set", "field", "organisation", organisation.trim()]);
 			glue.command(["set", "field", "url", ppObj["Reference"].trim()]);
 		}); 	 	
 	}
 	var ppID = ppObj["Name_primer_or_probe"].trim().replace(" ", "_");
 	glue.logInfo("Loading primer/probe "+ppID);
 	glue.command(["create", "custom-table-row", "cov_primer_probe", ppID]);
	var rawSequence = ppObj["Sequence"].trim();
	var processedSequenceFwd = processRawSequence(rawSequence);
	var processedSequenceRev;
	glue.inMode("module/covFastaUtility", function() {
		processedSequenceRev = glue.command(["reverse-complement", "string", "--fastaString", processedSequenceFwd]).reverseComplementFastaResult.reverseComplement;
	});
	
	var sequenceFwdRegex = sequenceToRegex(processedSequenceFwd);
	var sequenceRevRegex = sequenceToRegex(processedSequenceRev);
	
 	glue.inMode("custom-table-row/cov_primer_probe/"+ppID, function() {
 		glue.command(["set", "field", "display_name", ppObj["Name_primer_or_probe"].trim()]);
 		glue.command(["set", "field", "sequence_fwd", processedSequenceFwd]);
 		glue.command(["set", "field", "sequence_fwd_regex", sequenceFwdRegex]);
 		glue.command(["set", "field", "sequence_rev", processedSequenceRev]);
 		glue.command(["set", "field", "sequence_rev_regex", sequenceRevRegex]);
 		glue.command(["set", "link-target", "cov_primer_probe_assay", "custom-table-row/cov_primer_probe_assay/"+assayID]);
 	});
 	var fwdVariationName = "cov_pp_seq_match_anywhere_fwd:"+ppID;
 	var revVariationName = "cov_pp_seq_match_anywhere_rev:"+ppID;
	glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/whole_genome", function() {
		glue.command(["create", "variation", fwdVariationName, "-t", "nucleotideRegexPolymorphism", "--nucleotide", 1, 29903]);
		glue.inMode("variation/"+fwdVariationName, function() {
			glue.command(["set", "metatag", "REGEX_NT_PATTERN", sequenceFwdRegex]);
		});
		glue.command(["create", "variation", revVariationName, "-t", "nucleotideRegexPolymorphism", "--nucleotide", 1, 29903]);
		glue.inMode("variation/"+revVariationName, function() {
			glue.command(["set", "metatag", "REGEX_NT_PATTERN", sequenceToRegex(processedSequenceRev)]);
		});
	});
	var fwdHitObjs;
	var revHitObjs;
	glue.inMode("alignment/AL_GISAID_UNCONSTRAINED/member/cov-gisaid/EPI_ISL_402125", function() {
		fwdHitObjs = glue.tableToObjects(glue.command(["variation", "scan", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", "whole_genome", 
			"--whereClause", "name = '"+fwdVariationName+"'", "--showMatchesAsTable"]));
		revHitObjs = glue.tableToObjects(glue.command(["variation", "scan", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", "whole_genome", 
			"--whereClause", "name = '"+revVariationName+"'", "--showMatchesAsTable"]));
	});
 	glue.inMode("custom-table-row/cov_primer_probe/"+ppID, function() {
 		var refHits = fwdHitObjs.length + revHitObjs.length;
 		glue.command(["set", "field", "ref_hits", refHits]);
 		if(refHits == 1) {
 			if(fwdHitObjs.length == 1) {
 	 	 		glue.command(["set", "field", "sequence_to_scan", processedSequenceFwd]);
 	 	 		glue.command(["set", "field", "fwd_orientation", true]);
 	 	 		glue.command(["set", "field", "ref_start", fwdHitObjs[0].queryNtStart]);
 	 	 		glue.command(["set", "field", "ref_end", fwdHitObjs[0].queryNtEnd]);
 	 	 		glue.command(["set", "field", "length", processedSequenceFwd.length]);
 			} else {
 	 	 		glue.command(["set", "field", "sequence_to_scan", processedSequenceRev]);
 	 	 		glue.command(["set", "field", "fwd_orientation", false]);
 	 	 		glue.command(["set", "field", "ref_start", revHitObjs[0].queryNtStart]);
 	 	 		glue.command(["set", "field", "ref_end", revHitObjs[0].queryNtEnd]);
 	 	 		glue.command(["set", "field", "length", processedSequenceRev.length]);
 			}
 			
 		}
 	});

});

/*
 
 wuhan-hu-1 for testing
 
TCAGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCGTGCGGCACAGGCACTAGTACTGATGTCGTATACAGGGCTTTTGACATCTACAATGATAAAGTAGCTGGTTTTGCTAAATTCCTAAAAACTAATTGTTGTCGCTTCCAAGAAAAGGACGAAGATGACAATTTAATTGATTCTTACTTTGTAGTTAAGAGACACACTTTCTCTAACTACCAACATGAAGAAACAATTTATAATTTACTTAAGGATTGTCCAGCTGTTGCTAAACATGACTTCTTTAAGTTTAGAATAGACGGTGACATGGTACCACATATATCACGTCAACGTCTTACTAAATACACAATGGCAGACCTCGTCTATGCTTTAAGGCATTTTGATGAAGGTAATTGTGACACATTAAAAGAAATACTTGTCACATACAATTGTTGTGATGATGATTATTTCAATAAAAAGGACTGGTATGATTTTGTAGAAAACCCAGATATATTACGCGTATACGCCAACTTAGGTGAACGTGTACGCCAAGCTTTGTTAAAAACAGTACAATTCTGTGATGCCATGCGAAATGCTGGTATTGTTGGTGTACTGACATTAGATAATCAAGATCTCAATGGTAACTGGTATGATTTCGGTGATTTCATACAAACCACGCCAGGTAGTGGAGTTCCTGTTGTAGATTCTTATTATTCATTGTTAATGCCTATATTAACCTTGACCAGGGCTTTAACTGCAGAGTCACATGTTGACACTGACTTAACAAAGCCTTACATTAAGTGGGATTTGTTAAAATATGACTTCACGGAAGAGAGGTTAAAACTCTTTGACCGTTATTTTAAATATTGGGATCAGACATACCACCCAAATTGTGTTAACTGTTTGGATGACAGATGCATTCTGCATTGTGCAAACTTTAATGTTTTATTCTCTACAGTGTTCCCACCTACAAGTTTTGGACCACTAGTGAGAAAAATATTTGTTGATGGTGTTCCATTTGTAGTTTCAACTGGATACCACTTCAGAGAGCTAGGTGTTGTACATAATCAGGATGTAAACTTACATAGCTCTAGACTTAGTTTTAAGGAATTACTTGTGTATGCTGCTGACCCTGCTATGCACGCTGCTTCTGGTAATCTATTACTAGATAAACGCACTACGTGCTTTTCAGTAGCTGCACTTACTAACAATGTTGCTTTTCAAACTGTCAAACCCGGTAATTTTAACAAAGACTTCTATGACTTTGCTGTGTCTAAGGGTTTCTTTAAGGAAGGAAGTTCTGTTGAATTAAAACACTTCTTCTTTGCTCAGGATGGTAATGCTGCTATCAGCGATTATGACTACTATCGTTATAATCTACCAACAATGTGTGATATCAGACAACTACTATTTGTAGTTGAAGTTGTTGATAAGTACTTTGATTGTTACGATGGTGGCTGTATTAATGCTAACCAAGTCATCGTCAACAACCTAGACAAATCAGCTGGTTTTCCATTTAATAAATGGGGTAAGGCTAGACTTTATTATGATTCAATGAGTTATGAGGATCAAGATGCACTTTTCGCATATACAAAACGTAATGTCATCCCTACTATAACTCAAATGAATCTTAAGTATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTGGCACAACATGTTAAAAACTGTTTATAGTGATGTAGAAAACCCTCACCTTATGGGTTGGGATTATCCTAAATGTGATAGAGCCATGCCTAACATGCTTAGAATTATGGCCTCACTTGTTCTTGCTCGCAAACATACAACGTGTTGTAGCTTGTCACACCGTTTCTATAGATTAGCTAATGAGTGTGCTCAAGTATTGAGTGAAATGGTCATGTGTGGCGGTTCACTATATGTTAAACCAGGTGGAACCTCATCAGGAGATGCCACAACTGCTTATGCTAATAGTGTTTTTAACATTTGTCAAGCTGTCACGGCCAATGTTAATGCACTTTTATCTACTGATGGTAACAAAATTGCCGATAAGTATGTCCGCAATTTACAACACAGACTTTATGAGTGTCTCTATAGAAATAGAGATGTTGACACAGACTTTGTGAATGAGTTTTACGCATATTTGCGTAAACATTTCTCAATGATGATACTCTCTGACGATGCTGTTGTGTGTTTCAATAGCACTTATGCATCTCAAGGTCTAGTGGCTAGCATAAAGAACTTTAAGTCAGTTCTTTATTATCAAAACAATGTTTTTATGTCTGAAGCAAAATGTTGGACTGAGACTGACCTTACTAAAGGACCTCATGAATTTTGCTCTCAACATACAATGCTAGTTAAACAGGGTGATGATTATGTGTACCTTCCTTACCCAGATCCATCAAGAATCCTAGGGGCCGGCTGTTTTGTAGATGATATCGTAAAAACAGATGGTACACTTATGATTGAACGGTTCGTGTCTTTAGCTATAGATGCTTACCCACTTACTAAACATCCTAATCAGGAGTATGCTGATGTCTTTCATTTGTACTTACAATACATAAGAAAGCTACATGATGAGTTAACAGGACACATGTTAGACATGTATTCTGTTATGCTTACTAATGATAACACTTCAAGGTATTGGGAACCTGAGTTTTATGAGGCTATGTACACACCGCATACAGTCTTACAG

*/

function sequenceToRegex(sequence) {
	var seqRegex = "";
	for(var i = 0; i < sequence.length; i++) {
		seqRegex += ntCharToRegex[sequence[i]];
	}
	return seqRegex;
}

function processRawSequence(rawSequence) {
	var processedSequence = rawSequence.replaceAll(" ", "");
	processedSequence = processedSequence.replaceAll("6FAM-", "");
	processedSequence = processedSequence.replaceAll("FAM-", "");
	processedSequence = processedSequence.replaceAll("-BBQ", "");
	processedSequence = processedSequence.replaceAll("-BHQ1", "");
	processedSequence = processedSequence.replaceAll("-BQH1", "");
	processedSequence = processedSequence.replaceAll("-BHQ", "");
	processedSequence = processedSequence.replaceAll("-TAMRA", "");
	for(i = 0; i < processedSequence.length; i++) {
		if("ACGTRYKMSWBDHVN".indexOf(processedSequence[i]) < 0) {
			throw new Error("Unrecognised sequence char '"+processedSequence[i]+"'");
		}
	}
	return processedSequence;
}
