glue.command(["multi-unset", "field", "sequence", "-a", "num_unique_snps"]);

//static map of concrete / ambiguous NT chars to the set of concrete chars which match them.
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

var nucleotideAlmt;

glue.inMode("module/covFastaAlignmentExporterSeqIdOnly", function() {
	nucleotideAlmt = glue.command(["export", "AL_GISAID_CONSTRAINED", 
			"-r", "REF_MASTER_WUHAN_HU_1", "-f", "coding_spanning_region", 
			"-a", 
			"-p"]);
});


var totalNumSeqs = nucleotideAlmt.nucleotideFasta.sequences.length;
var positionInfoArray = null;
var pass1processed = 0;

_.each(nucleotideAlmt.nucleotideFasta.sequences, function(membObj) {
	if(positionInfoArray == null) {
		positionInfoArray = [];
		for(var i = 0; i < membObj.sequence.length; i++) {
			positionInfoArray.push({});
		}
	}
	for(var i = 0; i < membObj.sequence.length; i++) {
		var membNt = membObj.sequence[i];
		if(membNt != '-' && membNt != 'N') {
			var positionInfo = positionInfoArray[i];
			if(membObj.id == "EPI_ISL_402125") {
				positionInfo.refNt = membNt;
			} else {
				var subChars = ntCharToSubChars[membNt];
				for(var j = 0; j < subChars.length; subChars++) {
					var subChar = subChars[j];
					var currentCount = positionInfo[subChar];
					if(currentCount == null) {
						positionInfo[subChar] = 1;
					} else {
						positionInfo[subChar]++;
					}
				}
			}
		}
	}
	pass1processed++;
	if(pass1processed % 500 == 0) {
		glue.logInfo("Pass 1/2 of counting unique SNPs, processed "+pass1processed+" of "+totalNumSeqs);
	}
});

var pass2processed = 0;

_.each(nucleotideAlmt.nucleotideFasta.sequences, function(membObj) {
	var num_unique_snps;
	if(membObj.id == "EPI_ISL_402125") {
		num_unique_snps = 0;
	} else {
		num_unique_snps = 0;
		for(var i = 0; i < membObj.sequence.length; i++) {
			var membNt = membObj.sequence[i];
			if(membNt != '-' && membNt != 'N') {
				var positionInfo = positionInfoArray[i];
				var subChars = ntCharToSubChars[membNt];
				var refNt = positionInfo.refNt;
				for(var j = 0; j < subChars.length; subChars++) {
					var subChar = subChars[j];
					var currentCount = positionInfo[subChar];
					if(subChar != refNt && currentCount == 1) {
						num_unique_snps++;
					}
				}
			}
		}
	}
	glue.inMode("sequence/cov-gisaid/"+membObj.id, function() {
		glue.command(["set", "field", "num_unique_snps", num_unique_snps]);
	});
	pass2processed++;
	if(pass2processed % 500 == 0) {
		glue.logInfo("Pass 2/2 of counting unique SNPs, processed "+pass2processed+" of "+totalNumSeqs);
	}

});
glue.logInfo("Counting unique SNPs complete");

