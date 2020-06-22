var correctionsJsonString = glue.command(["file-util", "load-string", "glue/covGisaidCharacterCorrectionRules.json"]).fileUtilLoadStringResult.loadedString;

var correctionsMap = JSON.parse(correctionsJsonString);
var seqStatsMap = {};

var fastaFiles = glue.getTableColumn(glue.command(["file-util", "list-files", "--directory", "sequences/single_fastas/"]), "fileName");

var totalNumFiles = fastaFiles.length;
var filesProcessed = 0;

_.each(fastaFiles, function(fastaFile) {
	var fastaString = glue.command(["file-util", "load-string", "sequences/single_fastas/"+fastaFile]).fileUtilLoadStringResult.loadedString;
	var fastaLines = fastaString.split(/\r\n|\r|\n/g);

	if(filesProcessed % 500 == 0) {
		glue.logInfo("Character correction / stats applied to "+filesProcessed+" of "+totalNumFiles);
	}
	filesProcessed++;
	
	var cleanedUpLines = [];

	var seqID = null;
	var seqStats = null;
	
	var ntPos = 1;

	_.each(fastaLines, function(fastaLine) {
		var trimmedLine = fastaLine.trim();
		if(trimmedLine.length == 0) {
			return;
		} else if(trimmedLine[0] == '>') {
			seqID = trimmedLine.substring(1);
			seqStats = {
					num_ns: 0,
					longest_n_run: 0,
					num_bivalent_ambigs: 0,
					num_trivalent_ambigs: 0,
					num_hyphens: 0,
					num_whitespace: 0,
					num_other: 0,
					uncorrected_illegals: false,
					illegalFastaChars: []
			};
			seqStatsMap[seqID] = seqStats;
		} else {
			var cleanedUpChars = [];
			var currentNRunLength;
			for(var i = 0; i < trimmedLine.length; i++) {
				var char = trimmedLine[i];
				var cleanedUpChar = null;
				if("nN".indexOf(char) >= 0) {
					seqStats.num_ns++;
					currentNRunLength++;
					if(currentNRunLength > seqStats.longest_n_run) {
						seqStats.longest_n_run = currentNRunLength;
					}
					cleanedUpChar = char.toUpperCase();
				} else {
					currentNRunLength = 0;
					if("uU".indexOf(char) >= 0) {
						cleanedUpChar = "T";
					} else if("acgtACGT".indexOf(char) >= 0) {
						cleanedUpChar = char.toUpperCase();
					} else if("rykmswRYKMSW".indexOf(char) >= 0) {
						seqStats.num_bivalent_ambigs++;
						cleanedUpChar = char.toUpperCase();
					} else if("bdhvBDHV".indexOf(char) >= 0) {
						seqStats.num_trivalent_ambigs++;
						cleanedUpChar = char.toUpperCase();
					} else {
						// unsual character in some sense.
						var type;
						if(char == "-") {
							seqStats.num_hyphens++;
							type = "hyphen";
						} else if(/\s/.test(char)) {
							seqStats.num_whitespace++;
							type = "whitespace";
						} else {
							seqStats.num_other++;
							type = "other";
						}
						var asciiCode = char.charCodeAt(0);
						// attempt to apply corrections from rules
						var correctionApplied = false;
						var sequenceCorrections = correctionsMap[seqID];
						if(sequenceCorrections != null) {
							var asciiCodeCorrection = sequenceCorrections[asciiCode];
							if(asciiCodeCorrection != null && 
									asciiCodeCorrection.position == "any" && 
									asciiCodeCorrection.correction == "delete") {
								correctionApplied = true;
							} 
						}
						// if no correction applied, unusual character will be skipped
						// if not hyphen, the 
						// uncorrected_illegal flag will be set to true.
						// we also make record of illegal character which could be output somewhere else.
						if(type != "hyphen" && !correctionApplied) {
							seqStats.uncorrected_illegals = true;
							var illegalFastaChar = _.find(seqStats.illegalFastaChars, function(ifc) {
								return ifc.asciiCode == asciiCode;
							});
							if(illegalFastaChar == null) {
								illegalFastaChar = {
									type: type,
									character: char,
									asciiCode: asciiCode,
									ntPositions: []
								};
								seqStats.illegalFastaChars.push(illegalFastaChar);
							}
							illegalFastaChar.ntPositions.push(ntPos);
						}
					}				
				}
				if(cleanedUpChar != null) {
					cleanedUpChars.push(cleanedUpChar);
				}
				ntPos++;
			}
			var cleanedUpLine = cleanedUpChars.join("");
			cleanedUpLines.push(cleanedUpLine);
		}
	});	
	
	var cleanedUpFastaString = ">"+seqID+"\n"+cleanedUpLines.join("\n")+"\n";
	glue.command(["file-util", "save-string", cleanedUpFastaString, "sources/cov-gisaid/"+seqID+".fasta"]);
});

var seqStatsLines = ["sequenceID\tnum_ns\tlongest_n_run\tnum_bivalent_ambigs\tnum_trivalent_ambigs\tnum_hyphens\tnum_whitespace\tnum_other\tuncorrected_illegals\n"];

_.each(_.pairs(seqStatsMap), function(pair) {
	var sequenceID = pair[0];
	var seqStats = pair[1];
	var line = sequenceID+"\t"+
		seqStats.num_ns+"\t"+
		seqStats.longest_n_run+"\t"+
		seqStats.num_bivalent_ambigs+"\t"+
		seqStats.num_trivalent_ambigs+"\t"+
		seqStats.num_hyphens+"\t"+
		seqStats.num_whitespace+
		"\t"+seqStats.num_other+
		"\t"+seqStats.uncorrected_illegals+
		"\n";
	seqStatsLines.push(line);
});

glue.command(["file-util", "save-string", seqStatsLines.join(""), "tabular/gisaidSeqCharacterStats.txt"]);

