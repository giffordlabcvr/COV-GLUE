var correctionsJsonString = glue.command(["file-util", "load-string", "glue/covGisaidCharacterCorrectionRules.json"]).fileUtilLoadStringResult.loadedString;

var correctionsMap = JSON.parse(correctionsJsonString);
var seqStatsMap = {};

var fastaFiles = glue.getTableColumn(glue.command(["file-util", "list-files", "--directory", "sequences/single_fastas/"]), "fileName");


_.each(fastaFiles, function(fastaFile) {
	var fastaString = glue.command(["file-util", "load-string", "sequences/single_fastas/"+fastaFile]).fileUtilLoadStringResult.loadedString;
	var fastaLines = fastaString.split(/\r\n|\r|\n/g);

	var cleanedUpLines = [];

	var seqID = null;
	var seqStats = {
			num_ns: 0,
			longest_n_run: 0,
			num_bivalent_ambigs: 0,
			num_trivalent_ambigs: 0,
			num_hyphens: 0,
			num_whitespace: 0,
			num_other: 0,
			illegalFastaChars: []
	};
	seqStatsMap[seqID] = seqStats;
	
	var ntPos = 1;

	_.each(fastaLines, function(fastaLine) {
		var trimmedLine = fastaLine.trim();
		if(trimmedLine.length == 0) {
			return;
		} else if(trimmedLine[0] == '>') {
			seqID = trimmedLine.substring(1);
			seqStats.sequenceID = seqID;
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
					if("acgtuACGTU".indexOf(char) >= 0) {
						cleanedUpChar = char.toUpperCase();
					} else if("rykmswRYKMSW".indexOf(char) >= 0) {
						seqStats.num_bivalent_ambigs++;
						cleanedUpChar = char.toUpperCase();
					} else if("bdhvBDHV".indexOf(char) >= 0) {
						seqStats.num_trivalent_ambigs++;
						cleanedUpChar = char.toUpperCase();
					} else {
						// illegal in some sense.
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
						var sequenceCorrections = correctionsMap[seqID].corrections;
						if(sequenceCorrections != null) {
							var asciiCodeCorrection = sequenceCorrections[asciiCode];
							if(asciiCodeCorrection != null && 
									asciiCodeCorrection.position == "any" && 
									asciiCodeCorrection.correction == "delete") {
								correctionApplied = true;
							} 
						}
						// if no correction applied, make record of illegal character
						if(!correctionApplied) {
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

var seqStatsLines = ["num_ns\tlongest_n_run\tnum_bivalent_ambigs\tnum_trivalent_ambigs\tnum_hyphens\tnum_whitespace\tnum_other\n"];

// seqStatsMap
