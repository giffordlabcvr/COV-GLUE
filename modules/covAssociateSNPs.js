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

function associateSNPs(whereClause) {
	var memberObjs;
	glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
		memberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", whereClause]));
	});
	var refAlmtRow = alignmentRow({"sequence.source.name": "cov-gisaid", "sequence.sequenceID": "EPI_ISL_402125"});
	glue.command(["new-context"]);
	var processed = 0;
	_.each(memberObjs, function(memberObj) {
		var membAlmtRow = alignmentRow(memberObj);
		var sequenceID = memberObj["sequence.sequenceID"];
		var sourceName = memberObj["sequence.source.name"];
		for(var i = 0; i < refAlmtRow.length; i++) {
			var membChar = membAlmtRow[i];
			var refChar = refAlmtRow[i];
			var refNt = i+266;  // 266 is start of coding region on reference.
			if(membChar != 'N' && membChar != '-' && membChar != refChar) {
				var subChars = ntCharToSubChars[membChar];
				for(var j = 0; j < subChars.length; subChars++) {
					var subChar = subChars[j];
					if(subChar != refChar) {
						var snpName = refChar+refNt+subChar;
						var numCreated = glue.command(["create", "custom-table-row", "--allowExisting", "cov_nt_mutation", snpName]).okResult.number;
						if(numCreated == 1) {
							glue.inMode("custom-table-row/cov_nt_mutation/"+snpName, function() {
								glue.command(["set", "field", "display_name", snpName]);
								glue.command(["set", "field", "reference_nt", refNt]);
								glue.command(["set", "field", "reference_base", refChar]);
								glue.command(["set", "field", "mutation_base", subChar]);
							});
						}
						var assocId = snpName+":"+sequenceID+":"+sourceName;
						glue.command(["create", "custom-table-row", "cov_nt_mutation_sequence", assocId]);
						glue.inMode("custom-table-row/cov_nt_mutation_sequence/"+assocId, function() {
							glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
							glue.command(["set", "link-target", "cov_nt_mutation", "custom-table-row/cov_nt_mutation/"+snpName]);
						});
					}
				}
			}
		}
		processed++;
		if(processed % 50 == 0) {
			glue.logInfo("Associating SNPs, processed "+processed+" of "+memberObjs.length);
			glue.command(["new-context"]);
		}

	});
	glue.logInfo("Associating SNPs, processed "+processed+" of "+memberObjs.length);
	glue.command(["new-context"]);
}

function findMissingSNPs(whereClause) {
	var snpObjs = glue.tableToObjects(glue.command(["list", "custom-table-row", "cov_nt_mutation", "id", "reference_nt"]));
	var refNtToSnpObjs = _.groupBy(snpObjs, function(snpObj) {return snpObj.reference_nt;});
	
	glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
		memberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", whereClause]));
	});
	var refAlmtRow = alignmentRow({"sequence.source.name": "cov-gisaid", "sequence.sequenceID": "EPI_ISL_402125"});
	glue.command(["new-context"]);
	var processed = 0;
	_.each(memberObjs, function(memberObj) {
		var membAlmtRow = alignmentRow(memberObj);
		var sequenceID = memberObj["sequence.sequenceID"];
		var sourceName = memberObj["sequence.source.name"];
		for(var i = 0; i < refAlmtRow.length; i++) {
			var refNt = i+266;  // 266 is start of coding region on reference.
			var snpObjs = refNtToSnpObjs[refNt];
			if(snpObjs != null) { // at least one SNP in the db at this location
				var membChar = membAlmtRow[i];
				var refChar = refAlmtRow[i];
				if((membChar == 'N' || membChar == '-') && refChar != 'N' && refChar != '-') {
					_.each(snpObjs, function(snpObj) {
						var assocId = snpObj.id+":"+sequenceID+":"+sourceName;
						glue.command(["create", "custom-table-row", "cov_nt_mutation_sequence_missing", assocId]);
						glue.inMode("custom-table-row/cov_nt_mutation_sequence_missing/"+assocId, function() {
							glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
							glue.command(["set", "link-target", "cov_nt_mutation", "custom-table-row/cov_nt_mutation/"+snpObj.id]);
						});
					});
				}
			}
		}
		processed++;
		if(processed % 50 == 0) {
			glue.logInfo("Finding missing SNPs, processed "+processed+" of "+memberObjs.length);
			glue.command(["new-context"]);
		}

	});
	glue.logInfo("Finding missing SNPs, processed "+processed+" of "+memberObjs.length);
	glue.command(["new-context"]);

}

function alignmentRow(memberObj) {
	var nucleotideAlmt;

	glue.inMode("module/covFastaAlignmentExporterSeqIdOnly", function() {
		nucleotideAlmt = glue.command(["export", "AL_GISAID_CONSTRAINED", 
				"-r", "REF_MASTER_WUHAN_HU_1", "-f", "coding_spanning_region", 
				"-w", 
				"sequence.sequenceID = '"+memberObj["sequence.sequenceID"]+"'"+
				" and sequence.source.name = '"+memberObj["sequence.source.name"]+"'", 
				"-p"]);
	});
	return nucleotideAlmt.nucleotideFasta.sequences[0].sequence;
}
