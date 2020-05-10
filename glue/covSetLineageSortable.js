

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", 
	"-w", "cov_glue_lineage != null", 
	"source.name", "sequenceID", "cov_glue_lineage"]));

var processed = 0;

_.each(seqObjs, function(seqObj) {
	var lineageBits = seqObj["cov_glue_lineage"].split(".");
	var lineageSortBits = [];
	_.each(lineageBits, function(lineageBit) {
		if(lineageBit.match(/^\d+$/g)) {
			lineageSortBits.push(pad(lineageBit, 3));
		} else {
			lineageSortBits.push(lineageBit);
		}
	});
	
	glue.inMode("sequence/"+seqObj["source.name"]+"/"+seqObj["sequenceID"], function() {
		glue.command(["set", "field", "--noCommit", "cov_glue_lineage_sortable", lineageSortBits.join("_")]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set sortable lineage for "+processed+"/"+seqObjs.length+" sequences. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});

glue.logInfo("Set sortable lineage for "+processed+"/"+seqObjs.length+" sequences. ");
glue.command(["commit"]);
glue.command(["new-context"]);

function pad(num, size) {
    var s = num+"";
    while (s.length < size) {
    	s = "0" + s;
    }
    return s;
}
