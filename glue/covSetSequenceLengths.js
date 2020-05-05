

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", "source.name = 'cov-gisaid'"]));

var processed = 0;

_.each(seqObjs, function(seqObj) {
	glue.inMode("sequence/"+seqObj["source.name"]+"/"+seqObj["sequenceID"], function() {
		var length = glue.command(["show", "length"]).lengthResult.length;
		glue.command(["set", "field", "--noCommit", "length", length]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set sequence length for "+processed+" sequences. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});

glue.logInfo("Set sequence length for "+processed+" sequences. ");
glue.command(["commit"]);
glue.command(["new-context"]);
