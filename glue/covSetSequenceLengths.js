

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", "source.name = 'cov-gisaid'"]));

_.each(seqObjs, function(seqObj) {
	glue.inMode("sequence/"+seqObj["source.name"]+"/"+seqObj["sequenceID"], function() {
		var length = glue.command(["show", "length"]).lengthResult.length;
		glue.command(["set", "field", "length", length]);
	});
});