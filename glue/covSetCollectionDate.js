

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", 
	"-w", "collection_year != null and collection_month != null and collection_month_day != null", 
	"source.name", "sequenceID", "collection_year", "collection_month", "collection_month_day"]));

var processed = 0;

_.each(seqObjs, function(seqObj) {
	var collectionDateString = seqObj["collection_month_day"]+"-"+seqObj["collection_month"]+"-"+seqObj["collection_year"];
	if(seqObj["collection_month_day"] < 10) {
		collectionDateString = "0" + collectionDateString;
	}
	glue.inMode("sequence/"+seqObj["source.name"]+"/"+seqObj["sequenceID"], function() {
		glue.command(["set", "field", "--noCommit", "collection_date", collectionDateString]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set collection date for "+processed+"/"+seqObjs.length+" sequences. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});

glue.logInfo("Set collection date for "+processed+"/"+seqObjs.length+" sequences. ");
glue.command(["commit"]);
glue.command(["new-context"]);
