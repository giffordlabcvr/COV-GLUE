

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", 
//	"-w", "source.name = 'cov-gisaid'", 
//	"-w", "sequenceID = 'EPI_ISL_417149'", 
	"-w", "gisaid_authors like '%etl al%'", 
	"source.name", "sequenceID", "m49_country.id", "gisaid_authors"]));

var processed = 0;

var shortAuthorVals = {};

_.each(seqObjs, function(seqObj) {
	var longAuthors = seqObj["gisaid_authors"];
	if(longAuthors == null) {
		return;
	}
	longAuthors = longAuthors.trim();
	longAuthors = longAuthors.replace("，", ", "); // some wierd character that combines a comman with a space.
	longAuthors = longAuthors.replace("#", ""); 
	longAuthors = longAuthors.replace("*", ""); 
	longAuthors = longAuthors.replace(" & ", ", "); 
	var shortAuthors;
	if(longAuthors.indexOf("Wellington SCL") >= 0) {
		shortAuthors = "Wellington SCL";
	} else if(longAuthors.indexOf("DCLS") >= 0 || longAuthors.toLowerCase().indexOf("division of consolidated") >= 0) {
		shortAuthors = "Virginia DCLS";
	} else if(longAuthors.endsWith("et al") > 0 && longAuthors.length < 30) { 
		shortAuthors = longAuthors+"."; // Chu et al
	} else if(longAuthors.endsWith("et al.") > 0 && longAuthors.length < 30) { 
		shortAuthors = longAuthors; 
	} else if(longAuthors.indexOf("etl al") > 0 && longAuthors.length < 30) { 
		shortAuthors = longAuthors.replace("etl al", "et al.");// Chu etl al
	} else if(longAuthors.indexOf("; ") > 0) {
		var first = longAuthors.split("; ")[0];
		if(first.indexOf(", ") > 0) {
			first = first.split(", ")[0]; 
		}
		if(first.indexOf(",") > 0) {
			first = first.split(",")[0]; 
		}
		if(first.indexOf(" ") > 0) {
			var bits = first.split(" ");
			var lastBit = bits[bits.length-1];
			shortAuthors = lastBit+" et al.";  // mix of commas and semicolons separating names 
		} else {
			shortAuthors = first+" et al."; // names separated by semicolon, surname given first
		}
	} else if(longAuthors.indexOf(", ") > 0) {
		var first = longAuthors.split(", ")[0];
		if(first.indexOf(",") > 0) {
			shortAuthors = first.split(",")[0]+" et al."; // surname given first, then comma with no spaces and initials
		} else {
			var bits = first.split(" ");
			var lastBit = bits[bits.length-1];
			if(lastBit.replace(/\./g, "").replace(/-/g, "").length <= 2 && lastBit.toUpperCase() == lastBit) {
				shortAuthors = bits[0]+" et al."; // e.g. Smith JB, Bloggs F
			} else {
				shortAuthors = lastBit+" et al.";
			}
		}
	} else if(longAuthors.indexOf(",") > 0) {
		var bits = longAuthors.split(",");
		shortAuthors = bits[0]+" et al.";
	} else if(longAuthors.indexOf(" ") > 0) {
		shortAuthors = longAuthors;
	} else if(longAuthors.toUpperCase() == longAuthors) {
		shortAuthors = longAuthors+" et al.";
	}
	if(shortAuthors == null || shortAuthors.trim().length == 0) {
		throw new Error("Not set for: "+longAuthors);
	}
	if(shortAuthors.match(/^[A-Z]+ et al\.$/g)) {
		shortAuthors = shortAuthors[0] + shortAuthors.toLowerCase().substring(1);
	}
	if(shortAuthors.match(/^[a-z]+ et al\.$/g)) {
		shortAuthors = shortAuthors[0].toUpperCase() + shortAuthors.substring(1);
	}
	
	//glue.logInfo(seqObj["sequenceID"], shortAuthors);
	shortAuthorVals[shortAuthors] = seqObj["sequenceID"];
	glue.inMode("sequence/"+seqObj["source.name"]+"/"+seqObj["sequenceID"], function() {
		glue.command(["set", "field", "--noCommit", "gisaid_authors_short", shortAuthors]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set short authors for "+processed+"/"+seqObjs.length+" sequences. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});

glue.logInfo("Set short authors for "+processed+"/"+seqObjs.length+" sequences. ");
glue.command(["commit"]);
glue.command(["new-context"]);

//glue.logInfo("shortAuthorVals", shortAuthorVals);
