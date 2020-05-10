// This file is run after importing PANG lineages while the COGUK-COV-GLUE extension is loaded.
// it exports sequences and metadata back into the COV-GLUE repo so that cov-coguk-refs can 
// be loaded independently of the extension.

glue.command(["delete", "source", "cov-coguk-refs"]);
glue.command(["create", "source", "cov-coguk-refs"]);
glue.command(["copy", "sequence", "cov-coguk-refs", "-w", "source.name = 'cov-coguk' and pang_representative = true"]);

glue.command(["file-util", "make-directories", "sources/cov-coguk-refs"]);

var fileObjs = glue.tableToObjects(
		glue.command(["file-util", "list-files", "-d", "sources/cov-coguk-refs"]));

_.each(fileObjs, function(fileObj) {
	if(fileObj.fileName.indexOf(".fasta") >= 0) {
		glue.command(["file-util", "delete-file", "sources/cov-coguk-refs/"+fileObj.fileName]);
	}
});

glue.command(["export", "source", "-p", "sources", "cov-coguk-refs"]);

var cogukRefsMetadata = glue.command(
		["list", "sequence", "-w", "source.name = 'cov-coguk' and pang_representative = true",
		 "sequenceID", "isolate", "place_sampled", "m49_country.id", 
		 "collection_month_day", "collection_month", "collection_year"]);

glue.inMode("module/tabularUtilityTab", function() {
	glue.command({"save-tabular": {
		"tabularData":cogukRefsMetadata, 
		"fileName":"tabular/cogukRefsMetadata.txt"
	}});
});