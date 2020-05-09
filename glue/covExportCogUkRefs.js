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
