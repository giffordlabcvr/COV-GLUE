/*
 *  EPI_ISL_500981 contains 
 *  nt_insertion with id NSP6:nca:11074:83412:11075
 *  insertion with id NSP6:ca:34:70:35
 *  replacements with ids: 
 *  'N:G:204:R', 'N:R:203:K', 'NSP12:P:323:L', 'NSP6:L:37:F', 'S:D:614:G'
 */
	
/*  EPI_ISL_465549 contains
 *  nt_deletion with id NSP1:nca:521:523
 *  deletion with id NSP1:ca:86:86
 *  replacements with ids: 
 *  'M:T:175:M', 'N:G:204:R', 'N:R:203:K', 'NSP12:P:323:L', 'NSP1:H:83:R', 'NSP1:M:85:T', 'NSP1:V:84:R', 'S:D:614:G'
 */

/* Tests
 * 1. If variants are not in the cache, sequences are not in the cache, then variants and associations are created.
	a. Restore DB
 	b. Delete sequence and cache files
	c. Delete associations, then variants from DB
	d. Run generate scripts with appropriate sequence tweak
	e. Check that variants and associations created.
 * 2. If variants are in the cache, sequences are not in the cache, association is created.
	a. Restore DB and cache files
 	b. Delete sequence cache files
	c. Delete associations from DB
	d. Run generate scripts with appropriate sequence tweak
	e. Check that associations created.
 * 3. If both are in the cache, variants and associations are created during populate-from-cache
	a. Restore DB and cache files
	b. Delete associations, then variants from DB
	c. Run populate scripts with appropriate sequence tweak
	d. Check that variation and associations created.
 */

// this is the check step
glue.inMode("sequence/cov-gisaid/EPI_ISL_500981", function() {
	glue.logInfo("EPI_ISL_500981 replacements", 
			glue.getTableColumn(glue.command(["list", "link-target", "cov_replacement_sequence", "id"]), "id"))
	glue.logInfo("EPI_ISL_500981 NT insertions", 
			glue.tableToObjects(glue.command(["list", "link-target", "cov_nt_insertion_sequence", "id"]), "id"))
	glue.logInfo("EPI_ISL_500981 insertions", 
			glue.tableToObjects(glue.command(["list", "link-target", "cov_insertion_sequence", "id"]), "id"))
});

glue.inMode("sequence/cov-gisaid/EPI_ISL_465549", function() {
	glue.logInfo("EPI_ISL_465549 replacements", 
			glue.getTableColumn(glue.command(["list", "link-target", "cov_replacement_sequence", "id"]), "id"))
	glue.logInfo("EPI_ISL_465549 NT deletions", 
			glue.tableToObjects(glue.command(["list", "link-target", "cov_nt_deletion_sequence", "id"]), "id"))
	glue.logInfo("EPI_ISL_465549 deletions", 
			glue.tableToObjects(glue.command(["list", "link-target", "cov_deletion_sequence", "id"]), "id"))
});
