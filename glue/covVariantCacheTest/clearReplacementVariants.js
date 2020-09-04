var replacements = ['N:G:204:R', 'N:R:203:K', 'NSP12:P:323:L', 'NSP6:L:37:F', 
	'S:D:614:G', 'M:T:175:M', 'NSP1:H:83:R', 
	'NSP1:M:85:T', 'NSP1:V:84:R'];

_.each(replacements, function(repID) {
	glue.command(["multi-unset", "link-target", "variation", "cov_replacement", "-w", "name = 'cov_aa_rpl:"+repID+"'"]);
	glue.command(["multi-delete", "cov_replacement", "-w", "id = '"+repID+"'"]);
	glue.command(["multi-delete", "variation", "-w", "name = 'cov_aa_rpl:"+repID+"'"]);
});