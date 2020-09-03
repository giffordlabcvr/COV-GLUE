glue.command(["multi-unset", "link-target", "variation", "cov_replacement", "-a"]);
glue.command(["multi-unset", "link-target", "cov_replacement_sequence", "cov_replacement", "-a"]);
glue.command(["multi-unset", "link-target", "cov_replacement_sequence", "sequence", "-a"]);

glue.command(["multi-delete", "cov_replacement", "-a"]);
glue.command(["multi-delete", "cov_replacement_sequence", "-a"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_rpl%'"]);
