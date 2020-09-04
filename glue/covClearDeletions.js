glue.command(["multi-unset", "link-target", "variation", "cov_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_deletion_sequence", "cov_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_deletion_sequence", "sequence", "-a"]);

glue.command(["multi-unset", "link-target", "cov_deletion", "cov_nt_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "variation", "cov_nt_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_nt_deletion_sequence", "cov_nt_deletion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_nt_deletion_sequence", "sequence", "-a"]);

glue.command(["multi-delete", "cov_deletion", "-a"]);
glue.command(["multi-delete", "cov_deletion_sequence", "-a"]);

glue.command(["multi-delete", "cov_nt_deletion", "-a"]);
glue.command(["multi-delete", "cov_nt_deletion_sequence", "-a"]);

glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_del%'"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_nt_del%'"]);
