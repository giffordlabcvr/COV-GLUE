glue.command(["multi-unset", "link-target", "variation", "cov_insertion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_insertion_sequence", "cov_insertion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_insertion_sequence", "sequence", "-a"]);

glue.command(["multi-unset", "link-target", "cov_insertion", "cov_nt_insertion", "-a"]);
glue.command(["multi-unset", "link-target", "variation", "cov_nt_insertion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_nt_insertion_sequence", "cov_nt_insertion", "-a"]);
glue.command(["multi-unset", "link-target", "cov_nt_insertion_sequence", "sequence", "-a"]);

glue.command(["multi-delete", "cov_insertion", "-a"]);
glue.command(["multi-delete", "cov_insertion_sequence", "-a"]);

glue.command(["multi-delete", "cov_nt_insertion", "-a"]);
glue.command(["multi-delete", "cov_nt_insertion_sequence", "-a"]);


glue.command(["multi-delete", "variation", "-w", "name like 'cov_aa_ins%'"]);
glue.command(["multi-delete", "variation", "-w", "name like 'cov_nt_ins%'"]);
