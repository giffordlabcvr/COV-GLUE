git pull
gluetools.sh -p log-level:FINEST -i run file covProject.glue
gluetools.sh -p log-level:FINEST -i project cov run file glue/covRegenerateGisaidAlignment.glue
gluetools.sh -p log-level:FINEST -i project cov run file glue/covGeneratePhylogeny.glue
git add alignments/AL_GISAID_UNCONSTRAINED.fna
git add trees/gisaidUnconstrainedUnrooted.tree
git add trees/gisaidMidpointRooted.tree
git commit -m "Alignment and tree update"
git push
