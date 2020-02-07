git pull
gluetools.sh -p log-level:FINEST -i run file covAlmtTreeUpdate.glue
git add trees/gisaidUnconstrainedUnrooted.tree
git add trees/gisaidMidpointRooted.tree
git commit -m "Alignment and tree update"
git push
