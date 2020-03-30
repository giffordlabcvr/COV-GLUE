regex=".*(EPI_ISL_[0-9]+).*"
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
    	if [[ $line =~ $regex ]]
        then
	        id="${BASH_REMATCH[1]}"
	        outfile="sequences/single_fastas/${id}.fasta"
	        echo ">${id}" > $outfile
	    else
	    	echo Bad header $line
	    	exit 1
        fi
    else
	    echo $line >> $outfile
    fi
done < sequences/gisaid_cov202sequences.fasta
