rm -rf sequences/single_fastas
mkdir -p sequences/single_fastas

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
	    elif [[ ${line} == '>||' ]]
	    then
			echo "Found known bad header >||"
	    else
	    	echo Bad header $line
	    	exit 1
        fi
    else
	    echo $line >> $outfile
    fi
done < sequences/gisaid_cov2020_sequences.fasta
