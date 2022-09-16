#!/bin/env bash

#####################################################################################################
#Script Name    : lco_xv_v0.1.sh
#Description    : Given an aligned set of sequences calculate the cross-validation log-likelihood
#                 for a single column at a time
#Author         : Michael A Sennett & Douglas Theobald
#####################################################################################################

#####################################################################################################

usage() {
	echo "Usage: $0 [ -a alignment fasta format ] [ -m model iqtree format ] \
[ -t tree, newick format ] [-N number of cores available] [-n cores for 1 iqtree run]" 1>&2
}

exit_abnormal() {
	usage
	exit 1
}

#####################################################################################################
#Input Variables:  a/ALGN  is the alignment that will be reconstructed,
#                          alignment must be in .a2m or.fasta format
#		   m/MODEL is the evolutionary model, must be in IQ-Tree format
#                  t/TREE  is a starting tree for the alignment, default is "None"
#		   N/TND   is the total number of threads available to run analysis, 
#			   so the process can be parallelized
#                  n/NDS   the number of threads for one iqtree run
#####################################################################################################

TREE="None"

while getopts ":a:m:t:h:n:N:" options; do
	case "${options}" in
		a)
			ALGN=${OPTARG}
			;;
		m)
			MODEL=${OPTARG}
			;;
		t)
			TREE=${OPTARG}
			;;
		n)
			NDS=${OPTARG}
			;;
		N)
			TND=${OPTARG}
			;;
		h)
			usage
			;;
		*)
			exit_abnormal
			;;
	esac
done


echo "$ALGN" "$TREE" "$TND" "$NDS" "$MODEL";
sleep 2;

#####################################################################################################

N=$(($TND/$NDS)); #determine numb of parallel iqtree runs
rm siteloglik.txt; #just in case this file already exists in folder

#####################################################################################################
#
#Functions required to run the LOCO-XV
#
#####################################################################################################

function do_LCO_tree()
{
    nice iqtree2 -redo -nt $NDS -s ${ALN} -st AA -m ${MODEL} -quiet;
}
function do_LCO_tree_starter()
{
    nice iqtree2 -redo -nt $NDS -s ${ALN} -t ${TREE} -st AA -m ${MODEL} -quiet;
}
function do_wsl()
{
    PRM=${ALN}.iqtree
    echo "$PRM"

    IFS='+' read -a ARR <<< "$MODEL"; #creates an array containing  each model param

    SM=${ARR[0]};

    REC_MDL=$(echo $SM);

    for p in ${ARR[@]}; do
	if [ $p == 'FO' ] || [ $p == 'F' ]; then
		FREQ=$(egrep -A 21 "State frequencies:" ${PRM} | tail -20 | awk '{print $3}'| paste -s -d ',')
		FR=$(echo "F{"$FREQ"}")
		REC_MDL="$REC_MDL+$FR";
		break
	elif [ $p == 'FQ' ]; then
		REC_MDL="$REC_MDL+FQ";
		break
	fi
    done

    for p in ${ARR[@]}; do
	if [ $p == 'I' ]; then
		INV=$p;
		PROP=$(egrep "Proportion of invariable sites:" ${PRM} | sed 's/[^0-9.]*//g');
		echo "$PROP";
		REC_MDL="$REC_MDL+$INV";
		break
	else
		PROP=-1
	fi
    done

    for p in ${ARR[@]}; do
	if [[ $p == G* ]]; then
		GAM=$p;
		ALPHA=$(egrep "alpha" ${PRM} | sed 's/[^0-9.]*//g');
		REC_MDL="$REC_MDL+$GAM";
		break
	else
		ALPHA=-1
	fi
    done
		
    echo "$REC_MDL";
    echo "$PROP";
    echo "$ALPHA";

    if [ $PROP != -1 ] && [ $ALPHA != -1 ]; then
	nice iqtree2 -redo --prefix ${ALN%.*}_cv -s ${ALGN} -te ${ALN}.treefile -st AA -nt $NDS -m $REC_MDL -a $ALPHA -i $PROP -wsl -blfix -wsl -quiet;
    elif [ $PROP != -1 ]; then
	nice iqtree2 -redo --prefix ${ALN%.*}_cv -s ${ALGN} -te ${ALN}.treefile -st AA -nt $NDS -m $REC_MDL -i $PROP -blfix -wsl -quiet;
    elif [ $ALPHA != -1 ]; then
	nice iqtree2 -redo --prefix ${ALN%.*}_cv -s ${ALGN} -te ${ALN}.treefile -st AA -nt $NDS -m $REC_MDL -a $ALPHA -blfix -wsl -quiet;
    else
	nice iqtree2 -redo --prefix ${ALN%.*}_cv -s ${ALGN} -te ${ALN}.treefile -st AA -nt $NDS -m $REC_MDL -blfix -wsl -quiet;
    fi

}

#####################################################################################################
#
#infer phylogenetic tree for each dataset that contains one removed column
#
#####################################################################################################

# seqcon command makes N alignment files, each with one of the N columns omitted (in order, first column to last Nth column)
# files named "test_sc_1.a2m", etc.
./ms_seqcon.py -X --MSA ${ALGN} --form1 fasta

LCO_ALN=$(echo "${ALGN}" | sed 's/[.].*//');

for ALN in ${LCO_ALN}_*.fasta
do
    # run LOOCV for a single deleted column, find ML params
    (
    	if [ $TREE == "None" ]; then
		echo "do_LCO_tree $ALN"
    		do_LCO_tree
    		sleep 3
	else
		echo "do_LCO_tree_starter $ALN"
		do_LCO_tree_starter
		sleep 3
	fi
    )&
    if (( $(jobs -r -p | wc -l) > $N));
    then
        wait -n;
    fi
done

#####################################################################################################
#
#infer the site log-likelihood for the removed column in OG alignment, using the LCO-XV parameters and tree
#
#####################################################################################################

for ALN in ${LCO_ALN}_*.fasta
do
    echo "$ALN -wsl calculation started"
    x=0
    FILE=${ALN}.treefile
    while ! [ -f "$FILE" ] && [ $x -lt 2 ]
    do
        sleep 1m;
	echo "$x min wait for $FILE";
	x=$((x+1));
    done
    echo "moving on to do_wsl for $FILE";
    (
        do_wsl
	sleep 5
    )&
    if (( $(jobs -r -p | wc -l) > $N));
    then 
	    wait -n;
    fi
done

#####################################################################################################
#
#pull the predictive site log-likelihood for the column removed
#
#####################################################################################################

TLH=0
for ALN in ${LCO_ALN}_*.fasta
do
    x=0
    FILE=${ALN%.*}_cv.sitelh;
    echo "$FILE site LL being pulled";
    while ! [ -f "$FILE" ] && [ $x -lt 2 ]
    do
        sleep 1m;
	echo "$x min wait for $FILE";
	x=$((x+1));
    done
    echo "moving on to print out log-likelihood for the $FILE site column"
    egrep 'Site_Lh' ${ALN%.*}_cv.sitelh | sed 's/Site_Lh *//' > ${ALN%.*}_cv2.sitelh
    TMP1=$(echo ${ALN} | sed 's/_sc_.*//g'); 
    COL=$(echo ${ALN#$TMP1} | perl -pe 's/([^0-9].)//g');
    echo "$TMP1 $COL";
    # print out the log-likelihood for the site $COL, save to file
    awk -v col=$COL '{ printf "%f\n", $col }' ${ALN%.*}_cv2.sitelh >> siteloglik.txt
done

#####################################################################################################
#
#clean up
#
#####################################################################################################

mkdir IQTree_Data;
mv ${LCO_ALN}_* IQTree_Data/.;

echo "finished leave-one-column-out cross-validation";

#####################################################################################################
