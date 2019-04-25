
DATADIR="data/demo2"
OUTPUTDIR="data/demo2_output"

ARGS="$DATADIR/args.txt"

JOB_NUM=3

function getline { echo $(sed "$1 q;d" $ARGS); }

OTU_POS=$((($JOB_NUM-1) * 5 + 1))
SAMPLE_POS=$(($OTU_POS + 1))
REP_POS=$(($SAMPLE_POS + 1))
SEED_POS=$(($REP_POS + 1))

OTU=$(getline $OTU_POS)
SAMPLE=$(getline $SAMPLE_POS)
REP=$(getline $REP_POS)
SEED=$(getline $SEED_POS)

FILE_STRUCT="otu_size-$OTU--sample_size-$SAMPLE--rep-$REP--seed-$SEED"
BASE="$DATADIR/$FILE_STRUCT"
TABLE="$BASE.biom"
TREE="$BASE.newick"

OUTPUT_FILE="$OUTPUTDIR/$FILE_STRUCT--output.txt"

(/usr/bin/time -l python benchmark/time_skbio_faith.py $TABLE $TREE) &> $OUTPUT_FILE

echo "seed-$SEED\n$(cat $OUTPUT_FILE)" > $OUTPUT_FILE
echo "rep-$REP\n$(cat $OUTPUT_FILE)" > $OUTPUT_FILE
echo "sample-$SAMPLE\n$(cat $OUTPUT_FILE)" > $OUTPUT_FILE
echo "otu-$OTU\n$(cat $OUTPUT_FILE)" > $OUTPUT_FILE
