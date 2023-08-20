while getopts ":i:" opt; do
  case ${opt} in
    i )
      iteration=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      ;;
  esac
done
shift $((OPTIND -1))


disease='Pancan'
#Provide path to wdir where you have loom obj saved
wdir=""

obj_path=$wdir"/data/"$disease"_obj.v2.loom"

#Provide path to the databases
data_path=""

#STEP 1: Run GRN

pyscenic grn $obj_path $data_path/JAPSAR2020_unique_motifs.txt -o $wdir/$iteration/adj.csv --num_workers 80 --sparse



#STEP 2-3: Regulon prediction aka cisTarget from CLI

pyscenic ctx $wdir/$iteration/adj.csv \
    $data_path/cisTarget_databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather $data_path/cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
    --annotations_fname $data_path/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname $obj_path \
    --output $wdir/$iteration/reg.csv \
    --mask_dropouts \
    --num_workers 60


###STEP 4: Cellular enrichment (aka AUCell) from CLI
output_obj_path=$disease"_pyscenic_output.loom"
pyscenic aucell \
    $obj_path \
    $wdir/$iteration/reg.csv \
    --output $wdir/$iteration/$output_obj_path \
    --num_workers 60

