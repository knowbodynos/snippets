#!/bin/bash

# To get kmer count tsv: kmercountexact.sh fastadump=f k=5 in=five_prime_utr_u1000d0.fa out=five_prime_utr_u1000d0_k5.tsv

# Print usage
usage() {
  echo -n "$(basename $0) [-udeoh]... [GTF File] [Fasta File] [Region type]

 Generate GTF and Fasta files for a subset region of genome.

 Options:
  -u    Number of upstream base pairs to include.
  -d    Number of downstream base pairs to include.
  -e    Specify exon number (if it exists).
  -o    Output directory.
  -h    Display this help and exit.
" 1>&2; exit 1;
}

# Set initial argument values

UPSTREAM=0
DOWNSTREAM=0
EXON=1
OUT_DIR=$(pwd)

while getopts ":u:d:e:o:h" o; do
    case "${o}" in
        u)
            UPSTREAM=${OPTARG}
            [ "${UPSTREAM}" -ge 0  ] || { echo -e "\nUpstream must be a non-negative integer.\n" && usage; }
            ;;
        d)
            DOWNSTREAM=${OPTARG}
            [ "${DOWNSTREAM}" -ge 0 ] || { echo -e "\nDownstream must be a non-negative integer.\n" && usage; }
            ;;
        e)
            EXON=${OPTARG}
            [ "${EXON}" -gt 0 ] || { echo -e "\Exon number must be a positive integer.\n" && usage; }
            ;;
        o)
            OUT_DIR="$(cd $(dirname ${OPTARG}); pwd)/$(basename ${OPTARG})"
            if [ ! -d "${OUT_DIR}" ]
            then
                mkdir -p ${OUT_DIR}
            fi
            ;;
        h)
            usage
            ;;
        *)
            break
            ;;
    esac
done
shift $((OPTIND-1))

# Positional arguments
GTF_PATH=$1
FASTA_PATH=$2
REGION_NAME=$3

# Check if arguments are set
if [ -z "${GTF_PATH}" ] || [ -z "${FASTA_PATH}" ] || [ -z "${REGION_NAME}" ] || [ "${UPSTREAM}" -lt 0 ] || [ "${DOWNSTREAM}" -lt 0 ] || [ "${EXON}" -le 0 ] || [ -z "${OUT_DIR}" ]
then
    usage
fi

# Subset GTF file
cat ${GTF_PATH} | awk '{line=$0}; $3 == region {print line}' region="${REGION_NAME}" > ${OUT_DIR}/${REGION_NAME}.gtf


# Subset exon number if specified
if [ "$(head -n 1 ${OUT_DIR}/${REGION_NAME}.gtf | grep 'exon_number' | wc -l)" -gt 0 ]
then
    REGION_NAME_AUG1=${REGION_NAME}_e1
    cat ${OUT_DIR}/${REGION_NAME}.gtf | awk '{line=$0}; $14 == exon_number {print line}' exon_number="\"${EXON}\";" > ${OUT_DIR}/${REGION_NAME_AUG1}.tmp.gtf
    rm ${OUT_DIR}/${REGION_NAME}.gtf
    mv ${OUT_DIR}/${REGION_NAME_AUG1}.tmp.gtf ${OUT_DIR}/${REGION_NAME_AUG1}.gtf
else
    REGION_NAME_AUG1=${REGION_NAME}
fi

if [ "${UPSTREAM}" -gt 0 ] || [ "${DOWNSTREAM}" -gt 0 ]
then
    REGION_NAME_AUG2=${REGION_NAME_AUG1}_u${UPSTREAM}d${DOWNSTREAM}
    CHROM_SIZES_PATH="$(dirname ${FASTA_PATH})/chrom.sizes"
    if [ ! -f "${CHROM_SIZES_PATH}" ]
    then
        if [ ! -f "${FASTA_PATH}.fai" ]
        then
            # Generate Fasta index
            samtools faidx ${FASTA_PATH}
        fi
        # Generate chrom.sizes file (i.e. bedtools genome file)
        cut -f 1,2 ${FASTA_PATH}.fai > ${CHROM_SIZES_PATH}
    fi
    # Include upstream/downstream basees in subset GTF file
    bedtools slop -i ${OUT_DIR}/${REGION_NAME_AUG1}.gtf -g ${CHROM_SIZES_PATH} -l ${UPSTREAM} -r ${DOWNSTREAM} -s > ${OUT_DIR}/${REGION_NAME_AUG2}.gtf
    rm ${OUT_DIR}/${REGION_NAME_AUG1}.gtf
else
    REGION_NAME_AUG2=${REGION_NAME_AUG1}
fi

# Generate subset Fasta file
bedtools getfasta -fi ${FASTA_PATH} -bed ${OUT_DIR}/${REGION_NAME_AUG2}.gtf -fo ${OUT_DIR}/${REGION_NAME_AUG2}.fa

# Generate subset tsv file
cat ${OUT_DIR}/${REGION_NAME_AUG2}.fa | paste -sd '\t\n' | tr -d '>' | awk '{gsub(/[:-]/," ",$1)}; 1' > ${OUT_DIR}/${REGION_NAME_AUG2}.tmp.tsv
echo -e "Geneid\tChr\tStart\tEnd\tStrand\tSeq" > ${OUT_DIR}/${REGION_NAME_AUG2}.tsv
cat ${OUT_DIR}/${REGION_NAME_AUG2}.gtf | awk '{$1=$1}; {gsub(/[\";]/,"",$10)}; {print $10,$7}' | paste - ${OUT_DIR}/${REGION_NAME_AUG2}.tmp.tsv | awk '{strand=$2}; {for(i=2;i<5;i++) $i=$(i+1)}; {$5=strand}; {print}'  | tr ' ' '\t' >> ${OUT_DIR}/${REGION_NAME_AUG2}.tsv
rm ${OUT_DIR}/${REGION_NAME_AUG2}.tmp.tsv
