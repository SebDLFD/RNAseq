# =========================================
# 0.0 First time setup (only once)
# =========================================

# BEFORE RUNNING ANYTHING FROM THIS SCRIPT, MAKE SURE THE FOLLOWING IS DONE

# Make sure the raw sequencing data are in the folder called 1_Raw_data in your RNAseq analysis main folder
# Make sure the 0_Scripts folder is in the RNAseq analysis folder with the script called setup.sh inside 

# In the Ubuntu terminal, install dos2unix:
sudo apt update
sudo apt install dos2unix

# change directory to your analysis folder:
cd /mnt/c/Users/Sebastien/Documents/Postdoc_files/Experiments_and_projects/3_Results/3_Whole_body_organ_by_organ_RNAseq/2_RNAseq_data_and_analysis/4_Full_resequencing_2nd_RNAseq_round/Primary_Bioinformatics

# Make sure that the setup.sh script has Unix Line endings and is executable :
dos2unix 0_Scripts/setup.sh
chmod +x 0_Scripts/setup.sh

# Run the setup.sh script to install Miniconda and your environement
./0_Scripts/setup.sh

# Download FilterUncorrectabledPEfastq.py script
cd 0_Scripts
wget https://raw.githubusercontent.com/harvardinformatics/TranscriptomeAssemblyTools/master/utilities/FilterUncorrectabledPEfastq.py
chmod +x FilterUncorrectabledPEfastq.py
cd ..

# Create a directory structure
mkdir -p 1_Raw_data
mkdir -p 2_FastQC
mkdir -p 3_MultiQC
mkdir -p 4_Rcorrector
mkdir -p 5_Filtered 5_Filtered/logs
mkdir -p 6_Trimmed
mkdir -p 7_MultiQC_trimmed
mkdir -p 8_HISAT2_index
mkdir -p 9_Aligned
mkdir -p 10_Alignment_QC 10_Alignment_QC/multiqc
mkdir -p 11_Merged_BAM 11_Merged_BAM/qc 11_Merged_BAM/rg_fixed 11_Merged_BAM/tmp
mkdir -p 12_Merged_BAM_MultiQC
mkdir -p 13_Final_MultiQC
mkdir -p 14_Counts

# =========================================
# 0.1 Pre-session step (every time you open a new terminal)
# =========================================
# Activate your working environment:
conda activate rnaseq

# Change directory to your analysis folder:
cd /mnt/c/Users/Sebastien/Documents/Postdoc_files/Experiments_and_projects/3_Results/2_Whole_body_organ_by_organ_RNAseq/2_RNAseq_data_and_analysis/4_Full_resequencing_2nd_RNAseq_round/Primary_Bioinformatics

# You can now copy/paste blocks into the terminal.


# =========================================
# 1. FastQC on RAW data (compressed, ending with .fq.gz)
# =========================================

fastqc -t 15 1_Raw_data/*.fq.gz -o 2_FastQC/


# =========================================
# 2. MultiQC on RAW FastQC reports
# =========================================

multiqc 2_FastQC/ -o 3_MultiQC/


# =========================================
# 3. Rcorrector on RAW reads
# =========================================

# Run Rcorrector for all the .gz files of the Rawdata folder, with a stop if an error occurs
set -euo pipefail

for r1 in 1_Raw_data/*_1.fq.gz; do
    [ -e "$r1" ] || { echo "No 1_Raw_data/*_1.fq.gz files found"; break; }

    r2="${r1/_1.fq.gz/_2.fq.gz}"
    if [ ! -e "$r2" ]; then
        echo "WARNING: Cannot find mate for $r1 (expected $r2), skipping"
        continue
    fi

    sample=$(basename "$r1" _1.fq.gz)
    echo "=== Correcting: $sample ==="

    run_rcorrector.pl \
        -t 16 \
        -od 4_Rcorrector \
        -1 "$r1" \
        -2 "$r2"
done


# =========================================
# 5. Filter "unfixable" reads (Rcorrector flags)
# =========================================

# Run FilterUncorrectabledPEfastq.py for all the files in the 4_Rcorrector folder
for r1 in 4_Rcorrector/*_1.cor.fq.gz; do
    [ -e "$r1" ] || { echo "No 4_Rcorrector/*_1.cor.fq.gz files found"; break; }
    
    r2="${r1/_1.cor.fq.gz/_2.cor.fq.gz}"
    sample=$(basename "$r1" _1.cor.fq.gz)

    if [ ! -e "$r2" ]; then
        echo "WARNING: Cannot find mate for $r1 (expected $r2), skipping"
        continue
    fi

    echo "=== Filtering unfixable reads for: $sample ==="

    python 0_Scripts/FilterUncorrectabledPEfastq.py -1 "$r1" -2 "$r2"
    if [ -e 4_Rcorrector/rmunfixable_None.log ]; then
  mv 4_Rcorrector/rmunfixable_None.log "5_Filtered/logs/rmunfixable_${sample}.log"
fi
done

# Renaming log files
cd 5_Filtered/logs/

for f in rmunfixable_Bg_*.log; do
  base=$(basename "$f")

  # 1) Extract Bg_<ORGAN> (e.g., Bg_INT, Bg_HEART)
  organ=$(echo "$base" | sed -E 's/^rmunfixable_(Bg_[^_]+).*/\1/')

  # 2) Extract lane at the end (e.g., L1, L4, L12)
  lane=$(echo "$base" | sed -E 's/.*_(L[0-9]+)\.log$/\1/')

  # If lane extraction fails, keep it as "LNA"
  if [[ "$lane" == "$base" ]]; then
    lane="LNA"
  fi

  new="rmunfixable_${organ}_${lane}.log"

  # Safety: avoid overwriting
  if [ -e "$new" ] && [ "$new" != "$f" ]; then
    echo "WARNING: $new already exists, skipping $f"
    continue
  fi

  echo "Renaming $f -> $new"
  mv "$f" "$new"
done

cd ..
cd ..

# Mooving the filtered reads under 5_Filtered
mv 4_Rcorrector/unfixrm_* 5_Filtered/

# Removing the unfixrm_ prefix for simpler data manipulation, original filtering logs are kept
for f in unfixrm_*.cor.fq.gz; do
    newname="${f#unfixrm_}" 
    mv "$f" "$newname"
done

cd ..


# =========================================
# 6. Trim galore on Filtered reads to remove adaptaters and low quality ends
# =========================================

# Run trim galore for all the files in the 5_Filtered folder and store results in 6_Trimmed
# 4 jobs are run in parallel, each using 8 threads, to adjust if needed
MAX_JOBS=4
CURRENT_JOBS=0

for r1 in 5_Filtered/*_1.cor.fq.gz; do
    r2="${r1/_1.cor.fq.gz/_2.cor.fq.gz}"
    sample=$(basename "$r1" _1.cor.fq.gz)

    echo "=== Launching Trim Galore for: $sample ==="

    trim_galore --paired \
        --fastqc \
        --fastqc_args "--threads 8" \
        -o 6_Trimmed \
        "$r1" "$r2" &

    ((CURRENT_JOBS++))

    if [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ]; then
        wait
        CURRENT_JOBS=0
    fi
done

wait  
echo "All Trim Galore jobs completed."


# =========================================
# 7. MultiQC on Trimmer reads
# =========================================

multiqc 6_Trimmed/ -o 7_MultiQC_trimmed/


# =========================================
# 8. HISAT2 index
# =========================================

# Define the main working folder
MAIN="/mnt/c/Users/Sebastien/Documents/Postdoc_files/Experiments_and_projects/3_Results/2_Whole_body_organ_by_organ_RNAseq/2_RNAseq_data_and_analysis/4_Full_resequencing_2nd_RNAseq_round/Primary_Bioinformatics"

# Download the genome.fna and genome.gtf files in NCBI or elsewere and place them in 8_HISAT2_index folder
# Create HISAT2 index
cd 8_HISAT2_index
hisat2-build -p 16 xgBioGlab_20251119.hap1.fa xgBioGlab_20251119
cd ..

# Defining Index prefix
HISAT2_INDEX="$MAIN/8_HISAT2_index/xgBioGlab_20251119"

# =========================================
# 9. HISAT2 alignment
# =========================================

# Align trimmed reads with HISAT2, output directly into samtools, no SAM written 
# (See at the end of the script if you want to keep the SAM files)
set -euo pipefail

for r1 in 6_Trimmed/*_1.cor_val_1.fq.gz; do
    [ -e "$r1" ] || break

    r2="${r1/_1.cor_val_1.fq.gz/_2.cor_val_2.fq.gz}"
    sample=$(basename "$r1" _1.cor_val_1.fq.gz)

    [[ -e "$r2" ]] || continue

    echo "=== Aligning + sorting: $sample ==="

    hisat2 -p 16 \
        -x "$HISAT2_INDEX" \
        -1 "$r1" \
        -2 "$r2" \
        --dta \
        2> "9_Aligned/${sample}.hisat2.log" | \
    samtools sort -@ 8 \
        -o "9_Aligned/${sample}.sorted.bam"

    samtools index "9_Aligned/${sample}.sorted.bam"

    samtools flagstat "9_Aligned/${sample}.sorted.bam" \
        > "10_Alignment_QC/${sample}.flagstat.txt"

    samtools idxstats "9_Aligned/${sample}.sorted.bam" \
        > "10_Alignment_QC/${sample}.idxstats.txt"
done

# =========================================
# 10. MultiQC on the alignment
# =========================================

multiqc 9_Aligned 10_Alignment_QC -o 10_Alignment_QC/multiqc


# =========================================
# 11. RG-aware lane merge
# =========================================
set -euo pipefail

THREADS=8

mapfile -t samples < <(
  ls 9_Aligned/*.sorted.bam \
    | xargs -n1 basename \
    | sed -E 's/_[^_]+_L[0-9]+\.sorted\.bam$//' \
    | sort -u
)

for s in "${samples[@]}"; do
  mapfile -t bams < <(ls "9_Aligned/${s}"_*_L*.sorted.bam 2>/dev/null || true)
  [ "${#bams[@]}" -gt 0 ] || continue

  echo "=== Sample: $s  (lanes: ${#bams[@]}) ==="

  rg_bams=()

  for bam in "${bams[@]}"; do
    fn="$(basename "$bam")"
    core="${fn%.sorted.bam}"    # remove .sorted.bam
    core="${core%.bam}"         # safety if needed

    # Parse RUNID and lane from: <sample>_<RUNID>_L<lane>.sorted.bam
    # Example core: Bg_AG_...-1A_23HJW5TL4_L1
    runid="$(echo "$core" | sed -E 's/^.*_([^_]+)_L[0-9]+$/\1/')"
    lane="$(echo "$core" | sed -E 's/^.*_(L[0-9]+)$/\1/')"

    # RG fields
    RGSM="$s"                 # sample
    RGLB="$s"                 # library (same prep/tube across lanes -> same LB)
    RGPL="ILLUMINA"
    RGPU="${runid}.${lane}"   # platform unit
    RGID="${runid}.${lane}"   # RG id

    out="11_Merged_BAM/rg_fixed/${core}.rg.bam"

    echo "  - RG: SM=$RGSM LB=$RGLB PU=$RGPU ID=$RGID from $fn"

    picard AddOrReplaceReadGroups \
      I="$bam" \
      O="$out" \
      RGSM="$RGSM" \
      RGLB="$RGLB" \
      RGPL="$RGPL" \
      RGPU="$RGPU" \
      RGID="$RGID" \
      VALIDATION_STRINGENCY=SILENT \
      TMP_DIR="11_Merged_BAM/tmp"

    samtools index "$out"
    rg_bams+=("$out")
  done

  if [ "${#rg_bams[@]}" -eq 1 ]; then
    ln -sf "../${rg_bams[0]}" "11_Merged_BAM/${s}.merged.sorted.bam"
  else
    samtools merge -@ "$THREADS" -f "11_Merged_BAM/${s}.merged.bam" "${rg_bams[@]}"
    samtools sort  -@ "$THREADS" -o "11_Merged_BAM/${s}.merged.sorted.bam" "11_Merged_BAM/${s}.merged.bam"
    rm -f "11_Merged_BAM/${s}.merged.bam"
  fi

  samtools index "11_Merged_BAM/${s}.merged.sorted.bam"
  samtools flagstat "11_Merged_BAM/${s}.merged.sorted.bam" > "11_Merged_BAM/qc/${s}.merged.flagstat.txt"
  samtools idxstats "11_Merged_BAM/${s}.merged.sorted.bam" > "11_Merged_BAM/qc/${s}.merged.idxstats.txt"
done

# Sanity checks
for s in $(ls 9_Aligned/*.sorted.bam | xargs -n1 basename | sed -E 's/_[^_]+_L[0-9]+\.sorted\.bam$//' | sort -u); do
  echo "$s"
  ls 9_Aligned/${s}_*_L*.sorted.bam | sed 's#^#  - #'
done

# FOR SINGLE LANE SAMPLES, DELETE THE 0 KB .BAM FILES AND COPY PASTE IN THE FOLDER 11_Merged_BAM, FILES FROM THE rg_fixed SUBFOLDER

# Check that the renaming will work
for bam in 11_Merged_BAM/*.bam; do
  base=$(basename "$bam" .bam)
  short=$(echo "$base" | cut -d'_' -f2)
  echo "$bam  -->  11_Merged_BAM/${short}.bam"
done

# Rename all files so that they are with the following format "sample name".bam or .bam.bai
set -euo pipefail

for bam in 11_Merged_BAM/*.bam; do
  base=$(basename "$bam" .bam)
  short=$(echo "$base" | cut -d'_' -f2)

  new="11_Merged_BAM/${short}.bam"

  # Safety check: avoid accidental overwrites
  if [ -e "$new" ] && [ "$bam" != "$new" ]; then
    echo "ERROR: $new already exists. Skipping $bam"
    continue
  fi

  echo "Renaming: $bam -> $new"
  mv "$bam" "$new"

  # Rename index if it exists
  if [ -e "${bam}.bai" ]; then
    mv "${bam}.bai" "${new}.bai"
  elif [ -e "${bam%.bam}.bai" ]; then
    mv "${bam%.bam}.bai" "${new}.bai"
  fi
done

# =========================================
# 12. MultiQC on Merged_BAM files
# =========================================

multiqc 11_Merged_BAM -o 12_Merged_BAM_MultiQC


# =========================================
# 13. Final MultiQC on all QC files
# =========================================

multiqc 2_FastQC 6_Trimmed 9_Aligned 10_Alignment_QC 11_Merged_BAM -o 13_Final_MultiQC/


# =========================================
# 14. featureCounts for gene-level counts
# =========================================
MAIN="/mnt/c/Users/Sebastien/Documents/Postdoc_files/Experiments_and_projects/3_Results/2_Whole_body_organ_by_organ_RNAseq/2_RNAseq_data_and_analysis/4_Full_resequencing_2nd_RNAseq_round/Primary_Bioinformatics"
GENOME_GFF="$MAIN/8_HISAT2_index/xgBioGlab_hap1_braker_20260209.gff3"

featureCounts -T 16 \
  -p -B -C \
  -F GFF \
  -t exon \
  -g gene_id \
  -a "$GENOME_GFF" \
  -o "14_Counts/featureCounts_gene_counts.txt" \
  11_Merged_BAM/*.bam \
  2> "14_Counts/featureCounts.log"
  
# After running this pipeline, you can delete most intermediate files but you should keep and can rename :
# 0_SCRIPTS
# 1_RAW DATA
# 2_COUNTS
# LOG FILES, FASTQC AND MULTIQC REPORTS:
#     2_FastQC on raw data
#     3_MultiQC on raw data
#     5_Filtered logs
#     6_Trimmed FastQC and reports
#     7_MultiQC_trimmed
#     9_Aligned Text files
#     10_Alignment_QC
#     11_Merged_BAM_QC
#     12_Merged_BAM_MultiQC
#     13_Final_MultiQC

# Alternatively you can run the finalize_cleanup.sh script that will do that for you
dos2unix 0_Scripts/finalize_cleanup.sh
chmod +x 0_Scripts/finalize_cleanup.sh
bash 0_Scripts/finalize_cleanup.sh #To check that everything is right
bash 0_Scripts/finalize_cleanup.sh --delete #To delete



