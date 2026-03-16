#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# FINALIZE + CLEANUP for your RNA-seq folder architecture
#
# Keeps:
#   - 0_Scripts/
#   - 1_Raw_data/
#   - 14_Counts/
#   - QC reports: 2_FastQC, 3_MultiQC, 7_MultiQC_trimmed, 10_Alignment_QC/multiqc, 12_Merged_BAM_MultiQC, 13_Final_MultiQC
#   - All logs/text files along the way:
#       *.log, *.txt, *.hisat2, *.html, *.json, *.zip, *.tsv, *.csv
#       (copied into a compact 99_Reports_and_logs/ folder)
#
# Deletes (optional with --delete):
#   - 4_Rcorrector/, 5_Filtered/, 6_Trimmed/, 8_HISAT2_index/, 9_Aligned/, 11_Merged_BAM/
#   - plus any other heavy intermediates not explicitly kept
#
# Usage:
#   bash finalize_cleanup.sh            # dry-run (shows what would happen)
#   bash finalize_cleanup.sh --delete   # actually performs deletion
###############################################################################

DO_DELETE=false
REPORT_DIR="99_Reports_and_logs"

for arg in "$@"; do
  case "$arg" in
    --delete) DO_DELETE=true ;;
    *)
      echo "Unknown argument: $arg"
      echo "Usage: $0 [--delete]"
      exit 1
      ;;
  esac
done

# --- Keep these in place (never removed)
KEEP_ALWAYS=("0_Scripts" "1_Raw_data" "14_Counts")

# --- Things you want copied into 99_Reports_and_logs
COPY_DIRS_FULL=("2_FastQC" "3_MultiQC" "7_MultiQC_trimmed" "10_Alignment_QC" "12_Merged_BAM_MultiQC" "13_Final_MultiQC")

echo ">>> Mode: $([[ $DO_DELETE == true ]] && echo DELETE || echo DRY-RUN)"
echo ">>> Reports will be copied to: ${REPORT_DIR}/"
echo

# --- Sanity checks
for d in "${KEEP_ALWAYS[@]}"; do
  if [ ! -d "$d" ]; then
    echo "ERROR: required folder not found: $d"
    exit 1
  fi
done

# --- Create report dir
if [ "$DO_DELETE" = false ]; then
  echo "[DRY-RUN] Would create: ${REPORT_DIR}/"
else
  mkdir -p "$REPORT_DIR"
fi

copy_dir_full() {
  local src="$1"
  local dest="$2"
  if [ ! -d "$src" ]; then
    echo "Skipping missing directory: $src/"
    return 0
  fi
  if [ "$DO_DELETE" = false ]; then
    echo "[DRY-RUN] Would copy directory: $src/  ->  $dest/"
  else
    mkdir -p "$dest"
    rsync -a "${src}/" "${dest}/"
  fi
}

copy_filtered_logs() {
  local src="5_Filtered/logs"
  local dest="${REPORT_DIR}/5_Filtered/logs"
  if [ ! -d "$src" ]; then
    echo "Skipping missing directory: $src/"
    return 0
  fi
  if [ "$DO_DELETE" = false ]; then
    echo "[DRY-RUN] Would copy directory: $src/  ->  $dest/"
  else
    mkdir -p "$dest"
    rsync -a "${src}/" "${dest}/"
  fi
}

copy_trimmed_no_fastq() {
  local src="6_Trimmed"
  local dest="${REPORT_DIR}/6_Trimmed"
  if [ ! -d "$src" ]; then
    echo "Skipping missing directory: $src/"
    return 0
  fi
  if [ "$DO_DELETE" = false ]; then
    echo "[DRY-RUN] Would copy 6_Trimmed excluding *.fq.gz  ->  $dest/"
  else
    mkdir -p "$dest"
    rsync -a --exclude="*.fq.gz" "${src}/" "${dest}/"
  fi
}

copy_hisat2_texts() {
  local src="9_Aligned"
  local dest="${REPORT_DIR}/9_Aligned"
  if [ ! -d "$src" ]; then
    echo "Skipping missing directory: $src/"
    return 0
  fi
  if [ "$DO_DELETE" = false ]; then
    echo "[DRY-RUN] Would copy: ${src}/*.hisat2.txt and ${src}/*.hisat2.log  ->  $dest/"
  else
    mkdir -p "$dest"
    # You asked for .hisat2.txt; your pipeline writes .hisat2.log, so we grab both.
    rsync -a "${src}/"*.hisat2.txt "$dest/" 2>/dev/null || true
    rsync -a "${src}/"*.hisat2.log "$dest/" 2>/dev/null || true
  fi
}

copy_merged_bam_qc() {
  local src="11_Merged_BAM/qc"
  local dest="${REPORT_DIR}/11_Merged_BAM/qc"
  if [ ! -d "$src" ]; then
    echo "Skipping missing directory: $src/"
    return 0
  fi
  if [ "$DO_DELETE" = false ]; then
    echo "[DRY-RUN] Would copy directory: $src/  ->  $dest/"
  else
    mkdir -p "$dest"
    rsync -a "${src}/" "${dest}/"
  fi
}

echo ">>> Copying requested QC/log artifacts..."

# Full-folder copies
for d in "${COPY_DIRS_FULL[@]}"; do
  copy_dir_full "$d" "${REPORT_DIR}/${d}"
done

# Specific copies
copy_filtered_logs
copy_trimmed_no_fastq
copy_hisat2_texts
copy_merged_bam_qc

echo
echo ">>> Copy step complete."

if [ "$DO_DELETE" = false ]; then
  echo ">>> DRY-RUN complete. If correct, run:"
  echo "    bash $0 --delete"
  exit 0
fi

echo
echo ">>> Deleting everything except: ${KEEP_ALWAYS[*]} and ${REPORT_DIR}/"

# Delete everything in the current directory except KEEP_ALWAYS and REPORT_DIR
for item in * .*; do
  # skip current/parent
  [[ "$item" == "." || "$item" == ".." ]] && continue

  # skip keep dirs + reports
  skip=false
  for k in "${KEEP_ALWAYS[@]}"; do
    [[ "$item" == "$k" ]] && skip=true
  done
  [[ "$item" == "$REPORT_DIR" ]] && skip=true

  # If not skip, remove
  if [ "$skip" = false ]; then
    if [ -d "$item" ]; then
      echo "Deleting directory: $item/"
      rm -rf "$item"
    elif [ -f "$item" ]; then
      echo "Deleting file: $item"
      rm -f "$item"
    fi
  fi
done

echo
echo ">>> Done."
echo "Kept in place:"
printf " - %s/\n" "${KEEP_ALWAYS[@]}"
echo "Reports copied to:"
echo " - ${REPORT_DIR}/"
