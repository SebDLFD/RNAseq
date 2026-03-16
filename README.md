# RNAseq
Scripts for primary and secondary bioinformatic analysis of RNAseq .fq.gz files.

## Primary Bioinformatic.
This part is composed of three scripts that run on bash.

### setup.sh
This script is designed to be used on Linux / Virtual machines such as Ubuntu on windows

1. Install miniconda if not already present
2. Create a RNA-seq environment if it does not already exists
3. Check and install required tools for primary RNAseq bioinformatic used in the RNAseq_Analysis.sh scirpt (rcorrector, fastqc, multiqc, trim-galore, hisat2, samtools, htseq, picard)

### RNAseq_Analysis.sh
This script is designed to be used on a wirtual machine on windows, so must be used by copying pasting chucks of code in a Linux terminal.

This script follows this pipeline :
1. FastQC on Raw data (in .fq.gz format)
2. MultiQC on the FastQC files generated.
3. rCorrector on on RAW reads
4. Filtering unfixable reads flaged by rCorrector using the FilterUncorrectabledPEfastq.py script from https://github.com/harvardinformatics/TranscriptomeAssemblyTools
5. Trim galore on Filtered reads
6. MultiQC on Trimmed reads
7. HISAT2 index build
8. HISAT2 alignment
9. MultiQC on the alignment
10. RG-aware lane merging
11. MultiQC on Merged_BAM files
12. Final MultiQC on all QC files
13. featureCounts for counts table generation

After running this pipeling, you can delete most intermediate files but you should keep : 0_SCRIPTS, 1_RAW DATA, 2_COUNTS, LOG FILES, FASTQC AND MULTIQC REPORTS (FastQC on raw data, MultiQC on raw data, Filtered logs, Trimmed FastQC and reports, MultiQC_trimmed, Aligned Text files, Alignment_QC, Merged_BAM_QC, Merged_BAM_MultiQC, Final_MultiQC).
Alternatively, you can use the finalize_cleanup.sh script to do all that for you.

### finalize_cleanup.sh
Cleanup the RNAseq folder architecture

Keeps:
1. 0_Scripts
2. 1_Raw_data
3. 14_Counts
4. QC reports: 2_FastQC, 3_MultiQC, 7_MultiQC_trimmed, 10_Alignment_QC/multiqc, 12_Merged_BAM_MultiQC, 13_Final_MultiQC
5. All logs/text files along the way: *.log, *.txt, *.hisat2, *.html, *.json, *.zip, *.tsv, *.csv (copied into a compact 99_Reports_and_logs/ folder)

Deletes (optional with --delete): 4_Rcorrector/, 5_Filtered/, 6_Trimmed/, 8_HISAT2_index/, 9_Aligned/, 11_Merged_BAM/ plus any other heavy intermediates not explicitly kept

Usage:
bash finalize_cleanup.sh            # shows what would happen without deleting
bash finalize_cleanup.sh --delete   # actually performs deletion


## Secondary Bioinformatics:
One script is used for this part, and works on R. It uses the featureCounts_gene_counts.txt file produced by featureCounts.
This script follows these different steps :
1. Base setup : reading the count table, cleaning sample column names, creating a mapping table ussing .gff
2. Secondary setup to compare APO samples from other samples by building a DESeq object a making a VST matrix for heatmaps
3. PCA with all the samples
4. Heatmap of all transcript differentially expressed
5. Visualization of top expressed genes in APO (expression-based)
6. Visualization of top APO markers (highest DE genes)
7. Specific gene list
8. GO enrichment with fgsea/GSEA
9. GO enrichment with GO_MWU (RBGOA)
10. WGCNA FOR TFs
