#!/bin/bash

###############################################################################
# This script :
# 1. Install miniconda if not already present
# 2. Create a RNA-seq environment if it does not already exists
# 3. Check and install required tools for RNAseq 
###############################################################################


set -euo pipefail

ENV_NAME="rnaseq"
TOOLS=("rcorrector" "fastqc" "multiqc" "trim-galore" "hisat2" "samtools" "htseq" "picard")

if ! command -v conda >/dev/null 2>&1; then
    echo ">>> conda not found. Installing Miniconda..."

    cd /tmp
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p "$HOME/miniconda3"
    rm Miniconda3-latest-Linux-x86_64.sh

    "$HOME/miniconda3/bin/conda" init bash
    eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"

    echo ">>> Miniconda installed in $HOME/miniconda3"
else
    echo ">>> conda already installed."
    eval "$(conda shell.bash hook)"
fi

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo ">>> Conda environment '$ENV_NAME' already exists."
else
    echo ">>> Creating environment '$ENV_NAME'..."
    conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda
    echo ">>> Environment '$ENV_NAME' created."
fi

echo ">>> Activating environment '$ENV_NAME'"
conda activate "$ENV_NAME"

for TOOL in "${TOOLS[@]}"; do
    if conda list | awk '{print $1}' | grep -qx "$TOOL"; then
        echo ">>> $TOOL already installed."
    else
        echo ">>> Installing $TOOL ..."
        conda install -y -c conda-forge -c bioconda "$TOOL"
        echo ">>> Installed $TOOL."
    fi
done

echo
echo ">>> Setup complete. RNA-seq environment ready."
echo ">>> Activate it anytime with:"
echo "    conda activate $ENV_NAME"
echo
echo "Installed tools:"
printf " - %s\n" "${TOOLS[@]}"
