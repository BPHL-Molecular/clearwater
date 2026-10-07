#!/bin/bash

# (Add cluster directives here like nodes, time, qos, etc.)
#SBATCH --account=bphl-umbrella   # user's group account on cluster otherwise remove it
#SBATCH --qos=bphl-umbrella       # user's group account on cluster otherwise remove it
#SBATCH --mail-user=USER@flhealth.gov  # User's email address
#SBATCH --job-name=clearwater
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          #This parameter shoulbe be equal to the number of samples if you want fastest running speed. However, the setting number should be less than the max cpu limit(150).
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=clearwater.%j.out
#SBATCH --error=clearwater.%j.err



# 1. Load required environment software modules
module load nextflow
module load apptainer

echo "Starting Nextflow pipeline inside a single allocation sandbox..."

# 1. Dynamically find the absolute path to your environment

CONDA_PATH="/PATH-TO/USERNAME/conda/envs/CLEARWATER"

# 2. Execute the workflow using local tracking limits to stay under the QOS submission radar
nextflow run clearwater_wf.nf \
    -params-file params.yaml \
    -profile apptainer \
    --process.withName:'datapreprop|firth_regression|pyProc.*'.container=null \
    --process.withName:'datapreprop|firth_regression|pyProc.*'.conda="$CONDA_PATH" \
    -resume 


NEXTFLOW_EXIT_CODE=$?

# 3. Crash early evaluation check
if [ $NEXTFLOW_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Nextflow pipeline failed with exit code $NEXTFLOW_EXIT_CODE."
    echo "Check .nextflow.log for details."
    exit $NEXTFLOW_EXIT_CODE
fi

echo "Pipeline finished successfully. Checking for report files..."

# 4. Storage sync safety cushion check (Up to 1 minute max polling fallback)
MAX_TRIES=12
TRY_COUNT=0

while [ $TRY_COUNT -lt $MAX_TRIES ]; do
    if find ./output/ -type f -name "report.txt" 2>/dev/null | grep -q .; then
        echo "Found report files! Proceeding to collation..."
        break
    fi
    echo "Waiting for storage sync (Attempt $((TRY_COUNT+1))/$MAX_TRIES)..."
    sleep 5
    TRY_COUNT=$((TRY_COUNT+1))
done

sleep 5

echo "Compiling final summary report..."

# 5. Seed a fresh cumulative file with proper tab-delimited structural headers

echo -e "sampleID\tspeciesID_mash\tnearest_neighbor_mash\tmash_distance\tAssembly_Completeness\tContamination_level\tGenome_size\tNcontig\tsubspeciesID_kraken\tkrakenSubsp_percent\tBest_Ref_Genome\tHighest_ANI_Percent\tSeqType\tGene_flaA\tGene_pilE\tGene_asd\tGene_mip\tGene_mompS\tGene_proA\tGene_neuA_neuAh" > ./output/leg_sum_report.txt

# 6. Stream file text directly while filtering out duplicate header strings natively
find ./output/ -type f -name "report.txt" -exec grep -h -v "^sampleID" {} + >> ./output/leg_sum_report.txt

echo "Summary report leg_sum_report.txt generated successfully!"


# =================================================================
# 6b. NEW: MASTER REPORT AUDIT & VALIDATION CHECK
# =================================================================
echo "=== Running Master Report Validation ==="

# Dynamically parse the input folder directory directly from your params.yaml file
INPUT_DIR=$(grep "^km_size:" -B 5 params.yaml | grep "^input:" | awk '{print $2}' | tr -d '"' | tr -d "'")

# Fallback to standard param location if the layout is simpler
if [ -z "$INPUT_DIR" ]; then
    INPUT_DIR=$(grep "^input:" params.yaml | awk '{print $2}' | tr -d '"' | tr -d "'")
fi

# Count how many raw fasta files were originally fed into the pipeline
EXPECTED_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.fasta" | wc -l)

# Count how many unique data rows are inside the cumulative file (excluding the header line)
ACTUAL_COUNT=$(tail -n +2 ./output/leg_sum_report.txt | awk -F'\t' '{print $1}' | sort -u | wc -l)

echo "Expected samples based on input directory: ${EXPECTED_COUNT}"
echo "Actual unique samples present in summary:  ${ACTUAL_COUNT}"

if [ "${EXPECTED_COUNT}" -eq "${ACTUAL_COUNT}" ]; then
    echo "✅ SUCCESS: All ${ACTUAL_COUNT} samples were successfully processed and aggregated!"
else
    echo "❌ ERROR: Sample mismatch detected! Pipeline dropped $((EXPECTED_COUNT - ACTUAL_COUNT)) sample(s)."
    echo "=== Missing Samples List ==="
    for fasta in "$INPUT_DIR"/*.fasta; do
        sample_name=$(basename "$fasta" .fasta)
        if ! grep -q "^${sample_name}[[:space:]]" ./output/leg_sum_report.txt; then
            echo "  ⚠️ Missing Sample ID: ${sample_name}"
        fi
    done
    exit 1
fi
# =================================================================


# =================================================================
# 7. INTERMEDIATE WORKSPACE PURGE & STORAGE EVACUATION
# =================================================================
echo "=== Validation passed. Cleaning up heavy scratch files ==="

# 1. Force remove the large intermediate work/ directory
rm -rf ./work

# 2. Evacuate container space overhead of apptainer download cache
apptainer cache clean -f 2>/dev/null || true

# 3. Stamp output directory with current runtime tag and archive it
dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt

echo "All pipeline operations complete! Final structured results archived to: ./output-$dt"
