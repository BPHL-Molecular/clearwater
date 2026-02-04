#!/usr/bin/env bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --mail-user=Tassy.Bazile@flhealth.gov
#SBATCH --job-name=clearwater
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          #This parameter shoulbe be equal to the number of samples if you want fastest running speed. However, the setting number should be less than the max cpu limit(150). 
#SBATCH --mem=96gb
#SBATCH --time=48:00:00
#SBATCH --output=clearwater.%j.out
#SBATCH --error=clearwater.%j.err

module load nextflow apptainer
#nextflow run clearwater_wf.nf -params-file params.yaml -with-trace > leg.txt # redirects the standard output of the Nextflow execution to a file named leg.txt.

############################################################################
nextflow run clearwater_wf.nf -params-file params.yaml -with-trace > leg.txt
############################################################################

# concatenaing all individual (per sampleID) reports into one final summary report
sort ./output/*/report.txt | uniq > ./output/leg_sum_report.txt
sed -i '/sampleID\tspeciesID/d' ./output/leg_sum_report.txt # perform a delete in the specific pattern('sampleID\tspeciesID')
sed -i '1i sampleID\tspeciesID_mash\tnearest_neighb_mash\tmash_distance\tsubspeciesID_kraken\tkrakenSubsp_percent\tSeqType\tGene_flaA\tGene_pilE\tGene_asd\tGene_mip\tGene_mompS\tGene_neuA_neuAH' ./output/leg_sum_report.txt

#workid="$(grep -o -m 1 '\[.*/.*\]' leg.txt | cut -d ' ' -f 1 | cut -d '[' -f 2 | cut -d ']' -f 1)" 
#w2="$(realpath ./work/$workid*)"
#cp -r $w2/* ./output
#rm leg.txt

#Adding date to output dir, remove work and cache to clear space.
#singularity cache clean -f 
dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt
#mv ./work ./work-$dt
rm -r ./work
#rm -r ./cache
