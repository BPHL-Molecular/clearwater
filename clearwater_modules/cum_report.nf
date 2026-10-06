process compile_report {
    container null
    tag { "Compiling individual sample reports" }
   
    input:
    // FIXED: stageAs: 'report_*.txt' dynamically renames duplicates 
    // (e.g., report_0.txt, report_1.txt) so they can live in the same directory safely
    path individual_reports, stageAs: 'report_*.txt'

    output:
    path "lp_sum_report.txt", emit: global_report

    script:
    """
    # 1. Seed a fresh cumulative file with proper tab-delimited structural headers
    echo -e "sampleID\tspeciesID_mash\tnearest_neighbor_mash\tmash_distance\tAssembly_Completeness\tContamination_level\tGenome_size\tNcontig\tsubspeciesID_kraken\tkrakenSubsp_percent\tBest_Ref_Genome\tHighest_ANI_Percent\tSeqType\tGene_flaA\tGene_pilE\tGene_asd\tGene_mip\tGene_mompS\tGene_proA\tGene_neuA_neuAh" > lp_sum_report.txt

    # 2. Append all dynamically named sample files into it while stripping duplicate headers safely
    cat report_*.txt | grep -v "^sampleID" >> lp_sum_report.txt
    """
}