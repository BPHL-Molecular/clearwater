
nextflow.enable.dsl=2

include { elgato4mlst } from './clearwater_modules/elgato.nf'
include { legsta } from './clearwater_modules/legsta.nf'
include { kraken_subs } from './clearwater_modules/kraken_s.nf'
include { proc_mash } from './clearwater_modules/leg_mash.nf'
include { pyProc1 } from './clearwater_modules/pyproc1.nf'
include { prokka } from './clearwater_modules/prokka.nf'
include { pyProc4 } from './clearwater_modules/pyproc4.nf'
include { pyProc5 } from './clearwater_modules/pyproc5.nf'
include { pyProc7 } from './clearwater_modules/pyproc7.nf'
include { ksnp4 } from './clearwater_modules/ksnp4.nf'
include { checkm_qc } from './clearwater_modules/checkm.nf'
include {pyProc3 } from './clearwater_modules/pyproc3.nf'
include { runFastANI } from './clearwater_modules/fastani.nf'
include { pyProc6 } from './clearwater_modules/parseANILog.nf'
include {compile_report } from './clearwater_modules/cum_report.nf'
include { data_preprop } from './clearwater_modules/data_preprop.nf'
include { firth_regression } from './clearwater_modules/firth_regression.nf'
include { cv_confusionMatrix } from './clearwater_modules/cv_cmatrix.nf'
include {iqtree_bootstrap } from './clearwater_modules/iqtree_bootstrap.nf'
include {annotate_tree } from'./clearwater_modules/annotate_tree.nf'
include {PANGENOME_ANALYSIS} from './clearwater_modules/roary_pangenome.nf'
include {PLOT_PANGENOME } from './clearwater_modules/pangenome_plot.nf'
include { EXTRACT_CLOUD_MAP } from './clearwater_modules/cloud_gene_map.nf'
include { SKLEARN_LASSO } from './clearwater_modules/lassoPy.nf'
include {KNN_CLASSIFICATION} from './clearwater_modules/knn_class.nf'
include	{PCA_VISUALIZATION} from './clearwater_modules/pca_visual.nf'
include {calculate_snp_matrix} from './clearwater_modules/csnp_dist_matrix.nf'
include {cluster_transmission_chains} from './clearwater_modules/cluster_transmission.nf'
// cgMLST Surveillance Module Imports (Brand New Module File)
include { run_chewbbaca_cgmlst; generate_cgmlst_report } from './clearwater_modules/cgmlst_surveillance.nf'
include { RUN_MOB_TYPER; MERGE_MOB_REPORTS } from './clearwater_modules/mobSuite_plasmid_surveillance.nf'
include {PLOT_CLOUD_PANGENOME} from './clearwater_modules/cloud_pangenome_plot.nf'
include {PLOT_EXCLUSIVE_CLOUD} from './clearwater_modules/exclusive_cloud_genes.nf


workflow {
    // 1. Core sample channel
    //sorts the files instantly and streams them one by one without creating a blocking barrier.  
    ch_samples = Channel
        .fromPath("${params.input}/*.fasta")
        .toSortedList()
        .flatten()
	
    //ch_samples.view()

    // Dynamic file mapping fallback if parameter is null, blank, or missing
    ref_dir_path = params.ref_genomes_dir ?: 'no_ref_directory'
    
    ch_ref_dir   = Channel.value(file(ref_dir_path))

    // 2. Define the static database file object
    kraken_db_file = file(params.krakenSub_db)

    // 3. Main pipeline execution (up to legsta)
    ch_legsta_out = elgato4mlst(ch_samples) | legsta

    // 4. Kraken and downstream processes
    ch_kraken_out = kraken_subs(kraken_db_file, ch_legsta_out)

    // 5. Connect kSNP4 batch stream directly to a collected samples channel
    ch_collected_fastas = ch_samples.collect()

    // a) run ksnps for core snp tree and core snp alignment
    ksnp4(ch_collected_fastas)

    // b) Bootstrsp tree
    iqtree_bootstrap(ksnp4.out.run_dir)
    
    // c. Calculate the absolute matrix using exact output emitter track
    calculate_snp_matrix( ksnp4.out.run_dir )

    // d. Compute transmission chains using Single-Linkage Clustering
    cluster_transmission_chains( calculate_snp_matrix.out.matrix_tsv )

    // e. cgMLST execution path
    run_chewbbaca_cgmlst( ch_collected_fastas )

    generate_cgmlst_report(
        run_chewbbaca_cgmlst.out.allele_matrix,
        run_chewbbaca_cgmlst.out.qc_report
    )

    // Check your dynamic parameters for metadata file. Cast to a real file path channel if it is a string path,
    // or substitute a placeholder if false or missing.
    def meta_file_ch = (params.metadata && params.metadata != "false" && params.metadata != "None") ? file(params.metadata) : file("EMPTY_METADATA_PLACEHOLDER")

    // f. Plasmid mobilization surveillance track
    // --- SPECIAL ISOLATION FOR MOB_SUITE ---
    // 1. Fork a new channel specifically for mob_typer, transforming paths into [sample_id, path] tuples
    ch_mob_input = ch_samples.map { fasta ->
    def sample_id = fasta.baseName
    return tuple(sample_id, fasta)
    }
    // 2. Feed the isolated tuple channel into parallelized process
    RUN_MOB_TYPER( ch_mob_input )

    // 3. Gather all individual outputs and send them to the compiler
    ch_all_mob_reports = RUN_MOB_TYPER.out.raw_txt.collect()

    ch_Mobmeta_file = Channel.value( meta_file_ch )

    //Combine reports channel with the above metadata channel
    // to ensure Nextflow synchronizes them safely into one step.
    ch_merged_input = ch_all_mob_reports.combine(ch_Mobmeta_file)

    MERGE_MOB_REPORTS( ch_merged_input )

    // 6. Run CheckM QC asynchronously using the same collected FASTA channel
    checkm_qc(ch_collected_fastas)

    // 7. Core pipeline execution chain (up to pyProc5)
    ch_prokka_out = ch_kraken_out
        | proc_mash
        | pyProc1
        | prokka

    // Path A: Continues the original linear workflow exactly as before
    ch_processed_samples = ch_prokka_out
        | pyProc4
        | pyProc5

    // Path B: Isolates and collects .gff files for the downstream pangenome process
    ch_collected_gffs = ch_prokka_out
        | map { kraken_rep, kraken_out, mash, pyout, elgato, fasta_file, prokka_dir ->
            // This safely extracts the specific GFF using Nextflow's file pointer
            return file("${prokka_dir}/${fasta_file.baseName}.gff")
          }
        | collect

    // 7.1. Run multi-sample pangenome computation
    PANGENOME_ANALYSIS( ch_collected_gffs )

    // 7.2. Generate the automated matrix heatmaps and pie charts
    PLOT_PANGENOME( PANGENOME_ANALYSIS.out.pgrun_dir )
    
    // 8. BRIDGE LOGIC: Merge the per-sample flow with the single global CheckM results folder
    ch_with_checkm = ch_processed_samples.combine( checkm_qc.out.results_dir )

    // 9. Single linear pipeline stream execution
    ch_ani_inputs = ch_with_checkm | pyProc3

    // Explicitly feed the sample stream and the reference genomes channel separately
    ch_ani_outputs = runFastANI(ch_ani_inputs, ch_ref_dir)

    // Pipe the outputs directly into the final report compilers
    ch_global_report = ch_ani_outputs | pyProc6 | pyProc7 

    // A. Isolate index 0 (report.txt) and collect into a single list:
    all_reports_ch = ch_global_report.map { tuple -> tuple[0] }.collect()

    // B. Compile the individual sample report rows into one master table
    compile_report(all_reports_ch)

    // Nextflow will wait until BOTH upstream files are fully ready before running this // added on 10132026
    // use dynamic metadata parameters above.

    EXTRACT_CLOUD_MAP(
        PANGENOME_ANALYSIS.out.pgrun_dir,
        compile_report.out.global_report,
        meta_file_ch  // <--- First use
    )

    // C. Run the preprocessing step using ONLY the single table file channel( global_report emitted by compile_report)
    data_preprop(compile_report.out.global_report)

    // D. Feed dataprepop's emitted channel directly into the logistic regression analysis (emit: imputed_report by data_preprop)
    // Firth biased-reduced multinomial logistic regression
    firth_regression(data_preprop.out.imputed_report) // Treats genes'allele counts as continous numerical variables

    // E. Execute internal validation using data and parameters from the global model
    cv_confusionMatrix(
        firth_regression.out.scaled_subspecies_data,
        firth_regression.out.model_object
    )

    // F. Pass metadata into new tree annotation process along with the tree and regression data
    // Pass the active outbreak carriers channel file into the tree annotation module
    annotate_tree(
        iqtree_bootstrap.out.bootstrap_tree,
        firth_regression.out.scaled_subspecies_data,
        meta_file_ch,  // <--- Second use
        cluster_transmission_chains.out.active_carriers_tsv // <--- NEW: Dynamic pipeline link
        )

    // G. Lasso multinomial logistric regression and kNN classification  
    SKLEARN_LASSO(data_preprop.out.imputed_report) // Treats genes'alle numbers as nominal categories

    // k-nearest neiggbors classification
    KNN_CLASSIFICATION(firth_regression.out.scaled_subspecies_data)

    // Principal Component Analyis
    PCA_VISUALIZATION(
        firth_regression.out.scaled_subspecies_data,
        KNN_CLASSIFICATION.out.predictions
    )

    //  =========================================================================
    //  H. CONDITIONAL BRANCH FOR THE ANNOTATED CLOUD HEATMAP MATRIX MODULE
    //  =========================================================================

    // Gather roary outputs, imputed records, and metadata tracks into a single execution array channel
    ch_cloud_plot_inputs = PANGENOME_ANALYSIS.out.pgrun_dir
        .combine( data_preprop.out.imputed_report )
        .combine( ch_Mobmeta_file ) // Use of existing resilient metadata wrapper channel

    // Dispatch parameters to execute the custom visualization script
    PLOT_CLOUD_PANGENOME( ch_cloud_plot_inputs )

    // =========================================================================
    //   I. COMPRESSED PIE DISTRIBUTION & UNIQUE GENE ROSTER TARGET
    // =========================================================================
    // Reuses the safe, synchronized file array channel to isolate data segments
    PLOT_EXCLUSIVE_CLOUD( ch_cloud_plot_inputs )

    
}

    
