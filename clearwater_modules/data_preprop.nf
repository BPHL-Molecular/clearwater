process data_preprop {
    container null
    tag { "Preprocessing Global Merged Imputed Data Report" }
    publishDir "${params.output}", mode: 'copy'

    conda 'conda-forge::pandas conda-forge::miceforest'

    input:
    path cumulative_report

    output:
    path "imputed_subspecies_report.txt", emit: imputed_report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np

    # Load the clean compiled multi-row dataset directly
    df = pd.read_csv("${cumulative_report}", sep="\\t")

    target_cols = ['Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh']
    existing_targets = [col for col in target_cols if col in df.columns]

    print(f"DEBUG: Total collective rows processed from master file: {len(df)}")

    if existing_targets:
        # Step 1: Clean spaces, then explicitly map complete gene copy absences to "0" strings first
        for col in existing_targets:
            df[col] = df[col].astype(str).str.strip()
            # Catches both lone dashes and MD- tags indicating structural deletion
            df[col] = df[col].replace({'-': '0', 'MD-': '0'})

        # Step 2: Define remaining missingness items (excluding the dash since it is now "0")
        missing_variants = ['N/A', 'None', 'null', 'NAT', 'MA?', '?', '']
        df[existing_targets] = df[existing_targets].replace(missing_variants, np.nan)

        # Step 3: Coerce safely to numeric. The "0" strings cleanly become structural numeric 0.0 values.
        for col in existing_targets:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        # Step 4: Run miceforest imputation ONLY over the remaining unresolved NaN fields
        if df[existing_targets].isnull().values.any():
            if len(df) > 15:
                try:
                    import miceforest as mf
                    # Imputation kernel runs smoothly because true 0.0 entries are preserved as complete historical entries
                    kernel = mf.ImputationKernel(df[existing_targets].astype(float), datasets=1, random_state=42)
                    kernel.mice(iterations=3)
                    df_imputed = kernel.complete_data(dataset=0)
                    for col in existing_targets:
                        df[col] = df_imputed[col]
                except Exception as e:
                    print(f"WARNING: miceforest failed due to: {e}. Falling back to median imputation.")
                    for col in existing_targets:
                        col_median = df[col].median()
                        if np.isnan(col_median):
                            col_median = 0.0
                        df[col] = df[col].fillna(col_median)
            else:
                for col in existing_targets:
                    col_median = df[col].median()
                    if np.isnan(col_median):
                        col_median = 0.0
                    df[col] = df[col].fillna(col_median)

    df.to_csv("imputed_subspecies_report.txt", sep="\t", index=False)
    """
}
