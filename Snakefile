from glob import glob
import pandas as pd

metadata = pd.read_csv(
    "data/metadata.tsv",
    header=0,
    sep='\t'
)

patient_ids = metadata.iloc[:, 0]

rule all:
    """Dummy rule to generate all outputs for samples in metadata."""
    input:
        expand(
            "final/{patient_id}/{patient_id}_seurat_w_classification.rds",
            patient_id=patient_ids
        )

rule clean:
    """Convenience rule for cleaning up run."""
    shell:
        "rm -rf out final"

rule make_assay:
    """Import 10X output into Seurat data object and export to disk."""
    input:
        h5="data/samples/{patient_id}/filtered_peak_bc_matrix.h5",
        barcodes="data/samples/{patient_id}/singlecell.csv",
        fragments="data/samples/{patient_id}/fragments.tsv.gz"
    output:
        seurat_rds="cache/samples/{patient_id}_seurat_obj.rds"
    shell:
        "scripts/make_assay.R {input} {output}"

rule run_chromvar:
    """Calculate gene activity with chromVAR and write to file."""
    input:
        rules.make_assay.output.seurat_rds
    output:
        rds_out = temp("cache/samples/{patient_id}_seurat_w_chromvar.rds"),
        chromvar_output="out/chromVAR/{patient_id}_tf_accessibility.tsv.gz"
    shell:
        "scripts/run_chromvar.R {input} {output}"

rule run_classifier:
    """Execute classification by decision tree on chromVAR output and generate proba matrix."""
    input:
        rules.run_chromvar.output.chromvar_output
    output:
        classifier_output="out/classifier/{patient_id}_classification.tsv"
    shell:
        "scripts/run_classifier.py {input} {output}"

rule integrate_classifier_results:
    """Add classifier output results to Seurat object as new assay."""
    input:
        rules.run_chromvar.output.rds_out,
        rules.run_classifier.output.classifier_output
    output:
        rds_out = "final/{patient_id}/{patient_id}_seurat_w_classification.rds"
    shell:
        "scripts/integrate_classifier_results.R {input} {output}"
