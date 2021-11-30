
from glob import glob
import pandas as pd

metadata = pd.read_csv(
    "data/metadata.tsv",
    header=0,
    sep='\t'
)

patient_ids = metadata.iloc[:, 0]

rule all:
    input:
        expand(
            "final/{patient_id}/{patient_id}_seurat_w_classification.rds",
            patient_id=patient_ids
        )

rule clean:
    shell:
        "rm -rf out final cache/samples/*_w_chromvar.rds"

rule make_assay:
    input:
        h5="data/samples/{patient_id}/filtered_peak_bc_matrix.h5",
        barcodes="data/samples/{patient_id}/singlecell.csv",
        fragments="data/samples/{patient_id}/fragments.tsv.gz"
    output:
        seurat_rds="cache/samples/{patient_id}_seurat_obj.rds"
    shell:
        "scripts/make_assay.R {input} {output}"

rule run_chromvar:
    input:
        rules.make_assay.output.seurat_rds
    output:
        rds_out = "cache/samples/{patient_id}_seurat_w_chromvar.rds",
        chromvar_output="out/chromVAR/{patient_id}_tf_accessibility.tsv.gz"
    shell:
        "scripts/run_chromvar.R {input} {output}"

rule run_classifier:
    input:
        rules.run_chromvar.output.chromvar_output
    output:
        classifier_output="out/classifier/{patient_id}_classification.tsv"
    shell:
        "scripts/run_classifier.py {input} {output}"

rule integrate_classifier_results:
    input:
        rules.run_chromvar.output.chromvar_output,
        rules.run_classifier.output.classifier_output
    output:
        rds_out = "final/{patient_id}/{patient_id}_seurat_w_classification.rds"
    script:
        "scripts/integrate_classifier_results.R"
