version 1.0

task call_outliers{
    input{
        File TPM_path
        File Count_path
        String OutputPrefix 
        Int Memory
    }
command <<<
    Rscript /tmp/identify_sample_outliers.R \
        --TPM_file ${TPM_path} \
        --count_file ${Count_path} \
        --prefix ${OutputPrefix}
    >>>

runtime {
        docker: 'evinpadhi/transcriptomic_outlier_qc:latest'        
        memory: "${Memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        cpu: "1"
        zones: ["us-central1-c"]
    }

output {
    File Outliers = "${OutputPrefix}_connectivity_outliers.tsv"
    File Zscores = "${OutputPrefix}_connectivity_scores.tsv"

    }
}

workflow transcriptomic_outliers_QC {
    call call_outliers 
    }
