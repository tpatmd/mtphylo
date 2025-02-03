params.input_file = "input.fasta"

workflow {
    Channel
        .fromPath(params.input_file, type: "file", checkIfExists: true)
        .set { input_ch }
    splitCDS(input_ch)
}

process splitCDS {
    publishDir "01-SplitCDS", mode: "copy"

    input:
    path input_fasta

    output:
    path "*.fasta"
    path "split-report.tsv"

    script:
    """
    split.py ${params.input_file}
    """
}
