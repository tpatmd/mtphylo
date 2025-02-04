params.input_file = "input.fasta"
params.pcg = 13
// GBlocks parameter
params.gf = 20
params.trimmer = "bmge"
params.partition = "codon"
params.merge = true

workflow {
    Channel
        .fromPath(params.input_file, type: "file", checkIfExists: true)
        .set { input_ch }
    splitCDS(input_ch)
    multipleAlignment(splitCDS.out.output_cds.flatten(), splitCDS.out.num_species.collect())
    concatenateAlignments(multipleAlignment.out.alignment_fasta.collect())
    phylogeneticAnalysis(concatenateAlignments.out.cat_fasta, concatenateAlignments.out.cat_partition)
}

process splitCDS {
    publishDir "01-SplitCDS", mode: "copy"

    input:
    path input_fasta

    output:
    path "*.fasta", emit: output_cds
    path "split-report.tsv"
    stdout emit: num_species

    script:
    """
    split.py ${params.input_file}
    """
}

process multipleAlignment {
    publishDir "02-MultipleAlignment", mode: "copy"

    input:
    path input_cds
    val num_species

    output:
    path "*.nt_cleanali.fasta", emit: alignment_fasta
    path "*"

    script:
    if (params.trimmer == "gblocks") {
        """
        translatorx.pl -i ${input_cds} -o ${input_cds.baseName} -p F -c 5 -t F -g "-b2=${params.gf} -b4=5 -b5=h"
        """
    }
    else if (params.trimmer == "bmge") {
        """
        translatorx.pl -i ${input_cds} -o ${input_cds.baseName} -p F -c 5 -t F
        bmge -i ${input_cds.baseName}.nt_ali.fasta -of ${input_cds.baseName}.nt_cleanali.fasta -t CODON
        """
    }
}

process concatenateAlignments {
    publishDir "03-ConcatenateAlignments", mode: "copy"

    input:
    path alignment_files

    output:
    path "concatenated_dataset.fasta", emit: cat_fasta
    path "concatenated_dataset.partition", emit: cat_partition

    script:
    """
    cat ${alignment_files} > tmp.fasta
    concat.py tmp.fasta ${params.pcg} ${params.partition}
    """
}

process phylogeneticAnalysis {
    publishDir "04-PhylogeneticAnalysis", mode: "copy"

    input:
    path alignment_file
    path partition_file

    output:
    path "*"

    script:
    if (params.merge) {
        """
        iqtree2 -s ${alignment_file} -p ${partition_file} -m MFP+MERGE -B 1000 -T AUTO --merge greedy
        """
    }
    else {
        """
        iqtree2 -s ${alignment_file} -p ${partition_file} -m MFP -B 1000 -T AUTO
        """
    }
}
