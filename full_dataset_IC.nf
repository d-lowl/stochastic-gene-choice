#!/bin/env nextflow
dataset_name = params.dataset.split("/")[-1].split("\\.")[0]
process filter_chromosome {
    echo true
    publishDir "$dataset_name/intermediate/", mode: 'copy', overwrite: 'true'
    input:
    path gff from params.gff
    output:
    file "chr*_filtered.csv" into filtered_chromosome

    """
    split_chromosomes.py $gff
    """
}

process get_chromosome_names {
    input:
    each chr_file from filtered_chromosome
    output:
    tuple stdout, chr_file into _named_filtered_chromosome

    """
    #!/bin/env python
    print("$chr_file".split("/")[-1].split("_")[0])
    """
}

(chromosomes, chromosomes2) = _named_filtered_chromosome.into(2)
process calculate_thresholds{
    publishDir "$dataset_name/intermediate/$type/", mode: 'copy', overwrite: 'true'
    input:
    path genes from params.dataset
    val type from params.type
    output:
    path "optimal_thresholds.csv" into optimal_thresholds
    path "dichotomised_genes.csv" into dichotomised_genes

    """
    dichotomise.py $genes dichotomised_genes.csv optimal_thresholds.csv $type
    """
}

process filter_genes {
    echo true
    publishDir "$dataset_name/intermediate/$type/", mode: 'copy', overwrite: 'true'
    input:
    path all_genes from params.dataset
    path dichotomised_genes from dichotomised_genes
    val type from params.type
    output:
    file "filtered_dichotomised_genes.csv" into filtered_dichotomised_genes

    """
    filter_genes.py $all_genes $dichotomised_genes filtered_dichotomised_genes.csv $type
    """
}

process calculate_IC{
    publishDir "$dataset_name/results/$type/", mode: 'copy', overwrite: 'true'
    input:
    tuple chr_name, chromosome from chromosomes
    path dichotomised_genes from filtered_dichotomised_genes
    val type from params.type
    output:
    tuple chr_name, "*.csv" into _named_stage1_IC

    """
    calculate_IC.py $dichotomised_genes $chromosome $chr_name
    """    
}

process calculate_shuffled_IC{
    input:
    val type from params.type
    path dichotomised_genes from filtered_dichotomised_genes
    val chromosomes from filtered_chromosome.collect()
    each n from Channel.of(1..params.shuffle_N)
    output:
    file "shuffled_ic.csv" into shuffled_IC

    """
    calculate_chromosomes_shuffled_IC.py $dichotomised_genes "$chromosomes"
    """
}

process collect_shuffled_IC{
    publishDir "$dataset_name/results/$type/", mode: 'copy', overwrite: 'true'
    input:
    val shuffled_ICs from shuffled_IC.collect()
    val type from params.type
    output:
    path "shuffled_IC.csv" into all_shuffled_IC
    """
    collect_chromosomes_shuffled_IC.py "$shuffled_ICs"
    """
}

(named_stage1_IC, named_stage1_IC_copy) = _named_stage1_IC.into(2)
chromosome_IC_join = chromosomes2.join(named_stage1_IC_copy)

// process calculate_bootstrapped_IC{
//     publishDir "$dataset_name/results/$type/"
//     input:
//     tuple chr_name, chromosome, stage1_IC from chromosome_IC_join
//     path dichotomised_genes from filtered_dichotomised_genes
//     val type from params.type
//     output:
//     tuple chr_name, "bootstrapped_*_IC.csv" into named_bootstrapped_IC
//     """
//     calculate_IC_with_bootstrapping.py $dichotomised_genes $stage1_IC $chromosome $chr_name
//     """
// }

// process stitch_exlusive_regions{
//     publishDir "$dataset_name/results/$type/"
//     input:
//     tuple chr_name, bootstrapped_IC from named_bootstrapped_IC
//     val type from params.type
//     output:
//     path "exclusive_regions_*.csv" into exclusive_regions
//     """
//     stitch_exclusive_regions.py $bootstrapped_IC $chr_name
//     """
// }

// process ctcf_comparison{
//     echo true
//     publishDir "$dataset_name/results/$type/"
//     input:
//     val chromosomes from filtered_chromosome.collect()
//     val exclusive_regions from exclusive_regions.collect()
//     path ctcf_neuronal from params.ctcf_neuronal
//     path ctcf_stochastic from params.ctcf_stochastic
//     val type from params.type
//     output:
//     file "*.csv" into ctcf_tests
//     """
//     ctcf_comparison.py "$chromosomes" "$exclusive_regions" $ctcf_neuronal $ctcf_stochastic 
//     """
// }

// process ga {
//     publishDir "$dataset_name/results/$type/ga/full/"
//     input:
//     // each chromosome from filtered_chromosome
//     val type from params.type
//     path dichotomised_genes from filtered_dichotomised_genes
//     output:
//     path "ga_results*" into ga_results
//     path "ga_log*.csv" into ga_log

//     """
//     ga.py $dichotomised_genes 80 0.0005 0.8 uniform
//     """
// }