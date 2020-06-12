#!/bin/env nextflow
dataset_name = params.dataset.split("/")[-1].split("\\.")[0]
familyset_name = params.familyset.split("/")[-1].split("\\.")[0]

process transform_familyset{
    publishDir "$dataset_name/$familyset_name/intermediate/", mode: 'copy', overwrite: 'true'
    input:
    path familyset from params.familyset
    path dataset from params.dataset
    output:
    path "familyset.csv" into familyset_csv
    """
    transform_familyset.py $familyset $dataset
    """
}

process dichotomise{
    publishDir "$dataset_name/$familyset_name/intermediate/", mode: 'copy', overwrite: 'true'
    input:
    path familyset from familyset_csv
    path dataset from params.dataset
    val type from params.type
    output:
    path "dichotomised_genes.csv" into dichotomised_dataset
    path "family_thresholds.csv" into family_thresholds
    """
    dichotomise_families.py $familyset $dataset $type
    """
}

process calculate_IC{
    publishDir "$dataset_name/$familyset_name/results/", mode: 'copy', overwrite: 'true'
    input:
    path familyset from familyset_csv
    path dataset from dichotomised_dataset
    val type from params.type
    output:
    path "family_IC.csv" into family_IC
    """
    calculate_family_IC.py $familyset $dataset
    """
}

// family_sizes.splitText().println()

process permute_families{
    input:
    path familyset from familyset_csv
    path dataset from params.dataset
    each n from Channel.of(1..params.shuffle_N)
    output:
    path "shuffled_families.csv" into shuffled_familisets
    """
    permute_families.py $familyset $dataset
    """
}

process calculate_shuffled_IC{
    input:
    each familyset from shuffled_familisets
    path dataset from dichotomised_dataset
    val type from params.type
    output:
    path "family_IC.csv" into shuffled_family_IC
    """
    calculate_family_IC.py $familyset $dataset
    """
}

process collect_shuffled_IC{
    publishDir "$dataset_name/$familyset_name/results/", mode: 'copy', overwrite: 'true'
    input:
    val shuffled_family_ICs from shuffled_family_IC.collect()
    output:
    path "shuffled_family_IC.csv" into shuffled_family_IC_combined
    """
    collect_shuffled_family_IC.py "$shuffled_family_ICs"
    """
}

// process calculate_shuffled_IC_minmax{
//     // publishDir "$dataset_name/$familyset_name/intermediate/", mode: 'copy', overwrite: 'true'
//     input:
//     path familyset from familyset_csv
//     path dataset from params.dataset
//     each size from family_sizes.splitText()
//     output:
//     path "permutations.csv" into permutations
//     """
//     calculate_family_permutations.py $familyset $dataset --minmax $size 
//     """
// }

// process collect_shuffled_IC_minmax{
//     publishDir "$dataset_name/$familyset_name/intermediate/minmax/", mode: 'copy', overwrite: 'true'
//     input:
//     val shuffled_family_ICs from permutations.collect()
//     output:
//     path "shuffled_family_IC.csv" into shuffled_family_IC
//     """
//     collect_shuffled_family_IC.py "$shuffled_family_ICs"
//     """
// }

// process calculate_critical_values_minmax{
//     publishDir "$dataset_name/$familyset_name/results/minmax/", mode: 'copy', overwrite: 'true'
//     input:
//     path shuffled_family_IC from shuffled_family_IC
//     output:
//     path "critical_values.csv" into critical_values
//     """
//     calculate_family_critical_values.py $shuffled_family_IC --minmax
//     """ 
// }

// process calculate_family_pvalues_minmax{
//     publishDir "$dataset_name/$familyset_name/results/minmax/", mode: 'copy', overwrite: 'true'
//     input:
//     path shuffled_family_IC from shuffled_family_IC
//     path family_IC from family_IC
//     output:
//     path "family_IC_pvalues.csv" into family_IC_pvalues
//     """
//     calculate_family_pvalues.py $shuffled_family_IC $family_IC --minmax
//     """ 
// }

// process calculate_shuffled_IC_single{
//     // publishDir "$dataset_name/$familyset_name/intermediate/", mode: 'copy', overwrite: 'true'
//     input:
//     path familyset from familyset_csv
//     path dataset from params.dataset
//     each size from family_sizes.splitText()
//     output:
//     path "permutations.csv" into permutations_single
//     """
//     calculate_family_permutations.py $familyset $dataset $size
//     """
// }

// process collect_shuffled_IC_single{
//     publishDir "$dataset_name/$familyset_name/intermediate/single/", mode: 'copy', overwrite: 'true'
//     input:
//     val shuffled_family_ICs from permutations_single.collect()
//     output:
//     path "shuffled_family_IC.csv" into shuffled_family_IC_single
//     """
//     collect_shuffled_family_IC.py "$shuffled_family_ICs"
//     """
// }

// process calculate_critical_values_single{
//     publishDir "$dataset_name/$familyset_name/results/single/", mode: 'copy', overwrite: 'true'
//     input:
//     path shuffled_family_IC from shuffled_family_IC_single
//     output:
//     path "critical_values.csv" into critical_values_single
//     """
//     calculate_family_critical_values.py $shuffled_family_IC
//     """ 
// }

// process calculate_family_pvalues_single{
//     publishDir "$dataset_name/$familyset_name/results/single/", mode: 'copy', overwrite: 'true'
//     input:
//     path shuffled_family_IC from shuffled_family_IC_single
//     path family_IC from family_IC
//     output:
//     path "family_IC_pvalues.csv" into family_IC_pvalues_single
//     """
//     calculate_family_pvalues.py $shuffled_family_IC $family_IC
//     """ 
// }