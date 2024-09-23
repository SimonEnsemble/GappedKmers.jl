module GappedKmers

using BioSequences, IterTools, CairoMakie, NPZ, JLD

include("list_of_gapped_kmers.jl")

export list_of_gapped_kmers, convert_to_regex, string_to_DNA_seq, number_of_gapped_kmers, featurizer, build_feature_matrix, feature_matrix_file_direct, feature_matrix_import, gapped_kmer_kernel, gapped_kmer_kernel_matrix, kernel_matrix_file_name, kernel_matrix_import, viz_kernel_matrix 

end
