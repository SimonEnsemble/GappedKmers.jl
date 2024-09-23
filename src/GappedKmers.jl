module GappedKmers

using BioSequences, IterTools, CairoMakie

include("list_of_gapped_kmers.jl")
include("featurizer.jl")
include("kernel.jl")

export list_of_gapped_kmers, convert_to_regex, string_to_DNA_seq, number_of_gapped_kmers, # list_of_gapped_kmers.jl
       featurizer, build_feature_matrix, # featurizer.jl
       gapped_kmer_kernel # kernel.jl
end
