module GappedKmers

using BioSequences, IterTools

include("list_of_gapped_kmers.jl")

export list_of_gapped_kmers, convert_to_regex, string_to_DNA_seq, number_of_gapped_kmers

end
