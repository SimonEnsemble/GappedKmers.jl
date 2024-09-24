"""
    gkmer_feature_info(seqs, ℓ, k, include_absent_gkmers=false)

return gapped-kmer info about a list of DNA sequences, in the form of a data frame.

# arguments
* `seqs::Vector{LongDNA{2}}`: list of DNA sequences
* `ℓ::Int`: length of subsequence
* `k::Int`: number of informative (non-wildcard) positions (nucleotides)
* `include_absent_gkmers::Bool`: include rows where no sequence had that gapped k-mer.
"""
function gkmer_feature_info(seqs::Vector{LongDNA{2}}, ℓ::Int, k::Int, include_absent_gkmers::Bool=false)
    # get list of gapped k-mers
    gkmers_strings = list_of_gapped_kmers(ℓ, k)
    gkmers = convert_to_regex.(gkmers_strings)
    
    # store in data frame
    data = DataFrame(
        "gapped k-mer" => gkmers_strings,
    )

    # compute features; append to data frame
    for seq in seqs
        insertcols!(data, String(seq) => _featurizer(seq, ℓ, k, gkmers))
    end
    
    if include_absent_gkmers
        return data
    else
        ids_keep = sum(Matrix(data)[:, 2:end], dims=2)[:] .!= 0
        return data[ids_keep, :]
    end
end
