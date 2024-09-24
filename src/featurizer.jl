"""
    featurizer(seq, ℓ, k)
    _featurizer(seq, ℓ, k, gkmers) # pass list of gapped k-mers for speed

outputs gapped k-mer feature vector of a DNA sequence, 
which lists the number of occurences of each possible gapped k-mer in the DNA sequence.
order follows that of `list_of_gapped_kmers` if the `gkmers` list is not provided.

## arguments
* `ℓ::Int`: length of subsequence
* `k::Int`: number of informative (non-wildcard) positions (nucleotides)
* `seq::LongDNA`: DNA sequence
* `gkmers::Vector{Regex{DNA}}`: list of gapped k-mers to search for (gives corresponding entries)
"""
function featurizer(seq::LongDNA, ℓ::Int, k::Int)
	# get list of gapped k-mers
    gkmers = convert_to_regex.(
        list_of_gapped_kmers(ℓ, k) # gives vector of strings
    )
    return _featurizer(seq, ℓ, k, gkmers)
end

function _featurizer(seq::LongDNA, ℓ::Int, k::Int, gkmers::Vector{BioSequences.RE.Regex{DNA}})
    @assert ℓ ≥ k
    n = length(seq)

	# feature_vector[i] counts # of occurances of gapped-kmer i
    feature_vector = zeros(Int, length(gkmers))
	
	# iterate over list of gapped kmers
    for (gk, gkmer) in enumerate(gkmers)
        for i in 1:(n + 1 - ℓ) # i: starting index of ℓ-mer in seq
            lmer = seq[i:i + (ℓ - 1)]
			# does this lmer contain the gapped-kmer?
			if occursin(gkmer, lmer)
				# we found this gapped k-mer!
				feature_vector[gk] += 1
			end
		end
    end
	
    # output feature vector
    return feature_vector
end

"""
	build_feature_matrix(DNA_seqs, ℓ, k)

creates matrix of feature vectors of a list of DNA sequences, with the feature vectors in the columns.
* rows" gapped k-mer's
* cols: DNA sequences
* entries (integers): counts

## arguments
* `DNA_seqs::Vector{LongSequence{DNAAlphabet{2}}}`: list of sequences for rows
* `ℓ::Int`: length of subsequence
* `k::Int`: number of informative (non-wildcard) positions (nucleotides)

## returns
* `feature_matrix::Array{Int}`
"""
function build_feature_matrix(DNA_seqs::Vector{LongSequence{DNAAlphabet{2}}}, ℓ::Int, k::Int)
	# get list of gapped k-mers. pre-compute for speed.
    gkmers = convert_to_regex.(
        list_of_gapped_kmers(ℓ, k) # gives vector of strings
    )

    n_sequences = length(DNA_seqs)
    n_features = length(gkmers)

    feature_matrix = zeros(Int, n_features, n_sequences)
	for s = 1:n_sequences
		feature_matrix[:, s] = _featurizer(DNA_seqs[s], ℓ, k, gkmers)
	end

	return feature_matrix
end

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
