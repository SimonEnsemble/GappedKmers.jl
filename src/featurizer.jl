"""
    featurizer(seq, ℓ, k)
    featurizer(seq, ℓ, k, gkmers) # pass list of gapped k-mers for speed

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
    return featurizer(seq, ℓ, k, gkmers)
end

function featurizer(seq::LongDNA, ℓ::Int, k::Int, gkmers::Vector{BioSequences.RE.Regex{DNA}})
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
* `k::Int`: number of informative positions (nucleotides)
* `ℓ::Int`: length of l-mer

## returns
* `feature_matrix::Array{Int}`
"""
function build_feature_matrix(DNA_seqs::Vector{LongSequence{DNAAlphabet{2}}}, ℓ::Int, k::Int)
	# get list of gapped k-mers
    gkmers = convert_to_regex.(
        list_of_gapped_kmers(ℓ, k) # gives vector of strings
    )

    n_sequences = length(DNA_seqs)
    n_features = length(gkmers)

    feature_matrix = zeros(Int, n_features, n_sequences)
	for s = 1:n_sequences
		feature_matrix[:, s] = featurizer(DNA_seqs[s], ℓ, k, gkmers)
	end

	return feature_matrix
end
