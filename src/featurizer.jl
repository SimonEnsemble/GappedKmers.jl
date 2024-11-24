"""
    x = gkmer_feature(seq, ℓ, k)
    x = _gkmer_feature(seq, ℓ, k, gkmers) # pass list of gapped k-mers for speed

constructs the gapped k-mer feature vector of a DNA sequence.
this vector lists the number of occurences of each possible gapped k-mer in the DNA sequence.
order follows that of `list_of_gapped_kmers` if the `gkmers` list is not provided.
i.e. `x[i]` is count of gapped k-mer `list_of_gapped_kmers(ℓ, k)[i]` in `seq`.

#### arguments
* `ℓ::Int`: length of subsequence in gkmer
* `k::Int`: number of informative (non-wildcard) positions (nucleotides) in gkmer
* `seq::LongDNA`: DNA sequence we wish to featurize
* `gkmers::Vector{Regex{DNA}}`: (optional) list of gapped k-mers to search for (gives corresponding entries)
"""
function gkmer_feature(seq::LongDNA, ℓ::Int, k::Int)
	# get list of gapped k-mers
    gkmers = convert_to_regex.(
        list_of_gapped_kmers(ℓ, k) # gives vector of strings
    )
    return _gkmer_feature(seq, ℓ, k, gkmers)
end

function _gkmer_feature(seq::LongDNA, ℓ::Int, k::Int, gkmers::Vector{BioSequences.RE.Regex{DNA}})
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
	feature_matrix = gkmer_feature_matrix(DNA_seqs, ℓ, k)

creates matrix of feature vectors of a list of DNA sequences, with the feature vectors in the columns.
* rows: gkmer's
* columns: DNA sequences
* entries (integers): counts of gkmer's in the DNA sequences

#### arguments
* `DNA_seqs::Vector{LongSequence{DNAAlphabet{2}}}`: list of DNA sequences we wish to featurize
* `ℓ::Int`: length of subsequence of gkmers
* `k::Int`: number of informative (non-wildcard) positions (nucleotides) of gkmers

#### returns
* `feature_matrix::Array{Int}`
"""
function gkmer_feature_matrix(DNA_seqs::Vector{LongSequence{DNAAlphabet{2}}}, ℓ::Int, k::Int)
	# get list of gapped k-mers. pre-compute for speed.
    gkmers = convert_to_regex.(
        list_of_gapped_kmers(ℓ, k) # gives vector of strings
    )

    feature_matrix = zeros(Int, length(gkmers), length(DNA_seqs))
	for s = 1:length(DNA_seqs)
		feature_matrix[:, s] = _gkmer_feature(DNA_seqs[s], ℓ, k, gkmers)
	end

	return feature_matrix
end

"""
    gkmer_feature_info(seqs, ℓ, k, include_absent_gkmers=false)

return gapped-kmer info about a list of DNA sequences, in the form of a data frame.

#### arguments
* `seqs::Vector{LongDNA{2}}`: list of DNA sequences to featurize
* `ℓ::Int`: length of subsequence of gkmers
* `k::Int`: number of informative (non-wildcard) positions (nucleotides) of gkmers
* `include_absent_gkmers::Bool=false`: include rows where no sequence had that gkmer
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
        insertcols!(data, String(seq) => _gkmer_feature(seq, ℓ, k, gkmers))
    end
    
    if include_absent_gkmers
        return data
    else
        ids_keep = sum(Matrix(data)[:, 2:end], dims=2)[:] .!= 0
        return data[ids_keep, :]
    end
end
