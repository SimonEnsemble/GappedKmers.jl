"""
    number_of_gapped_kmers(3, 2) # ℓ = 3, k = 2

return the number of possible gapped k-mers.

## arguments
* `ℓ::Int`: length of subsequence.
* `k::Int`: number of informative (non-wildcard) positions (nucleotides)
"""
number_of_gapped_kmers(ℓ::Int, k::Int) = 4^k * binomial(ℓ, k)

"""
    list_of_gapped_kmers(3, 2) # ℓ = 3, k = 2

create and return list of all possible gapped k-mers.

## arguments
* `ℓ::Int`: length of subsequence.
* `k::Int`: number of informative (non-wildcard) positions (nucleotides)

## returns
* `all_gapped_kmers::Vector{String}` the vector of gapped k-mers. the wildcard positions are denoted with a "-".
"""
function list_of_gapped_kmers(ℓ::Int, k::Int)
    @assert k ≤ ℓ

	# list of gapped kmers
	gkmers = []
	# list of possible positions
	nucleotides = "ATGC-"
	
	# function found on https://discourse.julialang.org/t/how-to-create-all-possible-k-mers-with-a-set-of-alphabets/83186
	for gkmer in join.(Iterators.product(Iterators.repeated(nucleotides, ℓ)...))
		# if any gkmer produced has ℓ - k wildcards '-', append to list of kmers
		if count(x -> x == '-', gkmer) == (ℓ - k)
			push!(gkmers, gkmer)
		end
	end

    return gkmers
	
	# replace wildcard '-' with "[ATCG]" for BioSequences.jl
	gkmers = [replace(gkmer, s"-" => s"[ATCG]") for gkmer in gkmers]

	# convert to Regex expressions for searching
	return [BioSequences.RE.Regex{DNA}(gkmer) for gkmer in gkmers]
end

"""
    convert_to_regex("A-T")

convert a gapped k-mer in the form of a `String` with "-" representing a wildcard position, into a BioSequences regular expression for the purpose of searching DNA sequences for this gapped k-mer.

# example
```julia
using BioSequences

gkmer = "A-T"
gkmer_regex = convert_to_regex(gkmer)

seq = dna"AGTTTCG"

occursin(gkmer_regex, seq) # true
```
"""
function convert_to_regex(gkmer::String)
    new_gkmer = replace(gkmer, s"-" => s"[ATCG]")
    return BioSequences.RE.Regex{DNA}(new_gkmer)
end

"""
	string_to_DNA_seq("ATTTT")

converts a DNA sequence represented as a `String` to an efficent 2-bit-per-nucleotide representation via `BioSequences.jl`, specifically, a `LongDNA{2}`.
"""
string_to_DNA_seq(DNA_seq::String) = LongDNA{2}(DNA_seq)

"""
    featurizer(seq, k, ℓ)
outputs gapped k-mer feature vector which lists number of occurences of each possible gapped k-mer

##arguments
* `k::Int`: number of informative positions (nucleotides)
* `ℓ::Int`: length of l-mer
* `seq::LongDNA`: DNA sequence
"""
function featurizer(seq::LongDNA, k::Int, ℓ::Int)
	n = length(seq)
	
	# get list of gapped k-mers
	# TODO: pre-compute this outside for speed, if we need it?
    gkmers = list_of_gapped_kmers(k, ℓ)

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
	build_feature_matrix(DNA_seqs, k, ℓ)
creates matrix of feature vectors with each column corresponding to each possible gapped k-mer and each row corresponding to DNA sequence being evaluated
##arguments
* `DNA_seqs::Vector{LongSequence{DNAAlphabet{2}}}`: list of sequences for rows
* `k::Int`: number of informative positions (nucleotides)
* `ℓ::Int`: length of l-mer
"""
function build_feature_matrix(
	DNA_seqs::Vector{LongSequence{DNAAlphabet{2}}},
	k::Int,
	ℓ::Int;
	file_save::Bool=false
)
	seq_feature_vector = featurizer(DNA_seqs[1], k, ℓ)
	matrix = zeros(length(seq_feature_vector), length(DNA_seqs))
	for s = 1:length(DNA_seqs)
		matrix[:, s] = featurizer(DNA_seqs[s], k, ℓ)
	end
	if file_save
		object_name = feature_matrix_file_direct(k, ℓ)
		save(object_name, "matrix", matrix)
	end
	return matrix
end

"""
	feature_matrix_file_direct(k, ℓ)
provides file directory with specified filename corresponding to k and ℓ values
"""
function feature_matrix_file_direct(k::Int, ℓ::Int)
	prefix = "data/feature matrices/feature_matrix_k_"
	return string(prefix, string(k), "_l_", string(ℓ), ".jld")
end

"""
	feature_matrix_import(k, ℓ)
imports feature matrix saved via the build_feature_matrix() function
"""
function feature_matrix_import(k::Int, ℓ::Int)
	object_name = feature_matrix_file_direct(k, ℓ)
	return load(object_name)["matrix"]
end

"""
	gapped_kmer_kernel(s₁, s₂, k, ℓ)
compute the un-normalized gapped k-mer kernel between two input DNA sequences.
this counts the number of pairings of gapped k-mers between the two input sequences.

* `s₁::LongDNA`: sequence 1
* `s₂::LongDNA`: sequence 2
* `k::Int`: number of informative positions (number of nucleotides)
* `ℓ::Int`: length of l-mer
"""
function gapped_kmer_kernel(s₁::LongDNA, s₂::LongDNA, k::Int, ℓ::Int)
	@assert k ≤ ℓ
	
	# compute length of the two sequences
	n₁ = length(s₁)
	n₂ = length(s₂)
	
	# stores counts of number # l-mer pairs with m mis-matching nucleotides.
	#    nb_lmer_pairs_w_m_mismatches[m + 1]: # l-mer pairs with m mismatches
	#    (can be anywhere for 0 to ℓ mismatching nucleotides...)
	nb_lmer_pairs_w_m_mismatches = [0 for m = 0:ℓ]
	
	# loop over pairs of l-mers in the two sequences
	for i in 1:(n₁+1-ℓ) # i : starting position of l-mer in s₁
		lmer₁ = s₁[i:i+(ℓ-1)]
		for j in 1:(n₂+1-ℓ) # j : starting position of l-mer in s₂
			lmer₂ = s₂[j:j+(ℓ-1)]
			m = count(!=, lmer₁, lmer₂) # count the number of mismatches
			nb_lmer_pairs_w_m_mismatches[m+1] += 1 # record finding
		end
	end
	
	# compute kernel
	kernel_value = 0 # an integer
	for m = 0:(ℓ-k) # m ≤ (ℓ - k) contributes to kernel
		# TODO: can make faster by pre-computing binomial(ℓ-m, k)
		kernel_value += nb_lmer_pairs_w_m_mismatches[m+1] * binomial(ℓ-m, k)
	end

	# return # pairs of gapped-k-mers between the two input sequences
	return kernel_value
end

"""
	gapped_kmer_kernel_matrix(seqs₁, seqs₂, k, ℓ)

creates Gram matrix of kernel results through looping over pairs of DNA sequences in our two arrays of sequences passed.

* `seqs₁::Vector{LongSequence{DNAAlphabet{2}}}`: list of sequences for rows
* `seqs₂::Vector{LongSequence{DNAAlphabet{2}}}`: list of sequences for columns
* `k::Int`: number of informative positions (nucleotides)
* `ℓ::Int`: length of l-mer
"""
function gapped_kmer_kernel_matrix(
	seqs₁::Vector{LongSequence{DNAAlphabet{2}}}, seqs₂::Vector{LongSequence{DNAAlphabet{2}}}, 
	k::Int,
	ℓ::Int;
	normalize::Bool=true, 
	file_save::Bool=false,
	symmetric::Bool=false
)
	# initialize Gram matrix
	matrix = zeros((length(seqs₁), length(seqs₂)))
	for (i, s₁) in enumerate(seqs₁)
		if ! symmetric # for normalization
			kᵢᵢ = gapped_kmer_kernel(s₁, s₁, k, ℓ)
		end
		for (j, s₂) in enumerate(seqs₂)
			if ! symmetric # for normalization
				kⱼⱼ = gapped_kmer_kernel(s₂, s₂, k, ℓ)
			end
			
			if symmetric && (i > j)
				continue # avoid computing twice
			end
			
			matrix[i, j] = gapped_kmer_kernel(s₁, s₂, k, ℓ)
			
			if symmetric
				matrix[j, i] = matrix[i, j] # account for symmetry
			else
				if normalize
					matrix[i, j] = matrix[i, j] / (sqrt(kᵢᵢ * kⱼⱼ))
				end
			end
		end
	end
	
	if normalize && symmetric
		# kᵢᵢ's are in the diagonal...
		D = diagm(1.0 ./ sqrt.(diag(matrix)))
		matrix = D * matrix * D
	end
	
	if file_save
		object_name = kernel_matrix_file_name(k, ℓ, normalize)
		npzwrite(object_name, matrix)
	end
	
	return matrix
end

"""
	kernel_matrix_file_name(k, ℓ, normalize)

returns filename (a convention) for storing the kernel matrix.
"""
function kernel_matrix_file_name(k::Int, ℓ::Int, normalize::Bool)
	if normalize
		prefix = "normalized_kernel_matrix_k_"
	else
		prefix = "kernel_matrix_k_"
	end
	return string(prefix, string(k), "_l_", string(ℓ), ".npy")
end

"""
	kernel_matrix_import(k, ℓ)
imports kernel matrix saved via gapped_kmer_kernel_matrix() function
"""
function kernel_matrix_import(k::Int, ℓ::Int; normalize::Bool=true)
	return npzread(joinpath("data", "kernel matrices", kernel_matrix_file_name(k, ℓ, normalize)))
end

"""
	viz_kernel_matrix(K)
outputs heatmap visualization of gapped kmer kernel matrix
"""
function viz_kernel_matrix(K::Matrix{Float64})
	fig, ax, hm = heatmap(K; colormap = :turku10,  
		axis=(; aspect=DataAspect(), xlabel="DNA seq.", ylabel="DNA seq."),
		colorrange=(0.0, 1.0)
	)
	Colorbar(fig[:, end+1], hm, label="similarity score")
	return fig
end