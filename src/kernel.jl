using LinearAlgebra

"""
    gapped_kmer_kernel(s₁, s₂, ℓ, k)

compute the gapped k-mer kernel between two input DNA sequences, which counts the number of pairings of gapped k-mers between the two input sequences.

* `s₁::LongDNA`: sequence 1
* `s₂::LongDNA`: sequence 2
* `ℓ::Int`: length of subsequence
* `k::Int`: number of informative (non-wildcard) positions (nucleotides)
"""
function gapped_kmer_kernel(s₁::LongDNA, s₂::LongDNA, ℓ::Int, k::Int)
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

gapped_kmer_kernel(s₁::String, s₂::String, ℓ::Int, k::Int) = gapped_kmer_kernel(
    string_to_DNA_seq(s₁), string_to_DNA_seq(s₂), ℓ, k
)

"""
    gapped_kmer_kernel_matrix(seqs₁, seqs₂, k, ℓ, normalize=true, symmetric=false)
    gapped_kmer_kernel_matrix(seqs, k, ℓ, normalize=true) # symmetric

creates Gram matrix giving kernel value between every pair of DNA sequences between the two input lists of DNA sequences.

# arguments
* `seqs₁::Vector{LongDNA{2}}` list of sequences for rows
* `seqs₂::Vector{LongDNA{2}`: list of sequences for columns
* `ℓ::Int`: length of l-mer
* `k::Int`: number of informative positions (nucleotides)
* `normalize::Bool`: normalize the kernel matrix
* `symmetric::Bool`: assume symmetric (saves computation)

# returns
`kernel_matrix::Matrix{Float64}`: where `kernel_matrix[i, j]` is gapped k-mer kernel between `seqs₁[i]` and `seqs₂[j]`.
"""
function gapped_kmer_kernel_matrix(
    seqs₁::Vector{LongDNA{2}},
    seqs₂::Vector{LongDNA{2}},
    ℓ::Int,
    k::Int;
    normalize::Bool=true, 
    symmetric::Bool=false
)
    # initialize Gram matrix
    matrix = zeros((length(seqs₁), length(seqs₂)))
    for (i, s₁) in enumerate(seqs₁)
        if ! symmetric # for normalization
            kᵢᵢ = gapped_kmer_kernel(s₁, s₁, ℓ, k)
        end
        for (j, s₂) in enumerate(seqs₂)
            if ! symmetric # for normalization
                kⱼⱼ = gapped_kmer_kernel(s₂, s₂, ℓ, k)
            end
            
            if symmetric && (i > j)
                continue # avoid computing twice
            end
            
            matrix[i, j] = gapped_kmer_kernel(s₁, s₂, ℓ, k)
            
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
    
    return matrix
end

gapped_kmer_kernel_matrix(seqs::Vector{LongDNA{2}}, ℓ::Int, k::Int; normalize::Bool=true) = gapped_kmer_kernel_matrix(seqs, seqs, ℓ, k, normalize=normalize, symmetric=true)
