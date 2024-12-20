using LinearAlgebra

"""
    gapped_kmer_kernel(s₁, s₂, ℓ, k, w=(i, j) -> 1.0)

compute the gapped k-mer kernel between two input DNA sequences.
this kernel counts the number of pairings of gapped k-mers between the two input sequences.

optionally, provide a [symmetric, positive] weight function `w` that weighs contributions
of gkmers in common according to their positions in the sequence.

#### arguments
* `s₁::LongDNA`: sequence 1
* `s₂::LongDNA`: sequence 2
* `ℓ::Int`: length of subsequence
* `k::Int`: number of informative (non-wildcard) positions (nucleotides)
* `w::Function`: (optional) weight of how each gkmer contributes based on location in the DNA sequence. 
"""
function gapped_kmer_kernel(s₁::LongDNA, s₂::LongDNA, ℓ::Int, k::Int; w::Function=(i, j) -> 1.0)
    @assert k ≤ ℓ
    
    kernel_value = 0.0
    
    # compute length of the two sequences
    n₁ = length(s₁)
    n₂ = length(s₂)
    
    # pre-compute contribution to kernel of l-mer pair based on # of matching positions
    # contribution[nb_matches - k + 1] gives contribution for nb_matches matches
    contribution = [1.0 * binomial(nb_matches, k) for nb_matches = k:ℓ]
    
    # loop over pairs of l-mers in the two sequences
    for i in 1:(n₁+1-ℓ) # i : starting position of l-mer in s₁
        lmer₁ = @view s₁[i:i+(ℓ-1)]
        for j in 1:(n₂+1-ℓ) # j : starting position of l-mer in s₂
            lmer₂ = @view s₂[j:j+(ℓ-1)]
            # count the number of matching nucleotides
            nb_matches = count(==, lmer₁, lmer₂)
            # this pair contributes if k or more matches
            if nb_matches >= k
                kernel_value += contribution[1 + nb_matches - k] * w(i, j)
            end
        end
    end
    
    # return # pairs of gapped-k-mers between the two input sequences
    return kernel_value
end

gapped_kmer_kernel(s₁::String, s₂::String, ℓ::Int, k::Int; w::Function=(i, j) -> 1.0) = gapped_kmer_kernel(
    string_to_DNA_seq(s₁), string_to_DNA_seq(s₂), ℓ, k, w=w
)

"""
    gapped_kmer_kernel_matrix(seqs₁, seqs₂, k, ℓ, normalize=true, symmetric=false)
    gapped_kmer_kernel_matrix(seqs, k, ℓ, normalize=true) # symmetric

computes kernel matrix giving kernel value between every pair of DNA sequences 
made between the two input lists of DNA sequences.

if one list of sequences is provided, compute the Gram matrix.

optionally, provide a [symmetric] weight function that weighs contributions
of gkmers in common according to their positions in the sequences.

#### arguments
* `seqs₁::Vector{LongDNA{2}}` list of sequences for rows
* `seqs₂::Vector{LongDNA{2}`: list of sequences for columns
* `ℓ::Int`: length of l-mer
* `k::Int`: number of informative positions (nucleotides)
* `normalize::Bool`: normalize the kernel matrix
* `symmetric::Bool`: assume symmetric (saves computation)
* `w::Function`: (optional) weight function dictating how each gkmer contributes based on location in the DNA sequence.

#### returns
`kernel_matrix::Matrix{Float64}`: where `kernel_matrix[i, j]` is gapped k-mer kernel between `seqs₁[i]` and `seqs₂[j]`.
"""
function gapped_kmer_kernel_matrix(
    seqs₁::Vector{LongDNA{2}},
    seqs₂::Vector{LongDNA{2}},
    ℓ::Int,
    k::Int;
    normalize::Bool=true, 
    symmetric::Bool=false,
    w::Function=(i, j) -> 1.0
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
            
            matrix[i, j] = gapped_kmer_kernel(s₁, s₂, ℓ, k, w=w)
            
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
