using BioSequences

"""
    ℓ = 3 # length of subsequence
    k = 2 # number of informative (non-wildcard) positions (nucleotides)
    number_of_gapped_kmers(ℓ, k)

returns the number (`Int`) of possible gapped k-mers.
"""
number_of_gapped_kmers(ℓ::Int, k::Int) = 4 ^ k * binomial(ℓ, k)

"""
    ℓ = 3 # length of subsequence
    k = 2 # number of informative (non-wildcard) positions (nucleotides)
    list_of_gapped_kmers(ℓ, k)

create and return a vector (`Vector{String}`) of all possible gapped k-mers.
the wildcard position in a gapped k-mer is represented by a "-" character.
"""
function list_of_gapped_kmers(ℓ::Int, k::Int)
    @assert k ≤ ℓ

    # list of possible nucleotides, incl. wildcard
    nucleotides = "ATGC-"
    
    # https://discourse.julialang.org/t/how-to-create-all-possible-k-mers-with-a-set-of-alphabets/83186
    # list of gapped kmers
    gkmers = String[]
    for gkmer in join.(Iterators.product(Iterators.repeated(nucleotides, ℓ)...))
        # if any gkmer produced has ℓ - k wildcards '-', append to list of kmers
        if count(x -> x == '-', gkmer) == (ℓ - k)
            push!(gkmers, gkmer)
        end
    end

    return gkmers
end

"""
    convert_to_regex("A-T") # returns a Regex{DNA}

converts a gapped k-mer in the form of a `String`, with the character "-" representing 
a wildcard position, into a `BioSequences` regular expression `Regex{DNA}` for the purpose 
of searching DNA sequences for this gapped k-mer.

### example
```julia
using BioSequences

gkmer = "A-T"
gkmer_regex = convert_to_regex(gkmer) # gives a Regex{DNA}

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

converts a DNA sequence represented as a `String` to an efficent 2-bit-per-nucleotide 
representation via `BioSequences.jl`---specifically, a `LongDNA{2}`.
"""
string_to_DNA_seq(DNA_seq::String) = LongDNA{2}(DNA_seq)
