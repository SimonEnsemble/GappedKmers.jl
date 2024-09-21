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
