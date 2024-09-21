using Test

using BioSequences

push!(LOAD_PATH, joinpath(pwd(), "src")) # run from main dir
using GappedKmers


###
#   test list of gapped k-mers
###
@testset "testing list of gapped k-mers and conversions" begin
    k = 3
    ℓ = 6
    gkmers = list_of_gapped_kmers(ℓ, k)
    @test length(gkmers) == length(unique(gkmers))
    @test length(gkmers) == number_of_gapped_kmers(ℓ, k)
    @test "A-A--A" in gkmers
    
    gkmer = "A-T"
    gkmer_regex = convert_to_regex(gkmer)

    seq₁ = string_to_DNA_seq("AGTTTCG")
    seq₂ = string_to_DNA_seq("GGCGAAACCT")

    @test occursin(gkmer_regex, seq₁)
    @test ! occursin(gkmer_regex, seq₂)
end
