using Test

using BioSequences, LinearAlgebra

push!(LOAD_PATH, joinpath(pwd(), "src")) # run from main dir
using GappedKmers

@testset "testing list of gapped k-mers and conversions" begin
    k = 3
    ℓ = 6
    gkmers = list_of_gapped_kmers(ℓ, k)
    @test length(gkmers) == length(unique(gkmers))
    @test length(gkmers) == number_of_gapped_kmers(ℓ, k)
    @test "A-A--A" in gkmers
    @test "T--T-G" in gkmers
    
    gkmer = "A-T"
    gkmer_regex = convert_to_regex(gkmer)

    seq₁ = string_to_DNA_seq("AGTTTCG")
    seq₂ = string_to_DNA_seq("GGCGAAACCT")

    @test occursin(gkmer_regex, seq₁)
    @test ! occursin(gkmer_regex, seq₂)
end

@testset "testing gapped-kmer featurizer" begin
    @test featurizer(dna"TCGGAG", 1, 1) == [1, 1, 3, 1] # [#A, #T, #G, #C]
	@test sum(featurizer(dna"ATCGG", 2, 2)) == 4
	@test sum(featurizer(dna"ATCGG", 2, 1)) == 4 * 2
    
    ℓ = 3
    k = 2
    my_gkmers = convert_to_regex.(["-AT", "T-T"])
	@test sum(GappedKmers._featurizer(dna"AAAAAAAA", ℓ, k, my_gkmers)) == 0
	@test length(GappedKmers._featurizer(dna"AAAAAAAA", ℓ, k, my_gkmers)) == length(my_gkmers)
	@test GappedKmers._featurizer(dna"TTT", ℓ, k, my_gkmers) == [0, 1]
	@test GappedKmers._featurizer(dna"TTTAT", ℓ, k, my_gkmers) == [1, 2]

    # from Kyrstin's slides
    gkmers = list_of_gapped_kmers(3, 2)
    x = featurizer(dna"ATTCGGTGC", 3, 2)
    @test x[gkmers .== "T-C"] == [2]
end

@testset "testing gapped-kmer kernel" begin
    ### test 1
    s₁ = string_to_DNA_seq("ATT")
	s₂ = string_to_DNA_seq("ATC")

	@test gapped_kmer_kernel(s₁, s₂, 1, 1) == 1 + 2
	@test gapped_kmer_kernel(s₁, s₂, 2, 2) == 1 # one match, AT
	@test gapped_kmer_kernel(s₁, s₂, 2, 1) == 4 # A*, T*, *T x2

	@test gapped_kmer_kernel(s₁, s₂, 3, 1) == 2 # two matches A**, *T*
	@test gapped_kmer_kernel(s₁, s₂, 3, 2) == 1 # one match, AT*
	@test gapped_kmer_kernel(s₁, s₂, 3, 3) == 0 # no full matches
    
    # symmetry
    @test gapped_kmer_kernel(s₁, s₂, 3, 1) == gapped_kmer_kernel(s₂, s₁, 3, 1)

    ### test 2
    s₁ = string_to_DNA_seq("TCGGAG")
	s₂ = string_to_DNA_seq("GTAC")
	
	@test gapped_kmer_kernel(s₁, s₂, 2, 1) == 6
	
	@test gapped_kmer_kernel(s₁, s₂, 3, 3) == 0
	@test gapped_kmer_kernel(s₁, s₂, 3, 2) == 1 # just G*A
	
    # symmetry
    @test gapped_kmer_kernel(s₁, s₂, 2, 1) == gapped_kmer_kernel(s₂, s₁, 2, 1)

    ### test 3
    s₁ = string_to_DNA_seq("TTTAAACC")
    s₂ = string_to_DNA_seq("TTTAATCC")
	@test gapped_kmer_kernel(s₁, s₂, 8, 8) == 0
	@test gapped_kmer_kernel(s₁, s₂, 8, 7) == 1

    ### weighted
    s₁ = string_to_DNA_seq("ATTCG")
    s₂ = string_to_DNA_seq("ATTGG")
    
    w(i) = i <= 2 ? 1.0 : 0.0
    @test gapped_kmer_kernel(s₁, s₂, 2, 2, w=(i, j) -> w(i) * w(j)) == 2
    
    w(i, j) = (i - j) == 0 ? 1.0 : 0.0
    @test gapped_kmer_kernel(s₁, s₂, 1, 1, w=w) == 4 
    @test gapped_kmer_kernel(s₁, s₂, 2, 2, w=w) == 2
    @test gapped_kmer_kernel(s₁, s₂, 3, 3, w=w) == 1
    @test gapped_kmer_kernel(s₁, s₂, 3, 2, w=w) == 3 + 1 + 1
end

@testset "testing Gram matrix" begin
    ℓ = 5
    k = 3
    seqs₁ = string_to_DNA_seq.([random_DNA_seq(10) for i = 1:10])
	seqs₂ = string_to_DNA_seq.([random_DNA_seq(10) for i = 1:12])
    
    # entries to test
    i = 3
    j = 1
    
    # feature vectors
    xᵢ = featurizer(seqs₁[i], ℓ, k)
    xⱼ = featurizer(seqs₂[j], ℓ, k)
    
    # un-normalized
	K = gapped_kmer_kernel_matrix(seqs₁, seqs₂, ℓ, k, normalize=false)
	@test size(K) == (length(seqs₁), length(seqs₂))

	@test K[i, j] ≈ gapped_kmer_kernel(seqs₁[i], seqs₂[j], ℓ, k)
    @test K[i, j] ≈ dot(xᵢ, xⱼ)
	
    # normalized
    K_n = gapped_kmer_kernel_matrix(seqs₁, seqs₂, ℓ, k, normalize=true)
    @test K_n[i, j] ≈ dot(xᵢ, xⱼ) / (norm(xⱼ) * norm(xᵢ))
    
    # symmetric cases
    xᵢ = featurizer(seqs₁[i], ℓ, k)
    xⱼ = featurizer(seqs₁[j], ℓ, k)
    
    K_s = gapped_kmer_kernel_matrix(seqs₁, ℓ, k, normalize=false)
    @test K_s[i, j] ≈ dot(xᵢ, xⱼ)
    
    K_sn = gapped_kmer_kernel_matrix(seqs₁, ℓ, k, normalize=true)
    @test K_sn[i, j] ≈ dot(xᵢ, xⱼ) / (norm(xⱼ) * norm(xᵢ))
end

@testset "testing gapped-kmer featurizer and kernel consistency" begin
    function _gkmer_feature_dot_product_match(s₁::String, s₂::String, ℓ::Int, k::Int)
        s₁′ = string_to_DNA_seq(s₁)
        s₂′ = string_to_DNA_seq(s₂)

        x₁ = featurizer(s₁′, ℓ, k)
        x₂ = featurizer(s₂′, ℓ, k)

        return dot(x₁, x₂) == gapped_kmer_kernel(s₁′, s₂′, ℓ, k)
    end

    @test _gkmer_feature_dot_product_match("AAATCGGCC", "ATTCGGACC", 5, 2)
	@test _gkmer_feature_dot_product_match("AAATCGGCC", "ATTCGGACC", 5, 4)
	@test _gkmer_feature_dot_product_match("AAATGCC", "ATGACC", 3, 2)
    for _ = 1:10
        n₁ = rand(3:10)
        n₂ = rand(3:10)

        ℓ = rand(1:min(n₁, n₂))
        k = rand(1:ℓ)

        s₁ = random_DNA_seq(n₁)
        s₂ = random_DNA_seq(n₂)

        @test _gkmer_feature_dot_product_match(s₁, s₂, ℓ, k)
     end
end

@testset "testing gapped kmer info" begin
	ℓ = 4
	k = 2
	seq_1 = "TTCGAGC"
	seq_2 = "TGGAAGGG"

	seqs = string_to_DNA_seq.([seq_1, seq_2])
    for cut_zeros in [true, false]
        data = gkmer_feature_info(seqs, ℓ, k, cut_zeros)
        @test dot(data[:, seq_1], data[:, seq_2]) == gapped_kmer_kernel(seq_1, seq_2, ℓ, k)
    end
end
