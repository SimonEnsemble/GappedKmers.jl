### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 6e289366-7790-11ef-2381-2d4900e853e4
begin
	import Pkg; Pkg.activate()

	using Revise
	
	push!(LOAD_PATH, joinpath(pwd(), "src"))
	using GappedKmers, BioSequences, Test, DataFrames, 
		  BenchmarkTools, LinearAlgebra, CairoMakie, PlutoUI
end

# ╔═╡ e8c54f94-1d7a-4ac7-a31e-9d23c7654270
TableOfContents()

# ╔═╡ c9d45f38-cdb8-4d40-a5e7-6d9cabe33119
md"# 🧬 `GappedKmers.jl`

our Julia package `GappedKmers.jl` is intended for featurizing DNA sequences using gapped k-mers, for machine learning tasks on DNA sequences. this Pluto notebook illustrates the capabilities of `GappedKmers.jl` and serves as documentation.

!!! reference
	to learn about gapped k-mers for describing DNA sequences, see:
	> M. Ghandi, D. Lee, Mohammad-Noori, M. Beer. Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features. _PLoS Computational Biology_. (2014) [link.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003711)
"

# ╔═╡ 4cfd08b9-4a99-4914-b29e-5bc57bb73653
md"!!! note
    we rely on [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl) to represent our DNA sequences and search for k-mers with mis-matches allowed."

# ╔═╡ 0caebb26-c638-4379-a139-e5c440c0d700
md"## list the gapped k-mers"

# ╔═╡ 36ed8d7b-6723-4869-a336-5a15d6ea4078
@doc list_of_gapped_kmers

# ╔═╡ 7c3a61bb-234f-49dc-ad11-685c8d8647ac
ℓ = 3 # length of subsequence

# ╔═╡ 1a970302-cab9-4a27-9f57-9e97d9647768
k = 2 # number of informative (non-wildcard positions)

# ╔═╡ 4e5b2c7c-6b43-4708-b733-b52ea35312a3
gkmers = list_of_gapped_kmers(ℓ, k)

# ╔═╡ 4ad81401-84d7-4043-bb47-9066c14ee1e2
@doc number_of_gapped_kmers

# ╔═╡ 1704623d-e0f2-4d44-bf91-20b4d6b5969f
number_of_gapped_kmers(ℓ, k)

# ╔═╡ 6e89231f-cd6c-4e1a-b421-a6640b2cd43d
md"## gapped k-mer feature vector"

# ╔═╡ f9224e80-0f99-4731-9801-3e5308c4e789
@doc string_to_DNA_seq

# ╔═╡ 573bca78-cdf5-4ea4-8db2-561bf9dbb7cf
seq₁ = string_to_DNA_seq("ATTGGT")

# ╔═╡ eedb257f-5240-44e4-b8b7-c72f58a4676a
seq₂ = string_to_DNA_seq("GTGGT")

# ╔═╡ 3656998e-f76c-475b-9d01-86484af2f488
@doc gkmer_feature

# ╔═╡ a88fa7db-d8ea-46d2-9637-2b9dff33f36e
x₁ = gkmer_feature(seq₁, ℓ, k)

# ╔═╡ 3efd8b19-ab3d-4da5-a0b3-3781a82d6594
x₂ = gkmer_feature(seq₂, ℓ, k)

# ╔═╡ 398f2a54-1083-4d74-a99c-2df7c1495dcb
data = gkmer_feature_info([seq₁, seq₂], ℓ, k, false)

# ╔═╡ 76f3133f-ce3c-4621-b88c-1b62786953de
md"## gapped k-mer kernel"

# ╔═╡ 5a302c7b-c71b-40dc-9f2c-f718b0095758
@doc gapped_kmer_kernel

# ╔═╡ 2165106a-2f4a-42fb-998b-f349857798e7
gapped_kmer_kernel(seq₁, seq₂, ℓ, k)

# ╔═╡ 5d053946-5666-41aa-ae73-3c69090d2173
dot(x₁, x₂) # gives same as dot product between features

# ╔═╡ be742de2-cf23-4947-8d32-497fa0f391b6
md"## kernel matrix"

# ╔═╡ d0eb6410-3a30-448e-b8cf-ae54e7047844
seqs₁ = string_to_DNA_seq.([random_DNA_seq(10) for i = 1:10])

# ╔═╡ a9e67ae5-548a-421d-b339-dfbba6e09815
seqs₂ = string_to_DNA_seq.([random_DNA_seq(10) for i = 1:12])

# ╔═╡ 06e56447-24f5-469c-bdc7-2fd1e62d17c7
@doc gapped_kmer_kernel_matrix

# ╔═╡ 8d3252bf-5c72-41f6-93c6-474e9662c468
K = gapped_kmer_kernel_matrix(seqs₁, seqs₂, ℓ, k, normalize=false)

# ╔═╡ e15cc2cc-8f40-4d2b-a4ee-a13d1c03307f
begin
	i = 5
	j = 2
	K[i, j]
end

# ╔═╡ 81cf013e-6b9d-4ec7-b185-d652c0df90b4
gapped_kmer_kernel(seqs₁[i], seqs₂[j], ℓ, k)

# ╔═╡ a54b25dc-928e-4eb6-ad7b-ae2d5c62d040
begin
	fig, ax, hm = heatmap(K; colormap=:turku10,  
		axis=(; aspect=DataAspect(), xlabel="DNA seq.", ylabel="DNA seq.",
		xticks=1:length(seqs₁), yticks=1:length(seqs₂))
	)
	Colorbar(fig[:, end+1], hm, label="similarity score")
	fig
end

# ╔═╡ a083e6db-dcf2-42dc-90cf-23654f5be505
md"### Gram matrix"

# ╔═╡ 97c40301-4311-4bd0-b029-d1caf4db3701
G = gapped_kmer_kernel_matrix(seqs₁, ℓ, k, normalize=false)

# ╔═╡ 742409e0-6735-4316-abe6-a1f3ccc0f19e
X = gkmer_feature_matrix(seqs₁, ℓ, k)

# ╔═╡ 5f7bbe4f-05fd-42aa-9bff-0aed310e3972
all(X' * X .== G)

# ╔═╡ 658d68b1-eeba-4b63-842c-b76cf10c08e2
md"### benchmarking"

# ╔═╡ 2e3f544d-2578-4f0d-9002-75395fb92f8e
md"using `BenchmarkTools.jl` to measure the elapsed time to compute the Gram matrix between all pairs of 1000 length-10 DNA sequences"

# ╔═╡ e6de911e-f596-4a41-9dd7-d526d65850bd
begin
	BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60
	benchmarking_seqs = string_to_DNA_seq.([random_DNA_seq(10) for i = 1:1000])
end

# ╔═╡ 148ab4c2-ae8b-4772-b7eb-eaae8b5fe3ad
@benchmark gapped_kmer_kernel_matrix(
		benchmarking_seqs, benchmarking_seqs, ℓ, k
	)

# ╔═╡ Cell order:
# ╠═6e289366-7790-11ef-2381-2d4900e853e4
# ╠═e8c54f94-1d7a-4ac7-a31e-9d23c7654270
# ╟─c9d45f38-cdb8-4d40-a5e7-6d9cabe33119
# ╟─4cfd08b9-4a99-4914-b29e-5bc57bb73653
# ╟─0caebb26-c638-4379-a139-e5c440c0d700
# ╠═36ed8d7b-6723-4869-a336-5a15d6ea4078
# ╠═7c3a61bb-234f-49dc-ad11-685c8d8647ac
# ╠═1a970302-cab9-4a27-9f57-9e97d9647768
# ╠═4e5b2c7c-6b43-4708-b733-b52ea35312a3
# ╠═4ad81401-84d7-4043-bb47-9066c14ee1e2
# ╠═1704623d-e0f2-4d44-bf91-20b4d6b5969f
# ╟─6e89231f-cd6c-4e1a-b421-a6640b2cd43d
# ╠═f9224e80-0f99-4731-9801-3e5308c4e789
# ╠═573bca78-cdf5-4ea4-8db2-561bf9dbb7cf
# ╠═eedb257f-5240-44e4-b8b7-c72f58a4676a
# ╠═3656998e-f76c-475b-9d01-86484af2f488
# ╠═a88fa7db-d8ea-46d2-9637-2b9dff33f36e
# ╠═3efd8b19-ab3d-4da5-a0b3-3781a82d6594
# ╠═398f2a54-1083-4d74-a99c-2df7c1495dcb
# ╟─76f3133f-ce3c-4621-b88c-1b62786953de
# ╠═5a302c7b-c71b-40dc-9f2c-f718b0095758
# ╠═2165106a-2f4a-42fb-998b-f349857798e7
# ╠═5d053946-5666-41aa-ae73-3c69090d2173
# ╟─be742de2-cf23-4947-8d32-497fa0f391b6
# ╠═d0eb6410-3a30-448e-b8cf-ae54e7047844
# ╠═a9e67ae5-548a-421d-b339-dfbba6e09815
# ╠═06e56447-24f5-469c-bdc7-2fd1e62d17c7
# ╠═8d3252bf-5c72-41f6-93c6-474e9662c468
# ╠═e15cc2cc-8f40-4d2b-a4ee-a13d1c03307f
# ╠═81cf013e-6b9d-4ec7-b185-d652c0df90b4
# ╠═a54b25dc-928e-4eb6-ad7b-ae2d5c62d040
# ╟─a083e6db-dcf2-42dc-90cf-23654f5be505
# ╠═97c40301-4311-4bd0-b029-d1caf4db3701
# ╠═742409e0-6735-4316-abe6-a1f3ccc0f19e
# ╠═5f7bbe4f-05fd-42aa-9bff-0aed310e3972
# ╟─658d68b1-eeba-4b63-842c-b76cf10c08e2
# ╟─2e3f544d-2578-4f0d-9002-75395fb92f8e
# ╠═e6de911e-f596-4a41-9dd7-d526d65850bd
# ╠═148ab4c2-ae8b-4772-b7eb-eaae8b5fe3ad
