### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 6e289366-7790-11ef-2381-2d4900e853e4
begin
	import Pkg; Pkg.activate()
	using Revise
	
	push!(LOAD_PATH, joinpath(pwd(), "src"))
	using GappedKmers, BioSequences, Test, DataFrames, BenchmarkTools, LinearAlgebra, CairoMakie
end

# ╔═╡ c9d45f38-cdb8-4d40-a5e7-6d9cabe33119
md"# `GappedKmers.jl`"

# ╔═╡ 0caebb26-c638-4379-a139-e5c440c0d700
md"## list the gapped k-mers"

# ╔═╡ 36ed8d7b-6723-4869-a336-5a15d6ea4078
@doc list_of_gapped_kmers

# ╔═╡ 7c3a61bb-234f-49dc-ad11-685c8d8647ac
ℓ = 3

# ╔═╡ 1a970302-cab9-4a27-9f57-9e97d9647768
k = 2

# ╔═╡ 4e5b2c7c-6b43-4708-b733-b52ea35312a3
list_of_gapped_kmers(ℓ, k)

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
@doc featurizer

# ╔═╡ a88fa7db-d8ea-46d2-9637-2b9dff33f36e
x₁ = featurizer(seq₁, ℓ, k)

# ╔═╡ 3efd8b19-ab3d-4da5-a0b3-3781a82d6594
x₂ = featurizer(seq₂, ℓ, k)

# ╔═╡ 398f2a54-1083-4d74-a99c-2df7c1495dcb
data = gkmer_feature_info([seq₁, seq₂], ℓ, k, false)

# ╔═╡ 76f3133f-ce3c-4621-b88c-1b62786953de
md"## gapped k-mer kernel"

# ╔═╡ 5a302c7b-c71b-40dc-9f2c-f718b0095758
@doc gapped_kmer_kernel

# ╔═╡ 2165106a-2f4a-42fb-998b-f349857798e7
gapped_kmer_kernel(seq₁, seq₂, ℓ, k)

# ╔═╡ 5d053946-5666-41aa-ae73-3c69090d2173
dot(x₁, x₂)

# ╔═╡ f605fcff-cfb1-4946-832e-9230c1ab6c28
seq_1 = string_to_DNA_seq(random_DNA_seq(10))

# ╔═╡ be742de2-cf23-4947-8d32-497fa0f391b6
md"## Gram (kernel) matrix"

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

# ╔═╡ Cell order:
# ╠═6e289366-7790-11ef-2381-2d4900e853e4
# ╟─c9d45f38-cdb8-4d40-a5e7-6d9cabe33119
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
# ╠═f605fcff-cfb1-4946-832e-9230c1ab6c28
# ╟─be742de2-cf23-4947-8d32-497fa0f391b6
# ╠═d0eb6410-3a30-448e-b8cf-ae54e7047844
# ╠═a9e67ae5-548a-421d-b339-dfbba6e09815
# ╠═06e56447-24f5-469c-bdc7-2fd1e62d17c7
# ╠═8d3252bf-5c72-41f6-93c6-474e9662c468
# ╠═e15cc2cc-8f40-4d2b-a4ee-a13d1c03307f
# ╠═81cf013e-6b9d-4ec7-b185-d652c0df90b4
# ╠═a54b25dc-928e-4eb6-ad7b-ae2d5c62d040
