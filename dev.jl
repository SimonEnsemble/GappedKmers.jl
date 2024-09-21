### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 6e289366-7790-11ef-2381-2d4900e853e4
begin
	import Pkg; Pkg.activate()
	using Revise
	push!(LOAD_PATH, joinpath(pwd(), "src"))
	using GappedKmers, BioSequences
end

# ╔═╡ 4e5b2c7c-6b43-4708-b733-b52ea35312a3
all_gkmers = list_of_gapped_kmers(3, 2)

# ╔═╡ 573bca78-cdf5-4ea4-8db2-561bf9dbb7cf
convert_to_regex.(all_gkmers)

# ╔═╡ cf52c674-08f9-4a2b-a7ad-3985aaf93eab
string_to_DNA_seq("TTAGGC")

# ╔═╡ 15f516bd-5a37-4961-b3a9-c3eacb1ba607
x = "ATTC"

# ╔═╡ Cell order:
# ╠═6e289366-7790-11ef-2381-2d4900e853e4
# ╠═4e5b2c7c-6b43-4708-b733-b52ea35312a3
# ╠═573bca78-cdf5-4ea4-8db2-561bf9dbb7cf
# ╠═cf52c674-08f9-4a2b-a7ad-3985aaf93eab
# ╠═15f516bd-5a37-4961-b3a9-c3eacb1ba607
