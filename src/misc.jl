"""
    random_DNA_seq(n)

sample a random DNA sequence of length `n`. returns a `String`.
"""
function random_DNA_seq(n::Int)
    nucleotides = ["A", "T", "C", "G"]
    seq = ""
    for i = 1:n
        seq *= rand(nucleotides)
    end
    return seq
end
