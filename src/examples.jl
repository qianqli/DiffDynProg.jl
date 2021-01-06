#=
EXAMPLES FOR DIFFERENTIABLE DTW

Created on 06/01/2021

@author: dimiboeckaerts
=#

# 0 - LIBRARIES
# --------------------------------------------------
include("./DiffDynProg.jl")
using DiffDynProg


# 1 - TIME SERIES EXAMPLE
# --------------------------------------------------
X = [3.0, 1.0, 0.0, 0.0, 1.0, 3.0, 5.0, 6.0, 6.0, 5.0, 3.0, 1.0, 0.0, 0.0, 1.0, 3.0]
Y = [0.0, 1.0, 3.0, 5.0, 6.0, 6.0, 5.0, 3.0]

theta = (X' .- Y).^2
mo = SquaredMax(5.0)

dtw = DTW(theta)
D, E = ∂DPW!(mo::MaxOperator, theta, dtw::DTW)
scheme = :Greens_4
plot(heatmap(theta, color=scheme, yflip=true), heatmap(D, color=scheme, yflip=true),
    heatmap(E, color=scheme, yflip=true))

# plot the Mensch fig
lay = @layout [a b]
plot(scatter(X), heatmap(E, color=scheme, yflip=true), seriestype=[:scatter, :heatmap])


# 2 - SIMPLE SEQUENCES
# --------------------------------------------------
"""
If we want to apply DTW to biological sequences, we're actually doing global
sequence alignment without gap costs. We do need a distance matrix as input,
and this can be based on the substitution matrix (but the inverse of it, in
sequence alignment we try to maximize a score and in DTW we try to minimize a
distance).

To do:
- implement function to compute distance matrix (with gap costs)
- sample alignments as a random walk
"""
s1 = "KLMSP"
s2 = "KPMMSQ"

dist_matrix = [0.0 2.0 3.0 4.0 5.0;
                2.0 1.0 2.0 3.0 3.0;
                3.0 2.0 0.0 2.0 3.0;
                4.0 3.0 1.0 1.0 2.0;
                5.0 4.0 3.0 1.0 1.0;
                6.0 5.0 4.0 3.0 2.0]

dtw = DTW(dist_matrix)
D, E = ∂DPW!(mo::MaxOperator, dist_matrix, dtw::DTW)
plot(heatmap(dist_matrix, color=scheme, yflip=true), heatmap(D, color=scheme, yflip=true),
    heatmap(E[1:end-1, 1:end-1], color=scheme, yflip=true))
