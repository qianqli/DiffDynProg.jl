module DiffDynProg

export MaxOperator, Max, LeakyMax, EntropyMax, SquaredMax, smoothmax, smoothmin, smoothmax_argmax, smoothmin_argmin
export maxᵧ, minᵧ, max_argmaxᵧ, min_argminᵧ
export gap_cost_matrix, gumbel_softmax
export DP, getD, getE, getQ
export dynamic_time_warping, ∂DTW
export needleman_wunsch, ∂NW, ∂NW_θ, ∂NW_gs, ∂NW_gt, ∂NW_all

using ChainRulesCore
import ChainRulesCore: rrule

include("maxoperators.jl")
include("dynprog.jl")
include("utils.jl")
include("dynamictimewarping.jl")
include("needlemanwunsch.jl")

end