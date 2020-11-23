using ProgressMeter
using StatsBase

mutable struct MinimizationTarget
    std_x::Float64
    std_y::Float64

    currents::Vector{Float64}
end

mutable struct CandidateSolution
    beam::Array{Float64, 2}
    error::Float64
end
Base.isless(a::CandidateSolution, b::CandidateSolution) = Base.isless(a.error, b.error)

include("minimization_algorithms.jl")


MinimizationTarget(beam::Array{Float64,2}) = MinimizationTarget(mean(beam[:,1].^2), mean(beam[:,2].^2), get_current_at_probes(beam))

mse(a, b) =  mean(@. (a - b)^2)

function mse(a::MinimizationTarget, b::MinimizationTarget)
    currents =  mean(@. (a.currents - b.currents)^2)
    return 100*(a.std_x - b.std_x)^2 + 100*(a.std_y - b.std_y)^2 + currents
end




function get_entropy(candidate::CandidateSolution)
    return std(candidate.beam[:,1]) * std(candidate.beam[:,2])
end


function mutate_beam(guess_beam)
    unit_squeeze = 0.3
    unit_rotation = pi/5
    unit_distance = 0.3 * beam_pipe_radius

    return vcat(
        [shear(guess_beam, randn() * unit_squeeze,i) for i in 1:4],
        [scale(guess_beam, randn() * unit_squeeze, i) for i in 1:4],
        [bend(guess_beam, randn()*0.08, i) for i in 1:4],
        [trapezium(guess_beam, randn()*0.04, i) for i in 1:4]
        )
end



function random_mutating_operation(guess_beam::AbstractArray{T, 2}, mutation_rate = 1.0) where T<:Real
    unit_squeeze = 0.3
    unit_rotation = pi/5
    unit_distance = 0.3 * beam_pipe_radius

    _beam = deepcopy(guess_beam)

    operations = [
        (shear!, unit_squeeze, 4),
        (expand!, unit_squeeze, 2),
        (scale!, unit_squeeze, 2),
        (bend!, 0.05, 4), (trapezium, 0.05, 4),
        (bend!, 0.03, 4)
    ]

    chosen_operations = rand(operations, 5)

    for op in chosen_operations
        # print("Operating")
        op[1](_beam, mutation_rate*(rand()+0.5)*op[2], rand(1:op[3]))
    end

    return _beam
end

function random_mutating_operation(a::CandidateSolution, mutation_rate = 1.0)
    beam = a.beam
    beam2 = random_mutating_operation(beam, mutation_rate)
    return CandidateSolution(beam2, a.error)
end
