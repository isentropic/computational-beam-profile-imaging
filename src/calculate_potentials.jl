import Pkg
Pkg.activate(joinpath(@__DIR__, "../"))


using LsqFit
using Random, Distributions
using LinearAlgebra


using Random, Distributions
normal_distr = Normal()
uniform_distr = Uniform()

Random.seed!(123)
gaussian_distr = Normal()

const N = 10000
const beam_pipe_radius = 2e-2  # 2cm

const histogram_bins = 50
const no_of_probes = 8  # number of probes
const probe_angles = range(0, length=no_of_probes, step=2*pi/no_of_probes)
const probe_locations = map(θ -> [ beam_pipe_radius * cos(θ), beam_pipe_radius * sin(θ) ], probe_angles)


include("beam_manipulations.jl")
include("minimization.jl")

include("beam_generation.jl")

Current(r, phi) = 1/8e4
J_image(r, phi, a, θ) = Current(r,phi) / (2pi * a) * (a^2 - r^2) / (a^2 + r^2 -2a*r*cos(θ-phi))


function get_current_at_probes(beam, probe_locations=probe_locations)
    noofparticles = size(beam,1)
    rs = zeros(noofparticles)
    phis = zeros(noofparticles)

    @fastmath @inbounds for i in 1:noofparticles
        rs[i] = sqrt(beam[i,1]^2 + beam[i,2]^2)
        phis[i] = atan(beam[i,2], beam[i,1])
    end

    return [get_current_at_probe(beam, rs, phis, noofparticles, probe) for probe in probe_locations]
end


function get_current_at_probe(beam, rs, phis, noofparticles, probe)
    a, θ = norm(probe), atan(probe[2], probe[1])
    answer = 0
    @fastmath @inbounds for i in 1:noofparticles
        answer += J_image(rs[i], phis[i], a, θ)
    end
    return answer
end

function get_current_at_probe_gaussian_beam(beam, probe)
    a, θ = norm(probe), atan(probe[2], probe[1])
    noofparticles = size(beam)[1]
    avgX = mean(beam[:,1])
    avgY = mean(beam[:,2])

    varX = var(beam[:,1])
    varY = var(beam[:,2])

    I_beam = sum([Current(0,0) for i in 1:noofparticles])

    first = 2 * (avgX/a * cos(θ) + avgY/a * sin(θ))
    second = 2 * ((varX^2 - varY^2 + avgX^2 - avgY^2)*cos(2θ) + 2avgX*avgY*sin(2θ))/a^2
    third = 2/a^3 * ( 3*(varX^2 - varY^2) + avgX^2 - avgY^2) *
        (avgX*cos(3θ) + avgY*sin(3θ))

    return I_beam*(1+first+second+third)/(2pi * a)
end



function get_potential_at_probe(probe, beam)
    rel_x = beam[:,1] .- probe[1]
    rel_y = beam[:,2] .- probe[2]
    potentials = @. 1/sqrt(rel_x^2 + rel_y^2)

    return sum(potentials)
end


function get_equal_moment_beams()
    noofbeams = 10
    sigma_xs = range(1e-3, stop=10e-3, length=noofbeams)
    sigma2_difference = 1e-7
    sigma_ys = @. sqrt(sigma_xs^2 - sigma2_difference)

    beams = [getBeam(sigma_xs[i], sigma_ys[i], 0, 1e5) for i in 1:noofbeams]
end


include("plottings.jl")
