using LsqFit
using Random, Distributions
using BlackBoxOptim
Random.seed!(123)
gaussian_distr = Normal()

width = 3e-3
N = 1e4
radius = 2e-2  # 2cm
no_of_probes = 8  # number of probes
probe_angles = range(0, length=no_of_probes, step=2*pi/no_of_probes)
probe_locations = map(θ -> [ radius * cos(θ), radius * sin(θ) ], probe_angles)


"""
a = major semiaxis
b = minor semiaxis

θ - active rotation angle

...
cos θ -sin θ

sin θ  cos θ
...
"""
function getBeam(a, b, θ, N=N)
    N = floor(Int, N)
    xs = a .* rand(gaussian_distr, N)
    ys = b .* rand(gaussian_distr, N)

    x_rot = @. xs * cos(θ) - ys * sin(θ)
    y_rot = @. xs * sin(θ) + ys * cos(θ)

    return [x_rot y_rot]
end


beam = getBeam(10e-3, 1e-3, pi/4)

function getPotentialAtAProbe(probe, beam)
    rel_x = beam[:,1] .- probe[1]
    rel_y = beam[:,2] .- probe[2]
    potentials = @. 1/sqrt(rel_x^2 + rel_y^2)

    return sum(potentials)
end

potentials = [getPotentialAtAProbe(probe, beam) for probe in probe_locations]

function getAllMultipolesModel(order, angles, p)
    dc = p[1]
    y = zeros(length(angles)) .+ dc
    for i in 1:(order)
        y += @. p[2i] * sin((i)*angles + p[2i+1])
    end
    return y
end

function fitAllMultipolesBBO(potentials, maxMultipole=4, angles=probe_angles; kwargs...)
    function getResidialSumMultipolesModel(p)
        y = potentials
        x = angles
        y_model = getAllMultipolesModel(maxMultipole, angles, p)
        return sum(@. (y - y_model)^2)
    end

    lb = zeros(1 + 2*maxMultipole)
    ub = ones(1 + 2*maxMultipole)
    ub .= ub * 2 *pi

    ub[1] = maximum(potentials)
    for i in 2:2:(1+maxMultipole*2)
        ub[i] = maximum(potentials) - minimum(potentials)
    end
    tupleBounds = [(lb[i], ub[i]) for i in 1:length(lb)]
    return bboptimize(getResidialSumMultipolesModel; SearchRange = tupleBounds, kwargs...)
end


function fitAllMultipoles(potentials, maxMultipole=4, angles=probe_angles; kwargs...)
    getAllMultipolesModelLocal(angles, p) = getAllMultipolesModel(maxMultipole, angles, p)
    p0 = zeros(1 + 2*maxMultipole)
    lb = zeros(1 + 2*maxMultipole)
    ub = ones(1 + 2*maxMultipole)
    ub .= ub * 2 *pi

    ub[1] = maximum(potentials)
    for i in 2:2:(1+maxMultipole*2)
        ub[i] = maximum(potentials) - minimum(potentials)
    end


    curve_fit(getAllMultipolesModelLocal, angles, potentials, p0, lower=lb, upper=ub; kwargs...)
end



using PyPlot
