using PyCall
using PyPlot;

import Plots
using PlotThemes
import StatsPlots

# Plots.theme(:seaborn_deep)

const plt = PyPlot
const P = Plots

sns = pyimport("seaborn")

# rc("font",**{"family":"sans-serif","sans-serif":["Helvetica"]})
## for Palatino and other serif fonts use:
# rc("font",**{"family":"serif","serif":["Computer Modern"]})
sns.set()
sns.set_style("whitegrid")

rcParams = PyDict(matplotlib."rcParams")
# rcParams["font.family"] = "serif"
# rcParams["font.serif"] = "Computer Modern Roman"
# rcParams["font.sans-serif"] = "Palantino"
# rcParams["text.usetex"] = true

rcParams["mathtext.fontset"] = "custom"
rcParams["mathtext.rm"] = "Bitstream Vera Sans"
rcParams["mathtext.it"] = "Bitstream Vera Sans:italic"
rcParams["mathtext.bf"] = "Bitstream Vera Sans:bold"
rcParams["mathtext.fontset"] = "stix"
rcParams["font.family"] = "STIXGeneral"
rcParams["axes.formatter.use_mathtext"] = true
rcParams["errorbar.capsize"] = 4

PyPlot.svg(false);


function PyPlot.plot(x::AbstractArray{T, 2}; kwargs...) where T<:Real
    fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(4, 3))

    edges = range(-beam_pipe_radius, stop=beam_pipe_radius, length=histogram_bins+1)

    Z, xedges, yedges = hist2D(x[:,2], x[:,1], bins=[edges, edges]; kwargs...)

    xlabel("x [m]")
    ylabel("y [m]")
    plt.pcolormesh(edges, edges, Z)
    plt.colorbar()
    tight_layout()
end

function PyPlot.plot(x::AbstractArray{T, 2}, y::AbstractArray{T, 2}; kwargs...) where T<:Real
    fig, axes = plt.subplots(nrows=1, ncols=2, sharey=true,figsize=(7.2, 3))

    edges = range(-beam_pipe_radius, stop=beam_pipe_radius, length=histogram_bins+1)

    sca(axes[1])
    Z, xedges, yedges = hist2D(x[:,2], x[:,1], bins=[edges, edges]; kwargs...)
    plt.pcolormesh(edges, edges, Z)
    plt.colorbar()
    xlabel("x [m]")
    ylabel("y [m]")
    title("Original")

    sca(axes[2])
    Z, xedges, yedges = hist2D(y[:,2], y[:,1], bins=[edges, edges]; kwargs...)
    plt.pcolormesh(edges, edges, Z)
    plt.colorbar()
    title("Reproduction")
    xlabel("x [m]")
    tight_layout()
end

function compare_plot(original_beam, fitted_results::Vector{CandidateSolution}; kwargs...)
    noofplots = length(fitted_results)+1
    fig, axes = plt.subplots(nrows=1, ncols=noofplots, sharey=true,figsize=(noofplots*3.15, 3))

    edges = range(-beam_pipe_radius, stop=beam_pipe_radius, length=histogram_bins+1)

    sca(axes[1])
    Z, xedges, yedges = hist2D(original_beam[:,2], original_beam[:,1], bins=[edges, edges]; kwargs...)
    plt.pcolormesh(edges, edges, Z)
    plt.colorbar()
    xlabel("x [m]")
    ylabel("y [m]")
    title("Original")
    for i in 2:noofplots
        sca(axes[i])
        y = fitted_results[i-1].beam
        error = fitted_results[i-1].error

        Z, xedges, yedges = hist2D(y[:,2], y[:,1], bins=[edges, edges]; kwargs...)
        plt.pcolormesh(edges, edges, Z)
        plt.colorbar()
        title("Reproduction $(i-1), error: $(round(error, sigdigits=5))")
        xlabel("x [m]")

    end
    tight_layout()
end
