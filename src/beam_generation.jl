


function is_within_aperture(beam)
    no_outside_aperture = count(@. sqrt(beam[:,1]^2 + beam[:,2]^2) > beam_pipe_radius)

    if no_outside_aperture > 0
        return false
    else
        return true
    end
end

function get_gauss_beam(a, b, θ=0.0, N=N)
    N = floor(Int, N)
    xs = a .* rand(gaussian_distr, N)
    ys = b .* rand(gaussian_distr, N)

    x_rot = @. xs * cos(θ) - ys * sin(θ)
    y_rot = @. xs * sin(θ) + ys * cos(θ)

    return [x_rot y_rot]
end

function get_random_gaussian_beam()
    max_beam_radius = beam_pipe_radius/3
    min_beam_radius = max_beam_radius/10.0

    radius_sampler = Uniform(min_beam_radius, max_beam_radius)

    x_axis = rand(radius_sampler)
    y_axis = rand(radius_sampler)

    beam = get_gauss_beam(x_axis, y_axis, 0)

    major_axis = maximum((x_axis, y_axis))

    if beam_pipe_radius - 2.8*major_axis > 0
        translation_sampler = Uniform(0, beam_pipe_radius - 3*major_axis)
        translate!(beam, rand(translation_sampler), 1)
        translate!(beam, rand(translation_sampler), 2)
    end

    theta = rand(uniform_distr) * 2 * pi
    rotate!(beam, theta, 1)



    return beam
end


function get_random_distorted_gaussian_beam()
    max_beam_radius = beam_pipe_radius/4
    min_beam_radius = max_beam_radius/10.0

    radius_sampler = Uniform(min_beam_radius, max_beam_radius)

    while true
        x_axis = rand(radius_sampler)
        y_axis = rand(radius_sampler)

        beam = get_gauss_beam(x_axis, y_axis, 0)

        major_axis = maximum((x_axis, y_axis))

        if beam_pipe_radius - 2.5*major_axis > 0
            translation_sampler = Uniform(0, beam_pipe_radius - 3*major_axis)
            translate!(beam, rand(translation_sampler), 1)
            translate!(beam, rand(translation_sampler), 2)
        end



        for i in 1:2
            bend!(beam, rand()*0.2, rand(1:4))
            theta = rand(uniform_distr) * 2 * pi
            rotate!(beam, theta, 1)
            trapezium!(beam, randn()*0.01, rand(1:4))

        end

        if is_within_aperture(beam)
            return beam
        end
    end


end
