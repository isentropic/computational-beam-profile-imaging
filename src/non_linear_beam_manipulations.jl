
function trapezium!(beam2, unit, direction)
    # 1,2 = CCW, CW
    x = @view beam2[:,1]
    y = @view beam2[:,2]

    N = size(beam2)[1]

    xcenter = mean(x)
    ycenter = mean(y)
    varx = std(x)
    vary = std(y)
    @fastmath @inbounds for i in 1:N
        if direction == 1
            y[i] = y[i] * (1.0 + unit)
            x[i] = x[i] * (1.0 + 5.0 * unit * (y[i] - ycenter)/vary )
        elseif direction == 2
            x[i] = x[i] * (1.0 + unit)
            y[i] = y[i] * (1.0 + 5.0 * unit * (x[i] - xcenter)/varx )
        elseif direction == 3
            y[i] = y[i] * (1.0 + unit)
            x[i] = x[i] * (1.0 - 5.0 * unit * (y[i] + ycenter)/vary )
        elseif direction == 4
            x[i] = x[i] * (1.0 + unit)
            y[i] = y[i] * (1.0 - 5.0 * unit * (x[i] + xcenter)/varx )
        else
            throw(ArgumentError)
        end
    end
    return nothing
end

function trapezium(beam, unit, direction)
    beam2 = deepcopy(beam)
    trapezium!(beam2, unit, direction)
    return beam2
end


function nonlinear_scale!(beam2, unit_squeeze, direction)
    # 1,2 = CCW, CW
    N = size(beam2)[1]

    x = @view beam2[:,1]
    y = @view beam2[:,2]

    xcenter = mean(x)
    ycenter = mean(y)

    @fastmath @inbounds for i in 1:N
        if direction == 1
            x[i] = x[i] + x[i] * (unit_squeeze) * (y[i] > ycenter)
        elseif direction == 2
            x[i] = x[i] + x[i] / (unit_squeeze) * (y[i] > ycenter)
        elseif direction == 3
            y[i] = y[i] + y[i] * (unit_squeeze) * (x[i] > xcenter)
        elseif direction == 4
            y[i] = y[i] + y[i] / (unit_squeeze) * (x[i] > xcenter)
        else
            throw(ArgumentError)
        end
    end

    return nothing
end

function nonlinear_scale(beam, unit_squeeze, direction)
    beam2 = deepcopy(beam)
    nonlinear_scale!(beam2, unit_squeeze, direction)
    return beam2
end


function nonlinear_shear!(beam2, unit_shear, direction)
    # 1,2 = CCW, CW
    x = @view beam2[:,1]
    y = @view beam2[:,2]

    N = size(beam2)[1]

    xcenter = mean(x)
    ycenter = mean(y)
    @fastmath @inbounds for i in 1:N
        if direction == 1
            x[i] = x[i] + unit_shear * (y[i]-ycenter) * (y[i] > ycenter)
        elseif direction == 2
            y[i] = y[i] + unit_shear * (x[i]-xcenter) * (x[i] > xcenter)
        elseif direction == 3
            x[i] = x[i] - unit_shear * (y[i]-ycenter) * (y[i] > ycenter)
        elseif direction == 4
            y[i] = y[i] - unit_shear * (x[i]-xcenter) * (x[i] > xcenter)
        else
            throw(ArgumentError)
        end
    end
    return nothing
end


function nonlinear_shear(beam, unit_shear, direction)
    beam2 = deepcopy(beam)
    nonlinear_shear!(beam2, unit_shear, direction)
    return beam2
end


function bend!(beam2, unit_shear, direction)
    # 1,2 = CCW, CW
    x = @view beam2[:,1]
    y = @view beam2[:,2]

    N = size(beam2)[1]

    xcenter = mean(x)
    ycenter = mean(y)
    @fastmath @inbounds for i in 1:N
        if direction == 1
            x[i] = x[i] + unit_shear * abs(y[i]-ycenter)
        elseif direction == 2
            y[i] = y[i] + unit_shear * abs(x[i]-xcenter)
        elseif direction == 3
            x[i] = x[i] - unit_shear * abs(y[i]-ycenter)
        elseif direction == 4
            y[i] = y[i] - unit_shear * abs(x[i]-xcenter)
        else
            throw(ArgumentError)
        end
    end
    return nothing
end

function bend(beam, unit_shear, direction)
    beam2 = deepcopy(beam)
    bend!(beam2, unit_shear, direction)
    return beam2
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
            translation_sampler = Uniform(0, beam_pipe_radius - 2.8*major_axis)
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
