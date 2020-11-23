function translate!(beam2, unit_distance, direction)
    # 1,2,3,4 = right, top, left, bot
    N = size(beam2)[1]
    @fastmath @inbounds for i in 1:N
        if direction == 1
            beam2[i,1] = beam2[i,1] + unit_distance
        elseif direction == 2
            beam2[i,2] = beam2[i,2] + unit_distance
        elseif direction == 3
            beam2[i,1] = beam2[i,1] - unit_distance
        elseif direction == 4
            beam2[i,2] = beam2[i,2] - unit_distance
        else
            throw(ArgumentError)
        end
    end
    return nothing
end


function scale!(beam2, unit_squeeze, direction)
    # 1,2 = CCW, CW
    N = size(beam2)[1]
    @fastmath @inbounds for i in 1:N
        if direction == 1
            beam2[i,1] = beam2[i,1] * (1+unit_squeeze)
        elseif direction == 2
            beam2[i,1] = beam2[i,1] / (1+unit_squeeze)
        elseif direction == 3
            beam2[i,2] = beam2[i,2] * (1+unit_squeeze)
        elseif direction == 4
            beam2[i,2] = beam2[i,2] / (1+unit_squeeze)
        else
            throw(ArgumentError)
        end
    end

    return nothing
end


function expand!(beam2, unit_squeeze, direction)
    # 1,2 = CCW, CW
    if direction == 1
        scale!(beam2, unit_squeeze, 1)
        scale!(beam2, unit_squeeze, 3)
    elseif direction == 2
        scale!(beam2, unit_squeeze, 2)
        scale!(beam2, unit_squeeze, 4)
    else
        throw(ArgumentError)
    end
end

function expand(beam, unit_shear, direction)
    beam2 = deepcopy(beam)
    expand!(beam2, unit_shear, direction)
    return beam2
end





function rotate!(beam2, unit_rotation, direction)
    # 1,2 = CCW, CW
    if direction == 1
    elseif direction == 2
        unit_rotation = -unit_rotation
    else
        throw(ArgumentError)
    end

    xcenter = mean(beam2[:,1] )
    ycenter = mean(beam2[:,2] )

    beam2 .= beam2 * [cos(unit_rotation) sin(unit_rotation); -sin(unit_rotation) cos(unit_rotation)]
    return nothing
end


function shear!(beam2, unit_shear, direction)
    # 1,2 = CCW, CW
    x = @view beam2[:,1]
    y = @view beam2[:,2]

    N = size(beam2)[1]

    xcenter = mean(x)
    ycenter = mean(y)
    @fastmath @inbounds for i in 1:N
        if direction == 1
            x[i] = x[i] + unit_shear * y[i]
        elseif direction == 2
            y[i] = y[i] + unit_shear * x[i]
        elseif direction == 3
            x[i] = x[i] - unit_shear * y[i]
        elseif direction == 4
            y[i] = y[i] - unit_shear * x[i]
        else
            throw(ArgumentError)
        end
    end
    return nothing
end

function shear(beam, unit_shear, direction)
    beam2 = deepcopy(beam)
    shear!(beam2, unit_shear, direction)
    return beam2
end


function rotate(beam, unit_rotation, direction)
    beam2 = deepcopy(beam)
    rotate!(beam2, unit_rotation, direction)
    return beam2
end


function scale(beam, unit_squeeze, direction)
    beam2 = deepcopy(beam)
    scale!(beam2, unit_squeeze, direction)
    return beam2
end


function skew_scale(beam, unit_squeeze, direction)
    beam2 = deepcopy(beam)

    rotate!(beam2, pi/4, 1)
    scale!(beam2, unit_squeeze, direction)
    rotate!(beam2, pi/4, 2)

    return beam2
end


function translate(beam, unit_distance, direction)
    beam2 = deepcopy(beam)
    translate!(beam2, unit_distance, direction)
    return beam2
end

function get_transition_space(guess_beam, unit_squeeze, unit_rotation, unit_distance)

    translations = [translate(guess_beam, unit_distance, i) for i in 1:4]
    rotations = [rotate(guess_beam, unit_rotation,i) for i in 1:2]
    shrinks = [scale(guess_beam, unit_squeeze,i) for i in 1:4]

    next_states = vcat(translations, squeezes, rotations, shrinks)


    next_state_potentials = [get_current_at_probes(beam) for beam in next_states]

    return next_states, next_state_potentials
end



include("non_linear_beam_manipulations.jl")
