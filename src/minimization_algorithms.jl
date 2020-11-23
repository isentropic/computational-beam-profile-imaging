
function adaptive_greedy_fit(true_potentials, error_function=mse)
    unit_squeeze = 0.01
    unit_rotation = pi/16
    unit_distance = 0.01 * beam_pipe_radius

    guess_beam = get_gauss_beam(0.002, 0.002, 0)

    damping_per_failure = 1.5
    @showprogress for i in 1:1000

        translations = [translate(guess_beam, unit_distance, i) for i in 1:4]
        rotations = [rotate(guess_beam, unit_rotation,i) for i in 1:2]
        shrinks = [shrink(guess_beam, unit_squeeze,i) for i in 1:2]
        squeezes = [squeeze(guess_beam, unit_squeeze, i) for i in 1:2]

        trials = vcat(translations, squeezes, rotations)


        potentials = [get_current_at_probes(beam) for beam in trials]
        errors = [error_function(true_potentials, potential) for potential in potentials]

        initial_error = error_function(get_current_at_probes(guess_beam), true_potentials)



        if initial_error < minimum(errors)
            # No better candidate
            unit_squeeze = unit_squeeze / damping_per_failure
            unit_rotation = unit_rotation / damping_per_failure
            unit_distance = unit_distance / damping_per_failure
            println(unit_squeeze)
        else
            least_error, least_error_index= findmin(errors)

            guess_beam = trials[least_error_index]
        end
    end

    return guess_beam
end




function simulated_annealing_fit(true_potentials; error_function=mse, live_plotting=false)
    guess_beam =  get_gauss_beam(0.0003, 0.005, rand() * 2 * pi)

    T_start = 10.0
    T_end = 0.5
    alpha = 0.9

    time_per_temperature = 50


    temperature_range = [T_start]

    while temperature_range[end] > T_end
        current_temp = temperature_range[end]
        push!(temperature_range, current_temp * alpha)
    end


    unit_squeeze = 0.01
    unit_rotation = 2*pi/100
    unit_distance = 0.01 * beam_pipe_radius
    step_size_drop_rate = 1

    unit_squeeze = unit_squeeze * step_size_drop_rate
    unit_rotation = unit_rotation * step_size_drop_rate
    unit_distance = unit_distance * step_size_drop_rate

    @showprogress for T in temperature_range
        if live_plotting
            plt.figure()
            plot(guess_beam)
        end



        transition_made = false
        next_states, next_state_potentials = get_transition_space(guess_beam, unit_squeeze, unit_rotation, unit_distance)
        next_state_errors = [error_function(true_potentials, potential) for potential in next_state_potentials]
        initial_error = error_function(get_current_at_probes(guess_beam), true_potentials)

        for time in 1:time_per_temperature

            if transition_made
                next_states, next_state_potentials = get_transition_space(guess_beam, unit_squeeze, unit_rotation, unit_distance)
                next_state_errors = [error_function(true_potentials, potential) for potential in next_state_potentials]
                initial_error = error_function(get_current_at_probes(guess_beam), true_potentials)
            end


            if minimum(next_state_errors) < initial_error
                # Do greedy descent
                least_error, least_error_index= findmin(next_state_errors)
                guess_beam .= next_states[least_error_index]
                transition_made = true
            else
                random_choice_index = sample(1:length(next_states))
                error_at_next_state = next_state_errors[random_choice_index]
                acceptance_prob = exp( - (error_at_next_state - initial_error) / (T) )

                if acceptance_prob >= rand()
                    transition_made = true
                    guess_beam .= next_states[random_choice_index]
                else
                    transition_made = false
                    continue
                end
            end
        end
    end

    return guess_beam
end

function greedy_fit(true_potentials, guess_beam, error_function=mse)
    maxiters = 1000

    unit_squeeze = 0.01
    unit_rotation = pi/16
    unit_distance = 0.01 * beam_pipe_radius
    initial_error = 0.0
    for j in 1:maxiters

        translations = [translate(guess_beam, unit_distance, i) for i in 1:4]
        rotations = [rotate(guess_beam, unit_rotation,i) for i in 1:2]
        shrinks = [scale(guess_beam, unit_squeeze,i) for i in 1:4]

        trials = vcat(translations, shrinks, rotations)


        potentials = [get_current_at_probes(beam) for beam in trials]
        errors = [error_function(true_potentials, potential) for potential in potentials]

        initial_error = error_function(get_current_at_probes(guess_beam), true_potentials)



        if initial_error < minimum(errors)
            # No better candidate
            break
        else
            least_error, least_error_index= findmin(errors)
            guess_beam = trials[least_error_index]
        end
    end

    return CandidateSolution(guess_beam, initial_error)
end

function greedy_fit(minimization_target::MinimizationTarget, guess_beam, error_function=mse)
    maxiters = 1000

    unit_squeeze = 0.01
    unit_rotation = pi/16
    unit_distance = 0.01 * beam_pipe_radius
    initial_error = 0.0
    for j in 1:maxiters
        trials = vcat(
            [translate(guess_beam, unit_distance, i) for i in 1:4],
            [shear(guess_beam, unit_squeeze, i) for i in 1:4],
            [scale(guess_beam, unit_squeeze, i) for i in 1:4],
            [bend(guess_beam, unit_squeeze, i) for i in 1:4]
            # [trapezium(guess_beam, unit_squeeze, i) for i in 1:4]
        )

        targets = Vector{MinimizationTarget}(undef, length(trials))
        for (i,beam) in enumerate(trials)
            targets[i] = MinimizationTarget(beam)
        end

        errors = [error_function(minimization_target, target) for target in targets]

        initial_error = error_function(minimization_target, MinimizationTarget(guess_beam))



        if initial_error < minimum(errors)
            # No better candidate
            break
        else
            least_error, least_error_index= findmin(errors)
            guess_beam = trials[least_error_index]
        end
    end

    return CandidateSolution(guess_beam, initial_error)
end

function swarm_fit(true_potentials, error_function=mse)
    parallel_instances = 500
    abs_tolerance = 100e-6

    initial_guesses = [get_random_gaussian_beam() for i in 1:parallel_instances]

    candidate_solutions = [CandidateSolution(zeros(10,10), 1.0) for i in 1:parallel_instances];
    Threads.@threads for i in 1:parallel_instances
        # print("$(Threads.threadid())")
        candidate_solutions[i] = greedy_fit(true_potentials, initial_guesses[i])
    end

    sort!(candidate_solutions)

    top_candidates = candidate_solutions[1:10]
    top_entropies = get_entropy.(top_candidates)

    top_entropy_ranks = sortperm(top_entropies)

    return [top_candidates[top_entropy_ranks][1]
        top_candidates[top_entropy_ranks][5]
        top_candidates[top_entropy_ranks][10]]
end



function genetic_mutation_fit(true_potentials, error_function=mse)
    parallel_instances = 30
    abs_tolerance = 100e-6

    initial_guesses = [get_random_gaussian_beam() for i in 1:parallel_instances]
    total_candidate_solutions = CandidateSolution[]
    for i in 1:5
        # println("Iteration: $(i)")

        candidate_solutions = [CandidateSolution(zeros(10,10), 1.0) for i in 1:30];
        Threads.@threads for i in 1:30
            # print("$(Threads.threadid())")
            candidate_solutions[i] = greedy_fit(true_potentials, initial_guesses[i])
        end

        sort!(candidate_solutions)
        [push!(total_candidate_solutions, candidate_solutions[i]) for i in 1:3]
        top_three_beams = [candidate_solution.beam for candidate_solution in candidate_solutions[1:3]]

        initial_guesses = vcat([mutate_beam(beam) for beam in top_three_beams]...)
    end

    sort!(total_candidate_solutions)
    return total_candidate_solutions
end

function genetic_mutation_fit(minimization_target::MinimizationTarget, error_function=mse)
    parallel_instances = 16*10
    abs_tolerance = 100e-6

    initial_guesses = [get_random_distorted_gaussian_beam() for i in 1:parallel_instances]
    total_candidate_solutions = CandidateSolution[]

    @showprogress for i in 1:3
        # println("Iteration: $(i)")

        candidate_solutions = Vector{CandidateSolution}(undef, parallel_instances);
        Threads.@threads for i in 1:parallel_instances
            # print("$(Threads.threadid())")
            candidate_solutions[i] = greedy_fit(minimization_target, initial_guesses[i])
        end

        sort!(candidate_solutions)
        [push!(total_candidate_solutions, candidate_solutions[i]) for i in 1:3]
        top_beams = [candidate_solution.beam for candidate_solution in candidate_solutions[1:10]]

        initial_guesses = vcat([mutate_beam(beam) for beam in top_beams]...)
    end

    sort!(total_candidate_solutions)
    return total_candidate_solutions
end



function stochastic_mutation_descent(minimization_target::MinimizationTarget, error_function=mse)

    parallel_instances = 500
    candidate_solutions = Vector{CandidateSolution}(undef, parallel_instances);
    Threads.@threads for i in 1:parallel_instances
            # print("$(Threads.threadid())")
            candidate_solutions[i] = greedy_fit(minimization_target, get_random_distorted_gaussian_beam())
    end
    sort!(candidate_solutions)
    init_guess = candidate_solutions[1]


    iterations = 100
    parallel_mutations = 20

    mutations = 0
    @showprogress for k in 1:iterations
        mutated_beams = Vector{CandidateSolution}(undef, parallel_mutations);

        Threads.@threads for i in 1:parallel_mutations
            mutated_beams[i] = greedy_fit(minimization_target, random_mutating_operation(init_guess, .985^i).beam)
        end

        sort!(mutated_beams)

        mutated = mutated_beams[1]

        if mutated.error < init_guess.error
            mutations += 1
            init_guess = mutated
        end
#         mutated = random_mutating_operation(init_guess, .998^i)
#         mutated = greedy_fit(minimization_target, mutated.beam)




    end

    print("Successful mutations: $mutations")

    return init_guess
end
