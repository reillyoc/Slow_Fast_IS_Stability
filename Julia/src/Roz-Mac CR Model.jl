using Distributed
using PyPlot
using Parameters
using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Statistics
using NLsolve
using Base.Threads
using DataFrames
using CSV


# # Model Parameters
@with_kw mutable struct CRPar
    r = 1.5
    k = 2.0
    a = 0.8
    e = 0.8
    b = 0.8
    m = 0.3

    noise = 0.00
end

##Roz-Mac CR logistic Resource
function CR_mod!(du, u, p, t)
    @unpack r, k, a, b, e, m = p
    R, C, = u
    
    du[1] = r * R * (1 - R / k) - (a * R * C)/(R + b)
    du[2] = (e * a * R * C)/(R + b) - m * C

    return
end

##Stochastic CR Roz-Mac CR logistic Resource
function CR_mod(u, p)
    du = similar(u)
    CR_mod!(du, u, p, 0.0)
    return du
end

#Add noise to Stochastic
function stoch_CR_Mod!(du, u, p2, t)
    @unpack  noise = p2
    R, C = u

    du[1] = noise * R
    du[2] = noise * C
    return du 
end

##Isocline plotting function
function retrieve_RC(vectordata, RC)
    newarray = zeros(size(vectordata)[1],size(vectordata)[2])

    for i in 1:size(vectordata)[1]
        for j in 1:size(vectordata)[2]
            newarray[i,j] = vectordata[i,j][RC]
        end
    end
    return newarray
end

##Plot Resource isocline
function res_isocline(R, p)
    @unpack r, k, a, b = p
    return (r/a) * (1 - R/k) * (R+b)/a
end

##Plot Consumer isocline
function cons_isocline(C, p)
    @unpack m, a, e, b = p
    return (m * b) / (e * a - m)
end

##Plot different isoclines based on amount of dampening (D)
let
    par = CRPar(a = 0.6)

    Rrange = 0.0:0.1:6.0
    Crange = 0.0:0.1:6.0
    resconrange = 0.0:0.005:6.0

    CR_mod_iso_plot = figure()
    plot(collect(Rrange), [res_isocline(R, par) for R in Rrange], color="navy")
    plot([cons_isocline(C, par) for C in Crange], collect(Crange),color="#FF9191")

    xlabel("Resource")
    ylabel("Consumer")
    xlim(0,6)
    ylim(0,6)

    #savefig("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Julia/Figures/Isoclines.pdf", dpi = 300)
    return CR_mod_iso_plot
end

#extinction threshold
threshold = 1e-12

#Condition to trigger when any state variable falls below the threshold
function condition(u, t, integrator)
    #Trigger if any population is positive but below the threshold, or if it is negative
    for i in 1:length(u)
        if u[i] < threshold
            return true  #Trigger the callback
        end
    end
    return false
end

#Effect function to set the state variable to zero or threshold value
function affect!(integrator)
    #Set state variables below the threshold to zero
    for i in 1:length(integrator.u)
        if integrator.u[i] < threshold
            integrator.u[i] = 0.0  
        end
    end
end

#Create the discrete callback
cb_extinction = DiscreteCallback(condition, affect!)


#Functions to calculate real and imaginary eigenvalues
#maximum(real.(eigvals(M))), where eigvals is from the standard library LinearAlgebra

λ_stability(M) = maximum(real.(eigvals(M)))
imag_eig(M) = maximum(imag(eigvals(M)))

calc_λ1(eq, par) = λ_stability(ForwardDiff.jacobian(eq -> CR_mod(eq, par), eq))

calc_λe(eq, par) = imag_eig(ForwardDiff.jacobian(eq -> CR_mod(eq, par), eq))


#Test that each component of local stability analysis works
par = CRPar()
tspan = (0.0, 1000.0)
u0 = [2.0, 1.0]

prob = ODEProblem(CR_mod!, u0, tspan, par)

sol = solve(prob, Tsit5(), callback = cb_extinction, reltol = 1e-8, abstol = 1e-8)

eq = nlsolve((du, u) -> CR_mod!(du, u, par, 0.0), sol.u[end]).zero

u0 = [1.5, 1.5]
t_span = (0.0, 1000.0)
p = CRPar(noise = 0.1, a = 3.5)

prob_stoch = SDEProblem(CR_mod!, stoch_CR_Mod!, u0, t_span, p)
sol_stoch = solve(prob_stoch, ImplicitEM(), callback = cb_extinction, reltol = 1e-6, abstol = 1e-6, maxiters = 1e8)

# Stochastic Time Series
let 
    u0 = [1.5, 1.5]
    t_span = (0.0, 1000.0)
    p = CRPar(noise = 0.01, a = 0.6)
    prob_stoch = SDEProblem(CR_mod!, stoch_CR_Mod!, u0, t_span, p)
    sol_stoch = solve(prob_stoch, ImplicitEM(), callback = cb_extinction, reltol = 1e-6, abstol = 1e-6, maxiters = 1e8)
    adapt_rozmacts = figure()
    plot(sol_stoch.t[1:end], sol_stoch[1, 1:end], label = "Resource")
    plot(sol_stoch.t[1:end], sol_stoch[2, 1:end], label = "Consumer")
    xlabel("time", fontsize = 14, fontname = "Arial")
    ylabel("Density", fontsize = 14, fontname = "Arial")
    legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=2)
    
    gca().tick_params(axis="both", which="major", labelsize=14) 
    gca().tick_params(axis="both", which="minor", labelsize=14) 

    ax = gca()

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)

    tight_layout()

    return adapt_rozmacts
end

\
###########################################################################
##### Investigate Aspects of stability over changing a (attack rate) ######
###########################################################################

##Real and Imaginary Eigenvalues
let
    avals = 0.58:0.01:1.05
    stab = fill(NaN, 2, length(avals))
    par = CRPar()

    for (i, a) in enumerate(avals)
        par.a = a
        print(par)
        prob = ODEProblem(CR_mod!, u0, tspan, par)
        sol = solve(prob, callback = cb_extinction, reltol = 1e-6, abstol = 1e-6)
        eq = nlsolve((du, u) -> CR_mod!(du, u, par, 0.0), sol.u[end]).zero
        stab[1, i] = calc_λ1(eq, par)
        stab[2, i] = calc_λe(eq, par)
    end

    CR_Damp_eq = figure(figsize = (8,6))
    #subplot(211)
    plot(avals, stab[1, :], color="black")
    axhline(0, color = "black", linestyle = "--")
    ylabel(L"\lambda_1",fontsize=18,fontweight=:bold, fontname = "Arial")
    
    gca().tick_params(axis="both", which="major", labelsize=16) 
    gca().tick_params(axis="both", which="minor", labelsize=16) 

    ax = gca()

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)

    #subplot(212)
    #plot(avals, stab[2, :])
    xlabel("Attack Rate (a)",fontsize=18,fontweight=:bold, fontname = "Arial")
    #ylabel(L"\lambda_ imag",fontsize=16,fontweight=:bold)
    
    tight_layout()

    #savefig("/Users/reillyoconnor/Desktop/R Projects/Fast-Slow-Stability/Julia/Figures/Eigenvalues_D0.jpeg", dpi = 300)

    return CR_Damp_eq
end

##Interaction Strength (Relative Energy Flux - Nillson & McCann 2015))
let
    avals = 0.58:0.001:1.05
    is_flux = fill(NaN, 1, length(avals))
    par = CRPar()

    for (i, a) in enumerate(avals)
        par.a = a
        print(par)
        prob = ODEProblem(CR_mod!, u0, tspan, par)
        sol = solve(prob, callback = cb_extinction, reltol = 1e-06, abstol = 1e-06)
        eq = nlsolve((du, u) -> CR_mod!(du, u, par, 0.0), sol.u[end]).zero

        R, C = eq 
    
        is = par.k * ((par.a * par.e) - par.m) / (par.m * par.b)

        is_flux[i] = is 
    end

    df_is_flux = DataFrame(rmax = avals, interaction_strength_flux = is_flux[1, :])
    CSV.write("/Users/reillyoconnor/Desktop/R Projects/Fast-Slow-Stability/Julia/CR_IS_Flux_rmax.csv", df_is_flux)

    CR_is_flux = figure(figsize = (8,6))

    scatter(avals, is_flux[1,:], color="black", s = 50, alpha = 0.15)
    ylabel("Interaction Strength (Per Capita Flux)",fontsize=18,fontweight=:bold, fontname = "Arial")
    xlabel("Attack Rate (a)",fontsize=18,fontweight=:bold, fontname = "Arial")
    
    gca().tick_params(axis="both", which="major", labelsize=16) 
    gca().tick_params(axis="both", which="minor", labelsize=16) 
    
    ax = gca()

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)
    xlim(0.5, 1.15)
    ylim(0.75, 5.25)

    tight_layout()

    savefig("/Users/reillyoconnor/Desktop/R Projects/Fast-Slow-Stability/Julia/Figures/IS_Flux_rmax.jpeg", dpi = 300)

    return CR_is_flux
end


/


let 
    amaxs = 0.58:0.01:1.1
    num_simulations = 20 # Number of simulations for averaging

    # Initialize arrays to store average results
    avg_stdhold = fill(0.0, length(amaxs), 1)
    avg_meanhold = fill(0.0, length(amaxs), 1)
    avg_cvhold = fill(0.0, length(amaxs), 1)
    #avg_k_cvhold = fill(0.0, length(amaxs), 1)
    avg_taylorhold = fill(0.0, length(amaxs), 1)

    # Time settings
    ts = range(0, 2000, length = 2000)
    t_span = (0.0, 2000.0)

    # Create directory for CSV files if it doesn't exist
    output_dir = "/Users/reillyoconnor/Desktop/R Projects/Slow_Fast_IS_Stability/Julia/Outputs/CR"
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    # Loop through each value of amax in parallel using threads
    @threads for i in 1:length(amaxs)


        # Temporary arrays to hold results for multiple simulations
        std_vals = zeros(num_simulations)
        mean_vals = zeros(num_simulations)
        cv_vals = zeros(num_simulations)
        #k_cv_vals = zeros(num_simulations)
        taylor_vals = zeros(num_simulations)

        # Run multiple simulations in parallel
        for j = 1:num_simulations
            # Set parameters for the current value of amax
            p = CRPar(a = amaxs[i], noise = 0.01)
            u0 = [1.0, 0.5]

            # Solve the stochastic problem
            prob_stoch = SDEProblem(CR_mod!, stoch_CR_Mod!, u0, t_span, p)
            sol_stoch = solve(prob_stoch, ImplicitEM(), callback = cb_extinction, reltol = 1e-6, abstol = 1e-6, maxiters = 1e8)
            grid_sol = sol_stoch(ts)

            # Calculate statistics for the consumer population (index 2)
            std_vals[j] = std(grid_sol[2, 1000:2000])
            mean_vals[j] = mean(grid_sol[2, 1000:2000])
            cv_vals[j] = std_vals[j] / mean_vals[j]
            #k_cv_vals[j] = sqrt((cv_vals[j]^2)/(1 + cv_vals[j]^2))
            taylor_vals[j] = (std_vals[j]^2) / mean_vals[j]

            # Create a DataFrame for the current run
            results_df = DataFrame(
                rmax_increase = [amaxs[i]],
                simulation_run = [j],
                mean = [mean_vals[j]],
                std = [std_vals[j]],
                cv = [cv_vals[j]],
                # k_cv = [k_cv_vals[j]],
                taylor = [taylor_vals[j]]
            )

            # Write to a CSV file for each run
            filename = "$(output_dir)/rmax_$(i)_run_$(j).csv"
            CSV.write(filename, results_df)
        end

        # Store the average values for each amax increment
        avg_stdhold[i] = mean(std_vals)
        avg_meanhold[i] = mean(mean_vals)
        avg_cvhold[i] = mean(cv_vals)
        # avg_k_cvhold[i] = mean(k_cv_vals)
        avg_taylorhold[i] = mean(taylor_vals)
    end

    println("Simulation results have been saved to the directory: $output_dir")
end


