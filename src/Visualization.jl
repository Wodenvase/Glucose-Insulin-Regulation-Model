"""
    Visualization

Plotting functions for glucose-insulin dynamics.

Provides functions for:
- Time series plots
- Phase portraits
- Bifurcation diagrams
- Comparison plots
"""

using Plots

"""
    plot_time_series(sol; title="Glucose-Insulin Dynamics", 
                     show_components=true, add_basal=true, params=nothing)

Plot time series of glucose, insulin, and insulin action.

# Arguments
- `sol`: ODE solution from DifferentialEquations.jl
- `title`: Plot title
- `show_components`: If true, plot all three state variables
- `add_basal`: If true, show total glucose (G + Gb) instead of G above basal
- `params`: BergmanParams for adding basal glucose

# Returns
Plots.jl plot object
"""
function plot_time_series(sol; title="Glucose-Insulin Dynamics", 
                          show_components=true, add_basal=true, params=nothing)
    
    t = sol.t
    G = [u[1] for u in sol.u]
    X = [u[2] for u in sol.u]
    I = [u[3] for u in sol.u]
    
    if add_basal && !isnothing(params)
        G_total = G .+ params.Gb
        ylabel_G = "Total Glucose (mg/dL)"
    else
        G_total = G
        ylabel_G = "Glucose above basal (mg/dL)"
    end
    
    if show_components
        p1 = plot(t, G_total, label="Glucose", lw=2, color=:blue,
                  xlabel="", ylabel=ylabel_G, legend=:topright, 
                  title=title)
        
        p2 = plot(t, I, label="Insulin", lw=2, color=:red,
                  xlabel="", ylabel="Insulin (μU/mL)", legend=:topright)
        
        p3 = plot(t, X, label="Insulin Action", lw=2, color=:green,
                  xlabel="Time (min)", ylabel="Insulin Action (min⁻¹)", 
                  legend=:topright)
        
        plot(p1, p2, p3, layout=(3,1), size=(800, 800))
    else
        plot(t, G_total, label="Glucose", lw=2, color=:blue,
             xlabel="Time (min)", ylabel=ylabel_G, legend=:topright,
             title=title, size=(800, 400))
    end
end

"""
    plot_phase_portrait(sol, params::BergmanParams; 
                        fixed_point=nothing, title="Phase Portrait (G vs I)")

Plot phase portrait in glucose-insulin space.

# Arguments
- `sol`: ODE solution
- `params`: Model parameters
- `fixed_point`: Optional FixedPoint to mark
- `title`: Plot title

# Returns
Plots.jl plot object
"""
function plot_phase_portrait(sol, params::BergmanParams; 
                             fixed_point=nothing, title="Phase Portrait (G vs I)")
    
    G = [u[1] for u in sol.u]
    I = [u[3] for u in sol.u]
    
    # Plot trajectory
    p = plot(G, I, lw=2, label="Trajectory", color=:blue,
             xlabel="Glucose above basal (mg/dL)", 
             ylabel="Insulin (μU/mL)",
             title=title, legend=:topright)
    
    # Mark initial and final points
    scatter!([G[1]], [I[1]], marker=:circle, markersize=8, 
             label="Start", color=:green)
    scatter!([G[end]], [I[end]], marker=:square, markersize=8,
             label="End", color=:red)
    
    # Mark fixed point if provided
    if !isnothing(fixed_point) && fixed_point.exists
        scatter!([fixed_point.state[1]], [fixed_point.state[3]], 
                 marker=:star, markersize=12, label="Fixed Point", 
                 color=:gold)
    end
    
    return p
end

"""
    plot_phase_plane_with_nullclines(sol, params::BergmanParams;
                                     D=0.0, I_ext=0.0, title="Phase Plane with Nullclines")

Plot phase plane with nullclines and flow field.

# Arguments
- `sol`: ODE solution
- `params`: Model parameters
- `D`: Glucose input function value (constant)
- `I_ext`: Exogenous insulin (constant)
- `title`: Plot title

# Returns
Plots.jl plot object
"""
function plot_phase_plane_with_nullclines(sol, params::BergmanParams;
                                          D=0.0, I_ext=0.0, 
                                          title="Phase Plane with Nullclines")
    
    G = [u[1] for u in sol.u]
    I = [u[3] for u in sol.u]
    
    # Create grid for nullclines
    G_range = range(minimum(G) - 10, maximum(G) + 10, length=100)
    I_range = range(max(0, minimum(I) - 5), maximum(I) + 10, length=100)
    
    # Glucose nullcline: dG/dt = 0
    # -k1*G - X*(G + Gb) + D = 0, where X = k3*I/k2
    G_nullcline_I = Float64[]
    G_nullcline_G = Float64[]
    for G_val in G_range
        # X = k3*I/k2, so: -k1*G - (k3*I/k2)*(G + Gb) + D = 0
        # Solve for I: I = k2*(D - k1*G) / (k3*(G + Gb))
        if abs(G_val + params.Gb) > 1e-6
            I_val = params.k2 * (D - params.k1 * G_val) / (params.k3 * (G_val + params.Gb))
            if I_val >= 0  # Physical constraint
                push!(G_nullcline_G, G_val)
                push!(G_nullcline_I, I_val)
            end
        end
    end
    
    # Insulin nullcline: dI/dt = 0
    # -k4*I + k5*max(G, 0) + I_ext = 0
    # I = k5*max(G, 0)/k4 + I_ext/k4
    I_nullcline_G = collect(G_range)
    I_nullcline_I = [params.k5 * max(G_val, 0.0) / params.k4 + I_ext / params.k4 
                     for G_val in I_nullcline_G]
    
    # Plot
    p = plot(G, I, lw=2, label="Trajectory", color=:black,
             xlabel="Glucose above basal (mg/dL)", 
             ylabel="Insulin (μU/mL)",
             title=title, legend=:topright)
    
    # Add nullclines
    plot!(G_nullcline_G, G_nullcline_I, lw=2, ls=:dash, 
          label="dG/dt = 0", color=:blue)
    plot!(I_nullcline_G, I_nullcline_I, lw=2, ls=:dash,
          label="dI/dt = 0", color=:red)
    
    # Mark trajectory points
    scatter!([G[1]], [I[1]], marker=:circle, markersize=8,
             label="Start", color=:green)
    scatter!([G[end]], [I[end]], marker=:square, markersize=8,
             label="End", color=:orange)
    
    return p
end

"""
    plot_bifurcation_diagram(bifurcation_results; 
                             variable_index=1,
                             param_name="Parameter",
                             title="Bifurcation Diagram")

Plot bifurcation diagram showing fixed points vs parameter.

# Arguments
- `bifurcation_results`: Results from bifurcation_diagram()
- `variable_index`: Which state variable to plot (1=G, 2=X, 3=I)
- `param_name`: Name of parameter for x-axis label
- `title`: Plot title

# Returns
Plots.jl plot object
"""
function plot_bifurcation_diagram(bifurcation_results; 
                                  variable_index=1,
                                  param_name="Parameter",
                                  title="Bifurcation Diagram")
    
    param_values = bifurcation_results.param_values
    fixed_points = bifurcation_results.fixed_points
    stable = bifurcation_results.stable
    
    # Extract variable values
    variable_values = [fp.exists ? fp.state[variable_index] : NaN 
                       for fp in fixed_points]
    
    # Separate stable and unstable
    stable_params = param_values[stable]
    stable_vals = variable_values[stable]
    unstable_params = param_values[.!stable]
    unstable_vals = variable_values[.!stable]
    
    # Variable name
    var_names = ["Glucose (mg/dL)", "Insulin Action (min⁻¹)", "Insulin (μU/mL)"]
    ylabel = var_names[variable_index]
    
    # Plot
    p = scatter(stable_params, stable_vals, label="Stable", 
                color=:blue, markersize=3, markerstrokewidth=0,
                xlabel=param_name, ylabel=ylabel, title=title,
                legend=:topright)
    scatter!(unstable_params, unstable_vals, label="Unstable",
             color=:red, markersize=3, markerstrokewidth=0)
    
    # Detect and mark Hopf bifurcations
    hopf_indices = detect_hopf_bifurcation(bifurcation_results)
    if !isempty(hopf_indices)
        hopf_params = param_values[hopf_indices]
        hopf_vals = variable_values[hopf_indices]
        scatter!(hopf_params, hopf_vals, marker=:star, markersize=10,
                 label="Hopf Bifurcation", color=:gold)
    end
    
    return p
end

"""
    plot_comparison(sols, labels; title="Comparison", add_basal=true, params=nothing)

Plot multiple solutions on the same axes for comparison.

# Arguments
- `sols`: Vector of ODE solutions
- `labels`: Vector of labels for each solution
- `title`: Plot title
- `add_basal`: If true, show total glucose
- `params`: BergmanParams for basal glucose

# Returns
Plots.jl plot object
"""
function plot_comparison(sols, labels; title="Comparison", 
                         add_basal=true, params=nothing)
    
    colors = [:blue, :red, :green, :purple, :orange, :brown]
    
    p = plot(title=title, xlabel="Time (min)", 
             ylabel=add_basal && !isnothing(params) ? "Total Glucose (mg/dL)" : "Glucose (mg/dL)",
             legend=:topright, size=(800, 500))
    
    for (i, (sol, label)) in enumerate(zip(sols, labels))
        t = sol.t
        G = [u[1] for u in sol.u]
        
        if add_basal && !isnothing(params)
            G = G .+ params.Gb
        end
        
        color = colors[mod1(i, length(colors))]
        plot!(t, G, label=label, lw=2, color=color)
    end
    
    return p
end

"""
    plot_daily_profile(sol, params::BergmanParams, meal_times;
                       title="Daily Glucose Profile")

Plot glucose profile throughout the day with meal markers.

# Arguments
- `sol`: ODE solution
- `params`: Model parameters
- `meal_times`: Vector of meal times in minutes
- `title`: Plot title

# Returns
Plots.jl plot object
"""
function plot_daily_profile(sol, params::BergmanParams, meal_times;
                            title="Daily Glucose Profile")
    
    t = sol.t
    G_total = [u[1] + params.Gb for u in sol.u]
    
    p = plot(t, G_total, lw=2, color=:blue, label="Glucose",
             xlabel="Time (min)", ylabel="Total Glucose (mg/dL)",
             title=title, legend=:topright, size=(1000, 400))
    
    # Add normal range shading
    hline!([70, 140], lw=1, ls=:dash, color=:green, 
           label="Normal Range", alpha=0.3)
    
    # Mark meal times
    for (i, t_meal) in enumerate(meal_times)
        vline!([t_meal], lw=2, ls=:dot, color=:red, 
               label=(i==1 ? "Meals" : ""))
    end
    
    return p
end

"""
    plot_eigenvalue_trajectory(bifurcation_results; title="Eigenvalue Trajectory")

Plot eigenvalue trajectories in complex plane as parameter varies.

# Arguments
- `bifurcation_results`: Results from bifurcation_diagram()
- `title`: Plot title

# Returns
Plots.jl plot object
"""
function plot_eigenvalue_trajectory(bifurcation_results; 
                                    title="Eigenvalue Trajectory")
    
    eigenvalues_list = bifurcation_results.eigenvalues
    
    p = plot(xlabel="Real Part", ylabel="Imaginary Part",
             title=title, legend=:topright, aspect_ratio=:equal)
    
    # Plot imaginary axis (stability boundary)
    ymax = 0.1
    plot!([0, 0], [-ymax, ymax], lw=2, color=:black, ls=:dash,
          label="Stability Boundary")
    
    # Extract eigenvalue trajectories
    n_eigenvalues = length(eigenvalues_list[1])
    
    for i in 1:n_eigenvalues
        real_parts = [real(λs[i]) for λs in eigenvalues_list if length(λs) >= i]
        imag_parts = [imag(λs[i]) for λs in eigenvalues_list if length(λs) >= i]
        
        plot!(real_parts, imag_parts, lw=2, marker=:circle, markersize=2,
              label="λ$i", alpha=0.7)
    end
    
    return p
end

"""
    plot_input_output_curve(param_range, base_params::BergmanParams, 
                            input_param::Symbol; 
                            output_var::Symbol=:G,
                            title="Input-Output Curve")

Plot steady-state output vs input parameter.

# Arguments
- `param_range`: Range of input parameter values
- `base_params`: Base parameters
- `input_param`: Parameter to vary (e.g., :D for glucose input)
- `output_var`: Output variable to plot (:G, :X, or :I)
- `title`: Plot title

# Returns
Plots.jl plot object
"""
function plot_input_output_curve(param_range, base_params::BergmanParams, 
                                 input_param::Symbol; 
                                 output_var::Symbol=:G,
                                 title="Input-Output Curve")
    
    # For external inputs (D, I_ext), use different approach
    if input_param == :D
        D_values = collect(param_range)
        output_values = Float64[]
        
        for D in D_values
            fp = find_fixed_point(base_params, D)
            if fp.exists
                idx = output_var == :G ? 1 : (output_var == :X ? 2 : 3)
                push!(output_values, fp.state[idx])
            else
                push!(output_values, NaN)
            end
        end
        
        xlabel_str = "Glucose Input D (mg/dL/min)"
        
    elseif input_param == :I_ext
        I_ext_values = collect(param_range)
        output_values = Float64[]
        
        for I_ext in I_ext_values
            fp = find_fixed_point(base_params, 0.0, I_ext=I_ext)
            if fp.exists
                idx = output_var == :G ? 1 : (output_var == :X ? 2 : 3)
                push!(output_values, fp.state[idx])
            else
                push!(output_values, NaN)
            end
        end
        
        xlabel_str = "Exogenous Insulin (μU/mL)"
        D_values = I_ext_values  # For plotting
    else
        error("Unsupported input parameter for input-output curve")
    end
    
    ylabel_str = output_var == :G ? "Glucose (mg/dL)" :
                 (output_var == :X ? "Insulin Action (min⁻¹)" : "Insulin (μU/mL)")
    
    plot(D_values, output_values, lw=2, marker=:circle,
         xlabel=xlabel_str, ylabel=ylabel_str, title=title,
         legend=false, size=(800, 500))
end

"""
    plot_stability_region(param1_range, param2_range, 
                          param1_name::Symbol, param2_name::Symbol,
                          base_params::BergmanParams;
                          title="Stability Region")

Plot 2D stability region in parameter space.

# Arguments
- `param1_range`, `param2_range`: Parameter ranges
- `param1_name`, `param2_name`: Parameter names (:k1, :k2, etc.)
- `base_params`: Base parameters
- `title`: Plot title

# Returns
Plots.jl plot object
"""
function plot_stability_region(param1_range, param2_range, 
                               param1_name::Symbol, param2_name::Symbol,
                               base_params::BergmanParams;
                               title="Stability Region")
    
    p1_vals = collect(param1_range)
    p2_vals = collect(param2_range)
    
    n1 = length(p1_vals)
    n2 = length(p2_vals)
    
    stability_matrix = zeros(n2, n1)
    
    for (i, p1) in enumerate(p1_vals)
        for (j, p2) in enumerate(p2_vals)
            # Create modified parameters
            params_dict = Dict(
                :k1 => base_params.k1,
                :k2 => base_params.k2,
                :k3 => base_params.k3,
                :k4 => base_params.k4,
                :k5 => base_params.k5,
                :Gb => base_params.Gb
            )
            params_dict[param1_name] = p1
            params_dict[param2_name] = p2
            
            params = BergmanParams(
                k1=params_dict[:k1],
                k2=params_dict[:k2],
                k3=params_dict[:k3],
                k4=params_dict[:k4],
                k5=params_dict[:k5],
                Gb=params_dict[:Gb]
            )
            
            # Check stability
            fp = find_fixed_point(params, 0.0)
            if fp.exists
                analysis = analyze_stability(fp, params)
                stability_matrix[j, i] = analysis.is_stable ? 1.0 : 0.0
            else
                stability_matrix[j, i] = -1.0  # No fixed point
            end
        end
    end
    
    heatmap(p1_vals, p2_vals, stability_matrix,
            xlabel=String(param1_name), ylabel=String(param2_name),
            title=title, color=:viridis, colorbar_title="Stability",
            size=(700, 600))
end

export plot_time_series, plot_phase_portrait, plot_phase_plane_with_nullclines
export plot_bifurcation_diagram, plot_comparison, plot_daily_profile
export plot_eigenvalue_trajectory, plot_input_output_curve, plot_stability_region
