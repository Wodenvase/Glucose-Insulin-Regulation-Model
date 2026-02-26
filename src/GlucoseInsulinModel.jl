"""
# GlucoseInsulinModel

Main module for glucose-insulin regulation modeling using the Bergman minimal model.

This package provides tools for:
- Simulating glucose-insulin dynamics
- Analyzing stability and bifurcations
- Modeling diabetes and therapeutic interventions
- Visualizing results

# Usage

using Plots
include("src/GlucoseInsulinModel.jl")
using .GlucoseInsulinModel

# Define parameters (healthy adult)
params = BergmanParams()

# Oral glucose tolerance test
D_func = create_glucose_input(:ogtt)
I_ext_func = create_insulin_input(:none)

# Simulate
u0 = [0.0, 0.0, 0.0]  # Start at baseline
tspan = (0.0, 180.0)   # 3 hours
sol = simulate_bergman(u0, tspan, params, D_func, I_ext_func)

# Visualize
plot(sol)

# Modules
- `BergmanModel`: Core ODE system and parameters
- `AnalysisTools`: Fixed point and stability analysis
- `Visualization`: Plotting functions

# References
- Bergman et al. (1979) J Clin Invest 68:1456-1467
- Bergman et al. (1981) Endocr Rev 6:45-86
"""
module GlucoseInsulinModel

using DifferentialEquations
using Plots

# Include submodules (use absolute paths relative to this file)
include(joinpath(@__DIR__, "BergmanModel.jl"))
include(joinpath(@__DIR__, "AnalysisTools.jl"))
include(joinpath(@__DIR__, "Visualization.jl"))

# Submodules have been included above; their symbols are available in this module's namespace

export BergmanParams, bergman_dynamics!, bergman_dynamics
export insulin_sensitivity, disposition_index, time_constants
export glucose_disposal_rate, steady_state_glucose
export create_glucose_input, create_insulin_input

export FixedPoint, StabilityAnalysis
export find_fixed_point, find_multiple_fixed_points
export compute_jacobian, classify_stability, analyze_stability
export bifurcation_diagram, detect_hopf_bifurcation
export sensitivity_analysis

export plot_time_series, plot_phase_portrait, plot_phase_plane_with_nullclines
export plot_bifurcation_diagram, plot_comparison, plot_daily_profile
export plot_eigenvalue_trajectory, plot_input_output_curve, plot_stability_region

"""
    simulate_bergman(u0, tspan, params::BergmanParams, 
                     D_func, I_ext_func=t->0.0;
                     solver=Tsit5(), kwargs...)

Simulate the Bergman minimal model.

# Arguments
- `u0`: Initial state [G, X, I]
- `tspan`: Time span (t_start, t_end)
- `params`: Model parameters
- `D_func`: Glucose input function D(t)
- `I_ext_func`: Exogenous insulin function I_ext(t) (default: no insulin)
- `solver`: ODE solver (default: Tsit5)
- `kwargs`: Additional arguments for solve()

# Returns
ODE solution object

# Examples
# OGTT simulation
params = BergmanParams()
D = create_glucose_input(:ogtt)
I_ext = create_insulin_input(:none)

sol = simulate_bergman([0.0, 0.0, 0.0], (0.0, 180.0), params, D, I_ext)
"""
function simulate_bergman(u0, tspan, params::BergmanParams, 
                          D_func, I_ext_func=t->0.0;
                          solver=Tsit5(), kwargs...)
    
    # Create ODE problem
    p = (params, D_func, I_ext_func)
    prob = ODEProblem(bergman_dynamics!, u0, tspan, p)
    
    # Solve
    sol = solve(prob, solver; kwargs...)
    
    return sol
end

"""
    simulate_with_callback(u0, tspan, params::BergmanParams,
                           D_func, I_ext_func=t->0.0;
                           hypoglycemia_threshold=70.0,
                           hyperglycemia_threshold=250.0)

Simulate with event detection for hypo/hyperglycemia.

# Arguments
- `u0`: Initial state
- `tspan`: Time span
- `params`: Model parameters
- `D_func`: Glucose input function
- `I_ext_func`: Exogenous insulin function
- `hypoglycemia_threshold`: Low glucose threshold (mg/dL, total)
- `hyperglycemia_threshold`: High glucose threshold (mg/dL, total)

# Returns
Tuple (sol, events) where events is Dict with :hypoglycemia and :hyperglycemia times
"""
function simulate_with_callback(u0, tspan, params::BergmanParams,
                                D_func, I_ext_func=t->0.0;
                                hypoglycemia_threshold=70.0,
                                hyperglycemia_threshold=250.0)
    
    events = Dict(:hypoglycemia => Float64[], :hyperglycemia => Float64[])
    
    # Hypoglycemia callback
    function hypo_condition(u, t, integrator)
        G_total = u[1] + params.Gb
        return G_total - hypoglycemia_threshold
    end
    
    function hypo_affect!(integrator)
        push!(events[:hypoglycemia], integrator.t)
        @warn "Hypoglycemia detected at t=$(integrator.t) min"
    end
    
    hypo_cb = ContinuousCallback(hypo_condition, hypo_affect!)
    
    # Hyperglycemia callback
    function hyper_condition(u, t, integrator)
        G_total = u[1] + params.Gb
        return hyperglycemia_threshold - G_total
    end
    
    function hyper_affect!(integrator)
        push!(events[:hyperglycemia], integrator.t)
        @warn "Hyperglycemia detected at t=$(integrator.t) min"
    end
    
    hyper_cb = ContinuousCallback(hyper_condition, hyper_affect!)
    
    # Combine callbacks
    cb = CallbackSet(hypo_cb, hyper_cb)
    
    # Simulate
    p = (params, D_func, I_ext_func)
    prob = ODEProblem(bergman_dynamics!, u0, tspan, p)
    sol = solve(prob, Tsit5(), callback=cb)
    
    return sol, events
end

"""
    create_disease_params(disease_type::Symbol; severity=:moderate)

Create parameter set representing different disease states.

# Disease types
- `:healthy`: Normal parameters
- `:type1`: Type 1 diabetes (k5 = 0)
- `:type2_early`: Early Type 2 (reduced k3, compensatory k5)
- `:type2_late`: Late Type 2 (reduced k3 and k5)
- `:prediabetes`: Impaired glucose tolerance

# Severity
- `:mild`: Small parameter changes
- `:moderate`: Moderate changes (default)
- `:severe`: Large changes

# Returns
BergmanParams structure

# Examples
params_t1d = create_disease_params(:type1)
params_t2d = create_disease_params(:type2_late, severity=:severe)
"""
function create_disease_params(disease_type::Symbol; severity=:moderate)
    
    # Base healthy parameters
    k1_base = 0.028
    k2_base = 0.025
    k3_base = 2.5e-5
    k4_base = 0.05
    k5_base = 0.015
    Gb_base = 81.0
    
    # Severity factors
    if severity == :mild
        factor = 0.8
    elseif severity == :moderate
        factor = 0.6
    elseif severity == :severe
        factor = 0.4
    else
        error("Unknown severity: $severity")
    end
    
    if disease_type == :healthy
        return BergmanParams()
        
    elseif disease_type == :type1
        # Complete insulin deficiency
        return BergmanParams(k1=k1_base, k2=k2_base, k3=k3_base,
                            k4=k4_base, k5=0.0, Gb=Gb_base)
        
    elseif disease_type == :type2_early
        # Insulin resistance with compensation
        k3 = k3_base * factor
        k5 = k5_base * (1.5 - 0.5 * (1 - factor))  # Compensatory increase
        return BergmanParams(k1=k1_base, k2=k2_base, k3=k3,
                            k4=k4_base, k5=k5, Gb=Gb_base)
        
    elseif disease_type == :type2_late
        # Insulin resistance with β-cell failure
        k3 = k3_base * factor
        k5 = k5_base * factor
        return BergmanParams(k1=k1_base, k2=k2_base, k3=k3,
                            k4=k4_base, k5=k5, Gb=Gb_base)
        
    elseif disease_type == :prediabetes
        # Mild insulin resistance, normal secretion
        k3 = k3_base * 0.7
        return BergmanParams(k1=k1_base, k2=k2_base, k3=k3,
                            k4=k4_base, k5=k5_base, Gb=Gb_base)
        
    else
        error("Unknown disease type: $disease_type")
    end
end

"""
    simulate_ogtt(params::BergmanParams; duration=180.0, glucose_dose=75.0)

Simulate oral glucose tolerance test (OGTT).

Standard clinical test: 75g glucose over 5 minutes, measure for 2-3 hours.

# Arguments
- `params`: Model parameters
- `duration`: Test duration in minutes (default 180 = 3 hours)
- `glucose_dose`: Glucose dose in grams (default 75g)

# Returns
Tuple (sol, diagnosis) where diagnosis is :normal, :prediabetes, or :diabetes

# Diagnostic criteria (2-hour glucose)
- Normal: < 140 mg/dL
- Prediabetes: 140-199 mg/dL
- Diabetes: ≥ 200 mg/dL
"""
function simulate_ogtt(params::BergmanParams; duration=180.0, glucose_dose=75.0)
    
    # Create OGTT glucose input (75g over 5 min)
    # 1g carb ≈ 4 mg/dL, so 75g ≈ 300 mg/dL spread over 5 min
    D = create_glucose_input(:ogtt)
    I_ext = create_insulin_input(:none)
    
    # Simulate
    u0 = [0.0, 0.0, 0.0]
    tspan = (0.0, duration)
    sol = simulate_bergman(u0, tspan, params, D, I_ext)
    
    # Diagnose based on 2-hour glucose
    idx_120min = findfirst(t -> t >= 120.0, sol.t)
    if !isnothing(idx_120min)
        G_120 = sol.u[idx_120min][1] + params.Gb
        
        if G_120 < 140.0
            diagnosis = :normal
        elseif G_120 < 200.0
            diagnosis = :prediabetes
        else
            diagnosis = :diabetes
        end
    else
        diagnosis = :unknown
    end
    
    return sol, diagnosis
end

"""
    optimize_insulin_therapy(params::BergmanParams, target_glucose=100.0;
                             D_func, optimization_time=180.0)

Optimize basal-bolus insulin regimen to achieve target glucose.

Simple optimization: try different basal and bolus levels.

# Arguments
- `params`: Model parameters (typically Type 1 diabetes)
- `target_glucose`: Target glucose level (mg/dL total)
- `D_func`: Meal glucose input function
- `optimization_time`: Time to evaluate control (minutes)

# Returns
Dict with optimal basal and bolus insulin levels
"""
function optimize_insulin_therapy(params::BergmanParams, target_glucose=100.0;
                                  D_func, optimization_time=180.0)
    
    best_error = Inf
    best_basal = 0.0
    best_bolus = 0.0
    
    # Grid search
    for basal in 5.0:2.5:20.0
        for bolus in 20.0:10.0:80.0
            # Create insulin input
            I_ext = create_insulin_input(:basal_bolus, basal, [(0.0, bolus)])
            
            # Simulate
            u0 = [0.0, 0.0, 0.0]
            tspan = (0.0, optimization_time)
            sol = simulate_bergman(u0, tspan, params, D_func, I_ext)
            
            # Compute error (MSE from target)
            error = 0.0
            for u in sol.u
                G_total = u[1] + params.Gb
                error += (G_total - target_glucose)^2
            end
            error /= length(sol.u)
            
            if error < best_error
                best_error = error
                best_basal = basal
                best_bolus = bolus
            end
        end
    end
    
    return Dict(:basal => best_basal, :bolus => best_bolus, :error => best_error)
end

# Export main simulation functions
export simulate_bergman, simulate_with_callback
export create_disease_params, simulate_ogtt, optimize_insulin_therapy

end # module
