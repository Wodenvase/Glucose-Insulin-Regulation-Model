"""
    AnalysisTools

Tools for analyzing glucose-insulin dynamics:
- Fixed point analysis
- Stability analysis (Jacobian, eigenvalues)
- Bifurcation detection
- Sensitivity analysis
"""

using LinearAlgebra
using ForwardDiff
using NLsolve

"""
    FixedPoint

Structure representing a fixed point of the system.

# Fields
- `state`: [G*, X*, I*] equilibrium values
- `exists`: Whether fixed point was successfully found
- `residual`: Residual norm (should be near zero)
"""
struct FixedPoint
    state::Vector{Float64}
    exists::Bool
    residual::Float64
end

"""
    StabilityAnalysis

Results of stability analysis at a fixed point.

# Fields
- `fixed_point`: The fixed point being analyzed
- `jacobian`: Jacobian matrix at fixed point
- `eigenvalues`: Eigenvalues of Jacobian
- `is_stable`: Whether all eigenvalues have negative real part
- `stability_type`: Classification (stable_node, stable_spiral, etc.)
"""
struct StabilityAnalysis
    fixed_point::FixedPoint
    jacobian::Matrix{Float64}
    eigenvalues::Vector{ComplexF64}
    is_stable::Bool
    stability_type::Symbol
end

"""
    find_fixed_point(params::BergmanParams, D::Float64=0.0; 
                     initial_guess=[0.0, 0.0, 0.0])

Find fixed point of the system for constant glucose input D.

Solves:
    -k1*G - X*(G + Gb) + D = 0
    -k2*X + k3*I = 0
    -k4*I + k5*max(G, 0) = 0

# Arguments
- `params`: Model parameters
- `D`: Constant glucose input [mg/dL/min]
- `initial_guess`: Starting point for solver

# Returns
FixedPoint structure

# Examples
params = BergmanParams()

# Basal state (no input)
fp_basal = find_fixed_point(params, 0.0)

# With constant infusion
fp_infusion = find_fixed_point(params, 50.0)
"""
function find_fixed_point(params::BergmanParams, D::Float64=0.0; 
                          initial_guess=[0.0, 0.0, 0.0],
                          I_ext::Float64=0.0)
    
    function residual!(F, u)
        G, X, I = u
        k1, k2, k3, k4, k5, Gb = params.k1, params.k2, params.k3, 
                                  params.k4, params.k5, params.Gb
        
        # Fixed point conditions: du/dt = 0
        glucose_stimulus = max(G, 0.0)
        
        F[1] = -k1 * G - X * (G + Gb) + D
        F[2] = -k2 * X + k3 * I
        F[3] = -k4 * I + k5 * glucose_stimulus + I_ext
        
        return F
    end
    
    # Solve nonlinear system
    result = nlsolve(residual!, initial_guess, ftol=1e-9, iterations=1000)
    
    if converged(result)
        return FixedPoint(result.zero, true, result.residual_norm)
    else
        @warn "Fixed point solver did not converge"
        return FixedPoint(result.zero, false, result.residual_norm)
    end
end

"""
    find_multiple_fixed_points(params::BergmanParams, D::Float64=0.0;
                                n_guesses::Int=10)

Search for multiple fixed points using different initial guesses.

Systems with multiple fixed points can exhibit bistability.

# Arguments
- `params`: Model parameters
- `D`: Constant glucose input
- `n_guesses`: Number of random initial guesses

# Returns
Vector of unique FixedPoint structures
"""
function find_multiple_fixed_points(params::BergmanParams, D::Float64=0.0;
                                     n_guesses::Int=10, I_ext::Float64=0.0)
    
    fixed_points = FixedPoint[]
    
    # Try origin
    fp = find_fixed_point(params, D, initial_guess=[0.0, 0.0, 0.0], I_ext=I_ext)
    if fp.exists
        push!(fixed_points, fp)
    end
    
    # Try random initial guesses
    for i in 1:n_guesses
        guess = [rand() * 200.0 - 100.0,  # G: -100 to 100
                 rand() * 0.1,              # X: 0 to 0.1
                 rand() * 100.0]            # I: 0 to 100
        
        fp = find_fixed_point(params, D, initial_guess=guess, I_ext=I_ext)
        
        if fp.exists
            # Check if this is a new fixed point (not already found)
            is_new = true
            for existing_fp in fixed_points
                if norm(fp.state - existing_fp.state) < 1e-3
                    is_new = false
                    break
                end
            end
            
            if is_new
                push!(fixed_points, fp)
            end
        end
    end
    
    return fixed_points
end

"""
    compute_jacobian(fixed_point::FixedPoint, params::BergmanParams)

Compute Jacobian matrix at a fixed point.

The Jacobian J has elements J[i,j] = ∂(du[i]/dt)/∂u[j]

# Arguments
- `fixed_point`: Fixed point to analyze
- `params`: Model parameters

# Returns
3×3 Jacobian matrix
"""
function compute_jacobian(fixed_point::FixedPoint, params::BergmanParams)
    G_star, X_star, I_star = fixed_point.state
    k1, k2, k3, k4, k5, Gb = params.k1, params.k2, params.k3, 
                              params.k4, params.k5, params.Gb
    
    # Jacobian matrix (computed analytically)
    J = zeros(3, 3)
    
    # Row 1: dG/dt derivatives
    J[1, 1] = -k1 - X_star  # ∂(dG/dt)/∂G
    J[1, 2] = -(G_star + Gb)  # ∂(dG/dt)/∂X
    J[1, 3] = 0.0  # ∂(dG/dt)/∂I
    
    # Row 2: dX/dt derivatives
    J[2, 1] = 0.0  # ∂(dX/dt)/∂G
    J[2, 2] = -k2  # ∂(dX/dt)/∂X
    J[2, 3] = k3  # ∂(dX/dt)/∂I
    
    # Row 3: dI/dt derivatives
    # Note: d/dG[max(G, 0)] = Heaviside(G) ≈ 1 for G > 0, 0 for G < 0
    if G_star > 1e-6
        J[3, 1] = k5  # ∂(dI/dt)/∂G (when G > 0)
    else
        J[3, 1] = 0.0  # ∂(dI/dt)/∂G (when G ≤ 0)
    end
    J[3, 2] = 0.0  # ∂(dI/dt)/∂X
    J[3, 3] = -k4  # ∂(dI/dt)/∂I
    
    return J
end

"""
    classify_stability(eigenvalues::Vector{ComplexF64})

Classify fixed point stability based on eigenvalues.

# Eigenvalue criteria
- All real parts < 0: Stable
- Any real part > 0: Unstable
- Complex eigenvalues: Oscillatory (spiral)
- Real eigenvalues: Non-oscillatory (node)

# Returns
Symbol indicating stability type:
- `:stable_node`: All real, all negative
- `:stable_spiral`: Complex with negative real parts
- `:unstable_node`: All real, at least one positive
- `:unstable_spiral`: Complex with positive real parts
- `:saddle`: Mixed signs (real parts)
- `:center`: Purely imaginary (Re = 0)
"""
function classify_stability(eigenvalues::Vector{ComplexF64})
    real_parts = real.(eigenvalues)
    imag_parts = imag.(eigenvalues)
    
    all_stable = all(real_parts .< -1e-6)
    any_unstable = any(real_parts .> 1e-6)
    has_complex = any(abs.(imag_parts) .> 1e-6)
    
    if all_stable
        if has_complex
            return :stable_spiral
        else
            return :stable_node
        end
    elseif any_unstable
        if has_complex
            return :unstable_spiral
        else
            return :unstable_node
        end
    elseif any(real_parts .> 1e-6) && any(real_parts .< -1e-6)
        return :saddle
    else
        return :center
    end
end

"""
    analyze_stability(fixed_point::FixedPoint, params::BergmanParams)

Perform complete stability analysis at a fixed point.

Computes:
1. Jacobian matrix
2. Eigenvalues
3. Stability classification

# Arguments
- `fixed_point`: Fixed point to analyze
- `params`: Model parameters

# Returns
StabilityAnalysis structure

# Examples
params = BergmanParams()
fp = find_fixed_point(params)
analysis = analyze_stability(fp, params)

println("Stable: ", analysis.is_stable)
println("Type: ", analysis.stability_type)
println("Eigenvalues: ", analysis.eigenvalues)
"""
function analyze_stability(fixed_point::FixedPoint, params::BergmanParams)
    if !fixed_point.exists
        @warn "Cannot analyze stability of non-existent fixed point"
        return StabilityAnalysis(
            fixed_point,
            zeros(3, 3),
            ComplexF64[],
            false,
            :nonexistent
        )
    end
    
    # Compute Jacobian
    J = compute_jacobian(fixed_point, params)
    
    # Compute eigenvalues and ensure complex type for consistency
    λ_raw = eigvals(J)
    λ = ComplexF64.(λ_raw)

    # Classify stability
    stability_type = classify_stability(λ)
    is_stable = (stability_type == :stable_node || stability_type == :stable_spiral)

    return StabilityAnalysis(fixed_point, J, λ, is_stable, stability_type)
end

"""
    bifurcation_diagram(param_name::Symbol, param_range, 
                        base_params::BergmanParams;
                        D::Float64=0.0, I_ext::Float64=0.0)

Generate bifurcation diagram by varying a single parameter.

Tracks fixed points and their stability as parameter changes.

# Arguments
- `param_name`: Symbol for parameter to vary (:k1, :k2, :k3, :k4, :k5, :Gb)
- `param_range`: Range of parameter values
- `base_params`: Base parameter set
- `D`: Constant glucose input
- `I_ext`: Constant exogenous insulin

# Returns
Named tuple with:
- `param_values`: Parameter values scanned
- `fixed_points`: Fixed points at each parameter value
- `stability`: Stability at each point
- `eigenvalues`: Eigenvalues at each point

# Examples
params = BergmanParams()
results = bifurcation_diagram(:k5, 0.0:0.001:0.03, params)

# Plot bifurcation diagram
using Plots
plot(results.param_values, [fp.state[1] for fp in results.fixed_points])
"""
function bifurcation_diagram(param_name::Symbol, param_range, 
                             base_params::BergmanParams;
                             D::Float64=0.0, I_ext::Float64=0.0)
    
    param_values = collect(param_range)
    n_points = length(param_values)
    
    fixed_points = Vector{FixedPoint}(undef, n_points)
    stabilities = Vector{Bool}(undef, n_points)
    eigenvalues_list = Vector{Vector{ComplexF64}}(undef, n_points)
    stability_types = Vector{Symbol}(undef, n_points)
    
    for (i, param_value) in enumerate(param_values)
        # Create modified parameters
        if param_name == :k1
            params = BergmanParams(k1=param_value, k2=base_params.k2, k3=base_params.k3,
                                   k4=base_params.k4, k5=base_params.k5, Gb=base_params.Gb)
        elseif param_name == :k2
            params = BergmanParams(k1=base_params.k1, k2=param_value, k3=base_params.k3,
                                   k4=base_params.k4, k5=base_params.k5, Gb=base_params.Gb)
        elseif param_name == :k3
            params = BergmanParams(k1=base_params.k1, k2=base_params.k2, k3=param_value,
                                   k4=base_params.k4, k5=base_params.k5, Gb=base_params.Gb)
        elseif param_name == :k4
            params = BergmanParams(k1=base_params.k1, k2=base_params.k2, k3=base_params.k3,
                                   k4=param_value, k5=base_params.k5, Gb=base_params.Gb)
        elseif param_name == :k5
            params = BergmanParams(k1=base_params.k1, k2=base_params.k2, k3=base_params.k3,
                                   k4=base_params.k4, k5=param_value, Gb=base_params.Gb)
        elseif param_name == :Gb
            params = BergmanParams(k1=base_params.k1, k2=base_params.k2, k3=base_params.k3,
                                   k4=base_params.k4, k5=base_params.k5, Gb=param_value)
        else
            error("Unknown parameter: $param_name")
        end
        
        # Find fixed point
        initial_guess = i > 1 ? fixed_points[i-1].state : [0.0, 0.0, 0.0]
        fp = find_fixed_point(params, D, initial_guess=initial_guess, I_ext=I_ext)
        
        # Analyze stability
        if fp.exists
            analysis = analyze_stability(fp, params)
            fixed_points[i] = fp
            stabilities[i] = analysis.is_stable
            eigenvalues_list[i] = analysis.eigenvalues
            stability_types[i] = analysis.stability_type
        else
            fixed_points[i] = fp
            stabilities[i] = false
            eigenvalues_list[i] = ComplexF64[]
            stability_types[i] = :nonexistent
        end
    end
    
    return (param_values=param_values, 
            fixed_points=fixed_points, 
            stable=stabilities,
            eigenvalues=eigenvalues_list,
            stability_types=stability_types)
end

"""
    detect_hopf_bifurcation(bifurcation_results)

Detect Hopf bifurcation points in bifurcation diagram.

A Hopf bifurcation occurs when:
- Complex conjugate eigenvalues cross imaginary axis
- Re(λ) changes sign from negative to positive
- Im(λ) ≠ 0

# Arguments
- `bifurcation_results`: Results from bifurcation_diagram()

# Returns
Indices where Hopf bifurcations occur
"""
function detect_hopf_bifurcation(bifurcation_results)
    eigenvalues_list = bifurcation_results.eigenvalues
    n_points = length(eigenvalues_list)
    
    hopf_indices = Int[]
    
    for i in 2:n_points
        λ_prev = eigenvalues_list[i-1]
        λ_curr = eigenvalues_list[i]
        
        if isempty(λ_prev) || isempty(λ_curr)
            continue
        end
        
        # Check each pair of eigenvalues
        for (λp, λc) in zip(λ_prev, λ_curr)
            # Check for complex eigenvalues crossing imaginary axis
            if abs(imag(λp)) > 1e-6 && abs(imag(λc)) > 1e-6
                if real(λp) < 0 && real(λc) > 0
                    push!(hopf_indices, i)
                    break
                end
            end
        end
    end
    
    return unique(hopf_indices)
end

"""
    sensitivity_analysis(params::BergmanParams, perturbation::Float64=0.01)

Perform local sensitivity analysis around parameter values.

Computes how glucose response changes with small parameter perturbations.

# Arguments
- `params`: Nominal parameters
- `perturbation`: Relative perturbation (default 1%)

# Returns
Named tuple with sensitivity coefficients for each parameter
"""
function sensitivity_analysis(params::BergmanParams, perturbation::Float64=0.01)
    # Find nominal fixed point
    fp_nominal = find_fixed_point(params, 0.0)
    
    if !fp_nominal.exists
        @warn "Nominal fixed point not found"
        return nothing
    end
    
    G_nominal = fp_nominal.state[1]
    
    sensitivities = Dict{Symbol, Float64}()
    
    # Perturb each parameter
    for param_name in [:k1, :k2, :k3, :k4, :k5, :Gb]
        param_value = getfield(params, param_name)
        perturbed_value = param_value * (1 + perturbation)
        
        # Create perturbed parameters
        perturbed_params = if param_name == :k1
            BergmanParams(k1=perturbed_value, k2=params.k2, k3=params.k3,
                         k4=params.k4, k5=params.k5, Gb=params.Gb)
        elseif param_name == :k2
            BergmanParams(k1=params.k1, k2=perturbed_value, k3=params.k3,
                         k4=params.k4, k5=params.k5, Gb=params.Gb)
        elseif param_name == :k3
            BergmanParams(k1=params.k1, k2=params.k2, k3=perturbed_value,
                         k4=params.k4, k5=params.k5, Gb=params.Gb)
        elseif param_name == :k4
            BergmanParams(k1=params.k1, k2=params.k2, k3=params.k3,
                         k4=perturbed_value, k5=params.k5, Gb=params.Gb)
        elseif param_name == :k5
            BergmanParams(k1=params.k1, k2=params.k2, k3=params.k3,
                         k4=params.k4, k5=perturbed_value, Gb=params.Gb)
        else  # :Gb
            BergmanParams(k1=params.k1, k2=params.k2, k3=params.k3,
                         k4=params.k4, k5=params.k5, Gb=perturbed_value)
        end
        
        # Find perturbed fixed point
        fp_perturbed = find_fixed_point(perturbed_params, 0.0)
        
        if fp_perturbed.exists
            G_perturbed = fp_perturbed.state[1]
            
            # Normalized sensitivity: (ΔG/G) / (Δp/p)
            sensitivity = ((G_perturbed - G_nominal) / G_nominal) / perturbation
            sensitivities[param_name] = sensitivity
        else
            sensitivities[param_name] = NaN
        end
    end
    
    return sensitivities
end

export FixedPoint, StabilityAnalysis
export find_fixed_point, find_multiple_fixed_points
export compute_jacobian, classify_stability, analyze_stability
export bifurcation_diagram, detect_hopf_bifurcation
export sensitivity_analysis
