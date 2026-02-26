"""
    BergmanModel

Core Bergman minimal model implementation for glucose-insulin dynamics.

Implements the three-compartment minimal model:
- G(t): Plasma glucose above basal (mg/dL)
- X(t): Insulin action on glucose uptake (min⁻¹)
- I(t): Plasma insulin above basal (μU/mL)

References:
- Bergman et al. (1979) J Clin Invest 68:1456-1467
- Bergman et al. (1981) Endocr Rev 6:45-86
"""

using LinearAlgebra

"""
    BergmanParams

Parameters for the Bergman minimal model.

# Fields
- `k1`: Glucose effectiveness (insulin-independent uptake) [min⁻¹]
- `k2`: Insulin action decay rate [min⁻¹]
- `k3`: Insulin sensitivity parameter [min⁻² per μU/mL]
- `k4`: Insulin clearance rate [min⁻¹]
- `k5`: Pancreatic responsiveness (β-cell secretion rate) [min⁻¹ per mg/dL]
- `Gb`: Basal glucose concentration [mg/dL]

# Typical values (healthy adult)
- k1 = 0.028 min⁻¹
- k2 = 0.025 min⁻¹
- k3 = 2.5e-5 min⁻² per μU/mL
- k4 = 0.05 min⁻¹
- k5 = 0.015 min⁻¹ per mg/dL
- Gb = 81.0 mg/dL
"""
struct BergmanParams
    k1::Float64  # Glucose effectiveness
    k2::Float64  # Insulin action decay
    k3::Float64  # Insulin sensitivity
    k4::Float64  # Insulin clearance
    k5::Float64  # Pancreatic responsiveness
    Gb::Float64  # Basal glucose
end

"""
    BergmanParams(; kwargs...)

Constructor with keyword arguments and default healthy values.

# Examples
# Healthy adult
params = BergmanParams()

# Type 2 diabetes (insulin resistance)
params_t2d = BergmanParams(k3=1.0e-5, k5=0.010)

# Type 1 diabetes (no insulin secretion)
params_t1d = BergmanParams(k5=0.0)
"""
function BergmanParams(;
    k1::Float64 = 0.028,
    k2::Float64 = 0.025,
    k3::Float64 = 2.5e-5,
    k4::Float64 = 0.05,
    k5::Float64 = 0.015,
    Gb::Float64 = 81.0
)
    return BergmanParams(k1, k2, k3, k4, k5, Gb)
end

"""
    bergman_dynamics!(du, u, p, t)

In-place ODE function for the Bergman minimal model.

# Arguments
- `du`: Derivative vector [dG/dt, dX/dt, dI/dt]
- `u`: State vector [G, X, I]
- `p`: Tuple (params::BergmanParams, D_func, I_ext_func)
- `t`: Current time [min]

# State variables
- `u[1] = G`: Glucose above basal [mg/dL]
- `u[2] = X`: Insulin action [min⁻¹]
- `u[3] = I`: Insulin above basal [μU/mL]

# Model equations
dG/dt = -k1*G - X*(G + Gb) + D(t)
dX/dt = -k2*X + k3*I
dI/dt = -k4*I + k5*max(G - Gb, 0) + I_ext(t)

where:
- D(t) is the glucose input function [mg/dL/min]
- I_ext(t) is the exogenous insulin input [μU/mL/min]
"""
function bergman_dynamics!(du, u, p, t)
    # Unpack parameters and functions
    params, D_func, I_ext_func = p
    k1, k2, k3, k4, k5, Gb = params.k1, params.k2, params.k3, params.k4, params.k5, params.Gb
    
    # Unpack state
    G, X, I = u
    
    # Glucose input (meal/IV)
    D = D_func(t)
    
    # Exogenous insulin (therapy)
    I_ext = I_ext_func(t)
    
    # Glucose dynamics
    # - k1*G: Insulin-independent glucose uptake (brain, RBCs)
    # - X*(G + Gb): Insulin-dependent glucose uptake (muscle, adipose)
    # + D(t): Dietary glucose absorption or IV infusion
    du[1] = -k1 * G - X * (G + Gb) + D
    
    # Insulin action dynamics
    # - k2*X: Decay of insulin action (receptor internalization)
    # + k3*I: Insulin drives glucose uptake (insulin sensitivity)
    du[2] = -k2 * X + k3 * I
    
    # Insulin dynamics
    # - k4*I: Insulin clearance (hepatic/renal degradation)
    # + k5*max(G, 0): Glucose-dependent insulin secretion (β-cells)
    # + I_ext: Exogenous insulin (therapy)
    # Note: max(G - Gb, 0) would be threshold secretion, but since G is above basal,
    # we use max(G, 0) for simplicity when G can be negative
    glucose_stimulus = max(G, 0.0)  # Only secrete if above baseline
    du[3] = -k4 * I + k5 * glucose_stimulus + I_ext
    
    return nothing
end

"""
    bergman_dynamics(u, p, t)

Out-of-place version of the Bergman model dynamics.

Returns du = [dG/dt, dX/dt, dI/dt].
"""
function bergman_dynamics(u, p, t)
    du = similar(u)
    bergman_dynamics!(du, u, p, t)
    return du
end

"""
    insulin_sensitivity(params::BergmanParams)

Calculate the insulin sensitivity index SI.

SI = k3/k2 measures the effectiveness of insulin in lowering glucose.

# Units
(μU/mL)⁻¹ min⁻¹

# Typical values
- Healthy: SI = 5-10 × 10⁻⁴
- Obese: SI = 2-4 × 10⁻⁴
- T2D: SI = 1-2 × 10⁻⁴

# Returns
Insulin sensitivity index
"""
function insulin_sensitivity(params::BergmanParams)
    return params.k3 / params.k2
end

"""
    disposition_index(params::BergmanParams, insulin_response::Float64)

Calculate the disposition index (DI).

DI = SI × Insulin Response measures the adequacy of β-cell compensation
for insulin resistance.

# Arguments
- `params`: Model parameters
- `insulin_response`: Insulin secretion response (e.g., AUC or peak insulin)

# Typical values
- Healthy: DI > 1500
- Prediabetes: DI = 500-1500
- T2D: DI < 500

# Returns
Disposition index
"""
function disposition_index(params::BergmanParams, insulin_response::Float64)
    SI = insulin_sensitivity(params)
    return SI * insulin_response
end

"""
    time_constants(params::BergmanParams)

Calculate characteristic time constants for each compartment.

# Returns
Named tuple with:
- `τG`: Glucose time constant [min]
- `τX`: Insulin action time constant [min]
- `τI`: Insulin time constant [min]
"""
function time_constants(params::BergmanParams)
    τG = 1.0 / params.k1  # Fast: ~35 min
    τX = 1.0 / params.k2  # Slow: ~40 min
    τI = 1.0 / params.k4  # Intermediate: ~20 min
    
    return (τG=τG, τX=τX, τI=τI)
end

"""
    glucose_disposal_rate(G::Float64, X::Float64, params::BergmanParams)

Calculate total glucose disposal rate at given state.

# Components
- Insulin-independent: k1 * G
- Insulin-dependent: X * (G + Gb)

# Arguments
- `G`: Glucose above basal [mg/dL]
- `X`: Insulin action [min⁻¹]
- `params`: Model parameters

# Returns
Total glucose disposal rate [mg/dL/min]
"""
function glucose_disposal_rate(G::Float64, X::Float64, params::BergmanParams)
    insulin_independent = params.k1 * G
    insulin_dependent = X * (G + params.Gb)
    return insulin_independent + insulin_dependent
end

"""
    create_glucose_input(input_type::Symbol, args...)

Create glucose input function D(t).

# Input types
- `:none`: No input (D(t) = 0)
- `:bolus`: Single bolus at t=0
- `:infusion`: Constant infusion
- `:meal`: Gaussian meal absorption
- `:ogtt`: Oral glucose tolerance test (75g over 5 min)
- `:multiple_meals`: Multiple meals throughout day

# Examples
# OGTT (75g glucose over 5 min)
D = create_glucose_input(:ogtt)

# Single meal (60g carbs)
D = create_glucose_input(:meal, t_meal=0.0, carbs=60.0)

# Multiple meals
D = create_glucose_input(:multiple_meals, [(0, 60), (300, 80), (660, 100)])
"""
function create_glucose_input(input_type::Symbol, args...)
    if input_type == :none
        return t -> 0.0
        
    elseif input_type == :bolus
        # Single glucose bolus at t=0
        amount = length(args) > 0 ? args[1] : 300.0  # mg/dL
        return t -> t < 0.1 ? amount * 10 : 0.0  # Concentrated in first 0.1 min
        
    elseif input_type == :infusion
        # Constant glucose infusion
        rate = length(args) > 0 ? args[1] : 50.0  # mg/dL/min
        return t -> rate
        
    elseif input_type == :meal
        # Gaussian meal absorption
        t_meal = length(args) > 0 ? args[1] : 0.0
        carbs = length(args) > 1 ? args[2] : 75.0  # grams
        
        # Convert carbs to glucose appearance rate
        # 1g carb ≈ 4 mg/dL increase, absorbed over ~30 min
        peak_rate = carbs * 4.0 / (30.0 * sqrt(2π))
        σ = 15.0  # Standard deviation (absorption time)
        
        return t -> peak_rate * exp(-(t - t_meal)^2 / (2 * σ^2))
        
    elseif input_type == :ogtt
        # Standard OGTT: 75g glucose over 5 min
        return t -> (0.0 <= t <= 5.0) ? 300.0 : 0.0  # mg/dL/min
        
    elseif input_type == :multiple_meals
        # Multiple meals: args[1] = [(t1, carbs1), (t2, carbs2), ...]
        meals = args[1]
        σ = 15.0
        
        return function(t)
            total = 0.0
            for (t_meal, carbs) in meals
                peak_rate = carbs * 4.0 / (30.0 * sqrt(2π))
                total += peak_rate * exp(-(t - t_meal)^2 / (2 * σ^2))
            end
            return total
        end
        
    else
        error("Unknown input type: $input_type")
    end
end

"""
    create_insulin_input(input_type::Symbol, args...)

Create exogenous insulin input function I_ext(t).

# Input types
- `:none`: No exogenous insulin
- `:basal`: Constant basal insulin
- `:bolus`: Single bolus injection
- `:basal_bolus`: Basal + meal-time boluses

# Examples
# No insulin (healthy or T2D on oral meds)
I_ext = create_insulin_input(:none)

# Basal insulin only (long-acting)
I_ext = create_insulin_input(:basal, 10.0)  # 10 μU/mL

# Basal + bolus
I_ext = create_insulin_input(:basal_bolus, 10.0, [(0.0, 50.0), (300.0, 60.0)])
"""
function create_insulin_input(input_type::Symbol, args...)
    if input_type == :none
        return t -> 0.0
        
    elseif input_type == :basal
        # Constant basal insulin (long-acting)
        basal_level = length(args) > 0 ? args[1] : 10.0  # μU/mL
        return t -> basal_level * 0.05  # Scaled by clearance rate
        
    elseif input_type == :bolus
        # Single bolus at specified time
        t_bolus = length(args) > 0 ? args[1] : 0.0
        amount = length(args) > 1 ? args[2] : 50.0  # μU/mL
        
        # Rapid-acting insulin: peaks at 30-60 min, lasts 3-4 hours
        τ = 60.0  # Time constant
        return t -> (t >= t_bolus) ? (amount / τ) * exp(-(t - t_bolus) / τ) : 0.0
        
    elseif input_type == :basal_bolus
        # Basal + multiple boluses
        basal_level = args[1]
        boluses = args[2]  # [(t1, amount1), (t2, amount2), ...]
        
        basal_func = create_insulin_input(:basal, basal_level)
        
        return function(t)
            total = basal_func(t)
            for (t_bolus, amount) in boluses
                τ = 60.0
                if t >= t_bolus
                    total += (amount / τ) * exp(-(t - t_bolus) / τ)
                end
            end
            return total
        end
        
    else
        error("Unknown insulin input type: $input_type")
    end
end

"""
    steady_state_glucose(D::Float64, params::BergmanParams)

Calculate steady-state glucose level for constant input D.

For constant glucose infusion with no endogenous insulin (T1D),
the steady state is:

G* = D / k1

For healthy individuals with insulin feedback, this is an approximation.

# Arguments
- `D`: Constant glucose input [mg/dL/min]
- `params`: Model parameters

# Returns
Steady-state glucose above basal [mg/dL]
"""
function steady_state_glucose(D::Float64, params::BergmanParams)
    if params.k5 == 0.0
        # Type 1 diabetes: no feedback
        return D / params.k1
    else
        # With feedback: approximate (requires numerical solution)
        @warn "Steady state with feedback requires numerical solution"
        return D / params.k1  # Lower bound
    end
end

export BergmanParams, bergman_dynamics!, bergman_dynamics
export insulin_sensitivity, disposition_index, time_constants
export glucose_disposal_rate, steady_state_glucose
export create_glucose_input, create_insulin_input
