# Quick Start Guide

Get started with the Glucose-Insulin Regulation Model in 5 minutes.

## Installation

1. **Clone or navigate to the project**:
   ```bash
   cd /Users/dipantabhattacharyya/Insulin
   ```

2. **Start Julia with the project**:
   ```bash
   julia --project=.
   ```

3. **Install dependencies** (first time only):
   ```julia
   using Pkg
   Pkg.instantiate()
   ```

   This will install:
   - DifferentialEquations.jl (ODE solvers)
   - ForwardDiff.jl (automatic differentiation)
   - NLsolve.jl (nonlinear equation solving)
   - Plots.jl (visualization)
   - LinearAlgebra, Statistics (standard library)

## Basic Usage

### Example 1: Simulate Normal OGTT

```julia
using Plots
include("src/GlucoseInsulinModel.jl")
using .GlucoseInsulinModel

# Create healthy parameters
params = BergmanParams()

# OGTT protocol: 75g glucose over 5 minutes
D = create_glucose_input(:ogtt)
I_ext = create_insulin_input(:none)

# Simulate
u0 = [0.0, 0.0, 0.0]  # Start at baseline
tspan = (0.0, 180.0)   # 3 hours
sol = simulate_bergman(u0, tspan, params, D, I_ext)

# Visualize
p = plot_time_series(sol, params=params, add_basal=true)
display(p)
```

### Example 2: Type 1 Diabetes

```julia
# Create T1D parameters (no insulin secretion)
params_t1d = BergmanParams(k5=0.0)

# Simulate without treatment
sol_t1d = simulate_bergman(u0, tspan, params_t1d, D, I_ext)

# Compare with normal
plot_comparison([sol, sol_t1d], ["Normal", "Type 1 Diabetes"], 
                params=params, add_basal=true)
```

### Example 3: Insulin Therapy

```julia
# Add insulin therapy
basal = 10.0  # μU/mL
bolus = 50.0  # μU/mL at meal
I_therapy = create_insulin_input(:basal_bolus, basal, [(0.0, bolus)])

# Simulate T1D with treatment
sol_treated = simulate_bergman(u0, tspan, params_t1d, D, I_therapy)

# Visualize
plot_comparison([sol_t1d, sol_treated], 
                ["T1D (no treatment)", "T1D (with insulin)"],
                params=params_t1d, add_basal=true)
```

### Example 4: Bifurcation Analysis

```julia
# Analyze how k5 affects stability
k5_range = 0.0:0.001:0.03
results = bifurcation_diagram(:k5, k5_range, params)

# Plot bifurcation diagram
plot_bifurcation_diagram(results, variable_index=1,
                         param_name="k5 (Pancreatic Responsiveness)",
                         title="Bifurcation Diagram: Glucose vs k5")
```

## Running Complete Examples

Run the comprehensive example scripts:

```julia
# Example 1: Normal regulation
include("scripts/01_normal_regulation.jl")

# Example 2: Type 1 diabetes
include("scripts/02_type1_diabetes.jl")

# Example 3: Type 2 diabetes progression
include("scripts/03_type2_diabetes.jl")

# Example 6: Bifurcation analysis
include("scripts/06_bifurcation_analysis.jl")

# Example 7: Drug interventions
include("scripts/07_drug_interventions.jl")
```

Each script will:
- Print detailed analysis to console
- Save figures to `results/` directory
- Provide biological/clinical interpretations

## Key Functions

### Parameter Creation
```julia
# Healthy adult
params = BergmanParams()

# Custom parameters
params = BergmanParams(k1=0.028, k2=0.025, k3=2.5e-5, 
                       k4=0.05, k5=0.015, Gb=81.0)

# Disease states
params_t1d = create_disease_params(:type1)
params_t2d = create_disease_params(:type2_late, severity=:severe)
```

### Glucose Inputs
```julia
D = create_glucose_input(:none)           # No input
D = create_glucose_input(:ogtt)           # Standard OGTT
D = create_glucose_input(:meal, 0, 60)    # 60g meal at t=0
D = create_glucose_input(:multiple_meals, 
                         [(0, 60), (300, 80), (660, 100)])
```

### Insulin Inputs
```julia
I_ext = create_insulin_input(:none)                    # No insulin
I_ext = create_insulin_input(:basal, 10.0)            # Basal only
I_ext = create_insulin_input(:basal_bolus, 10.0, 
                             [(0, 50), (300, 60)])     # Basal + boluses
```

### Simulation
```julia
sol = simulate_bergman(u0, tspan, params, D_func, I_ext_func)

# With event detection
sol, events = simulate_with_callback(u0, tspan, params, D_func, I_ext_func,
                                     hypoglycemia_threshold=70.0,
                                     hyperglycemia_threshold=250.0)
```

### Analysis
```julia
# Fixed point analysis
fp = find_fixed_point(params, D=0.0)

# Stability analysis
analysis = analyze_stability(fp, params)
println("Stable: ", analysis.is_stable)
println("Eigenvalues: ", analysis.eigenvalues)

# Insulin sensitivity
SI = insulin_sensitivity(params)
```

### Visualization
```julia
# Time series
plot_time_series(sol, params=params, add_basal=true)

# Phase portrait
plot_phase_portrait(sol, params, fixed_point=fp)

# Bifurcation diagram
plot_bifurcation_diagram(results, variable_index=1)

# Comparison
plot_comparison(sols, labels, params=params)
```

## Parameter Reference

| Parameter | Meaning | Healthy Value | Units |
|-----------|---------|---------------|-------|
| k1 | Glucose effectiveness | 0.028 | min⁻¹ |
| k2 | Insulin action decay | 0.025 | min⁻¹ |
| k3 | Insulin sensitivity | 2.5×10⁻⁵ | min⁻² per μU/mL |
| k4 | Insulin clearance | 0.05 | min⁻¹ |
| k5 | Pancreatic responsiveness | 0.015 | min⁻¹ per mg/dL |
| Gb | Basal glucose | 81.0 | mg/dL |

## State Variables

| Variable | Meaning | Units |
|----------|---------|-------|
| G | Glucose above basal | mg/dL |
| X | Insulin action | min⁻¹ |
| I | Insulin above basal | μU/mL |

Total glucose = G + Gb

## Disease Parameters

### Type 1 Diabetes
- k5 = 0.0 (no insulin secretion)

### Type 2 Diabetes (Early)
- k3 ≈ 40% of normal (insulin resistance)
- k5 ≈ 150% of normal (compensation)

### Type 2 Diabetes (Late)
- k3 ≈ 30% of normal (severe resistance)
- k5 ≈ 50% of normal (β-cell failure)

## Troubleshooting

### Package Installation Issues
```julia
# Update package registry
using Pkg
Pkg.update()
Pkg.instantiate()
```

### Plotting Issues
```julia
# If plots don't display
using Plots
gr()  # Use GR backend
```

### Slow Simulations
```julia
# Use faster ODE solver for stiff systems
sol = simulate_bergman(u0, tspan, params, D, I_ext, solver=Rodas5())
```

## Next Steps

1. **Read the full README.md** for comprehensive background
2. **Run example scripts** to see complete analyses
3. **Explore parameter space** with bifurcation analysis
4. **Create your own experiments** modifying parameters
5. **Extend the model** (see README Extensions section)

## Getting Help

- Check the comprehensive docstrings: `?BergmanParams`, `?simulate_bergman`, etc.
- Read the example scripts for usage patterns
- See README.md for theoretical background

## Citation

If you use this model in research, please cite:
- Bergman, R. N., et al. (1979). *Journal of Clinical Investigation*, 68(6), 1456-1467.
- Bergman, R. N., et al. (1981). *Endocrine Reviews*, 6(1), 45-86.

---

**Happy modeling!** 🩺📊
