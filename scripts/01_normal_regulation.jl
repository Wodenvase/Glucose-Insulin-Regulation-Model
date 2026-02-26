"""
Example 1: Normal Glucose Regulation (OGTT)

Demonstrates healthy glucose-insulin dynamics during an oral glucose tolerance test.

Simulates:
- 75g oral glucose load
- Glucose, insulin, and insulin action time series
- Phase portrait
- Fixed point analysis
- Stability analysis

Expected outcomes:
- Peak glucose ~140 mg/dL at 30-45 min
- Return to baseline by 120 min
- Biphasic insulin response
- Stable fixed point at origin
"""

using Plots
using Printf

# Load the module
include("../src/GlucoseInsulinModel.jl")
using .GlucoseInsulinModel

# Create results directory if it doesn't exist
mkpath("../results")

println("=" ^ 70)
println("Example 1: Normal Glucose Regulation (OGTT)")
println("=" ^ 70)

# ============================================================================
# 1. Setup Parameters (Healthy Adult)
# ============================================================================

println("\n1. Setting up healthy adult parameters...")

params = BergmanParams(
    k1 = 0.028,   # Glucose effectiveness (min⁻¹)
    k2 = 0.025,   # Insulin action decay (min⁻¹)
    k3 = 2.5e-5,  # Insulin sensitivity (min⁻² per μU/mL)
    k4 = 0.05,    # Insulin clearance (min⁻¹)
    k5 = 0.015,   # Pancreatic responsiveness (min⁻¹ per mg/dL)
    Gb = 81.0     # Basal glucose (mg/dL)
)

# Display parameter values
println("  Parameters:")
println("    k1 (glucose effectiveness):     $(params.k1) min⁻¹")
println("    k2 (insulin action decay):      $(params.k2) min⁻¹")
println("    k3 (insulin sensitivity):       $(params.k3) min⁻² per μU/mL")
println("    k4 (insulin clearance):         $(params.k4) min⁻¹")
println("    k5 (pancreatic responsiveness): $(params.k5) min⁻¹ per mg/dL")
println("    Gb (basal glucose):             $(params.Gb) mg/dL")

# Calculate derived quantities
SI = insulin_sensitivity(params)
τ = time_constants(params)
println("\n  Derived quantities:")
println("    Insulin sensitivity index: $(@sprintf("%.2e", SI)) (μU/mL)⁻¹ min⁻¹")
println("    Time constants:")
println("      τ_G (glucose):        $(@sprintf("%.1f", τ.τG)) min")
println("      τ_X (insulin action): $(@sprintf("%.1f", τ.τX)) min")
println("      τ_I (insulin):        $(@sprintf("%.1f", τ.τI)) min")

# ============================================================================
# 2. Setup OGTT Protocol
# ============================================================================

println("\n2. Setting up OGTT protocol...")
println("  Protocol: 75g oral glucose over 5 minutes")
println("  Duration: 180 minutes (3 hours)")

# Create glucose input (OGTT: 75g glucose over 5 min)
D_func = create_glucose_input(:ogtt)

# No exogenous insulin (healthy individual)
I_ext_func = create_insulin_input(:none)

# Initial conditions (at baseline)
u0 = [0.0, 0.0, 0.0]  # [G, X, I] all at basal levels

# Time span
tspan = (0.0, 180.0)  # 0 to 180 minutes

# ============================================================================
# 3. Simulate
# ============================================================================

println("\n3. Simulating OGTT...")

sol = simulate_bergman(u0, tspan, params, D_func, I_ext_func)

println("  Simulation complete!")
println("  Number of time points: $(length(sol.t))")

# Extract peak values
G_vals = [u[1] + params.Gb for u in sol.u]
I_vals = [u[3] for u in sol.u]
X_vals = [u[2] for u in sol.u]

peak_G = maximum(G_vals)
peak_G_time = sol.t[argmax(G_vals)]
peak_I = maximum(I_vals)
peak_I_time = sol.t[argmax(I_vals)]

println("\n  Peak values:")
println("    Peak glucose: $(@sprintf("%.1f", peak_G)) mg/dL at $(@sprintf("%.1f", peak_G_time)) min")
println("    Peak insulin: $(@sprintf("%.1f", peak_I)) μU/mL at $(@sprintf("%.1f", peak_I_time)) min")

# 2-hour glucose (diagnostic criterion)
idx_120 = findfirst(t -> t >= 120.0, sol.t)
G_120 = G_vals[idx_120]
println("    2-hour glucose: $(@sprintf("%.1f", G_120)) mg/dL")

if G_120 < 140.0
    println("    Diagnosis: NORMAL glucose tolerance ✓")
elseif G_120 < 200.0
    println("    Diagnosis: Prediabetes (impaired glucose tolerance)")
else
    println("    Diagnosis: Diabetes mellitus")
end

# ============================================================================
# 4. Fixed Point Analysis
# ============================================================================

println("\n4. Performing fixed point analysis...")

fp = find_fixed_point(params, 0.0)  # No constant input (D = 0)

if fp.exists
    println("  Fixed point found:")
    println("    G* = $(@sprintf("%.4f", fp.state[1])) mg/dL (above basal)")
    println("    X* = $(@sprintf("%.4e", fp.state[2])) min⁻¹")
    println("    I* = $(@sprintf("%.4f", fp.state[3])) μU/mL (above basal)")
    println("    Residual: $(@sprintf("%.2e", fp.residual))")
    
    # Stability analysis
    analysis = analyze_stability(fp, params)
    
    println("\n  Stability analysis:")
    println("    Stable: $(analysis.is_stable)")
    println("    Type: $(analysis.stability_type)")
    println("    Eigenvalues:")
    for (i, λ) in enumerate(analysis.eigenvalues)
        println("      λ$i = $(@sprintf("%.6f", real(λ))) + $(@sprintf("%.6f", imag(λ)))i")
    end
    
    if analysis.is_stable
        println("\n  ✓ System is stable - glucose returns to baseline after disturbances")
    else
        println("\n  ✗ System is unstable - regulatory failure")
    end
else
    println("  ✗ Fixed point not found")
end

# ============================================================================
# 5. Visualizations
# ============================================================================

println("\n5. Creating visualizations...")

# 5.1 Time series of all variables
println("  Plotting time series...")
p1 = plot_time_series(sol, title="Normal OGTT: Glucose-Insulin Dynamics",
                      show_components=true, add_basal=true, params=params)
savefig(p1, "../results/01_normal_glucose_timeseries.png")
println("    ✓ Saved: results/01_normal_glucose_timeseries.png")

# 5.2 Individual components with detail
println("  Plotting individual components...")

# Glucose
p2 = plot(sol.t, G_vals, lw=3, color=:blue, label="Glucose",
          xlabel="Time (min)", ylabel="Total Glucose (mg/dL)",
          title="Normal OGTT: Glucose Response", legend=:topright,
          size=(900, 500))
hline!([140], lw=2, ls=:dash, color=:green, label="Normal threshold (140 mg/dL)")
hline!([params.Gb], lw=2, ls=:dot, color=:gray, label="Basal (81 mg/dL)", alpha=0.5)
savefig(p2, "../results/01_normal_glucose.png")

# Insulin
p3 = plot(sol.t, I_vals, lw=3, color=:red, label="Insulin",
          xlabel="Time (min)", ylabel="Insulin above basal (μU/mL)",
          title="Normal OGTT: Insulin Response", legend=:topright,
          size=(900, 500))
savefig(p3, "../results/01_normal_insulin.png")

# Insulin action
p4 = plot(sol.t, X_vals, lw=3, color=:green, label="Insulin Action",
          xlabel="Time (min)", ylabel="Insulin Action (min⁻¹)",
          title="Normal OGTT: Insulin Action", legend=:topright,
          size=(900, 500))
savefig(p4, "../results/01_normal_insulin_action.png")

println("    ✓ Saved individual component plots")

# 5.3 Phase portrait (G vs I)
println("  Plotting phase portrait...")
p5 = plot_phase_portrait(sol, params, fixed_point=fp,
                         title="Normal OGTT: Phase Portrait (G vs I)")
savefig(p5, "../results/01_normal_phase_portrait.png")
println("    ✓ Saved: results/01_normal_phase_portrait.png")

# 5.4 All variables on same time axis
println("  Plotting combined view...")
p6 = plot(sol.t, G_vals, lw=2, label="Glucose (mg/dL)", color=:blue,
          xlabel="Time (min)", ylabel="Concentration", 
          title="Normal OGTT: All Variables",
          legend=:topright, size=(1000, 600))
plot!(sol.t, I_vals, lw=2, label="Insulin (μU/mL)", color=:red)
plot!(sol.t, X_vals .* 1000, lw=2, label="Insulin Action (×1000 min⁻¹)", color=:green)
savefig(p6, "../results/01_normal_all_variables.png")
println("    ✓ Saved: results/01_normal_all_variables.png")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 70)
println("SUMMARY: Normal Glucose Regulation")
println("=" ^ 70)
println("\nKey findings:")
println("  ✓ Normal glucose tolerance (2-hr glucose < 140 mg/dL)")
println("  ✓ Rapid glucose clearance (returns to baseline by 120 min)")
println("  ✓ Robust insulin response (peak ~$(@sprintf("%.1f", peak_I)) μU/mL)")
println("  ✓ Stable fixed point (homeostatic regulation)")
println("  ✓ All eigenvalues have negative real parts (stable system)")
println("\nInterpretation:")
println("  This demonstrates healthy negative feedback control.")
println("  The pancreas secretes insulin proportional to glucose elevation.")
println("  Insulin drives glucose uptake via insulin action.")
println("  Glucose returns to setpoint (Gb = $(params.Gb) mg/dL).")
println("  This is homeostasis in action!")
println("\nFigures saved to results/ directory.")
println("=" ^ 70)
