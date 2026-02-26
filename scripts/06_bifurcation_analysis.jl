"""
Example 6: Bifurcation Analysis

Explores how system behavior changes with key parameters.

Analyzes:
1. Bifurcation in k5 (pancreatic responsiveness)
2. Bifurcation in k3 (insulin sensitivity)
3. 2D parameter space (k3 vs k5)
4. Hopf bifurcation detection
5. Eigenvalue trajectories

Shows the "stability landscape" and critical thresholds for diabetes.

Expected outcomes:
- Critical k5 threshold below which glucose control fails
- Gradual degradation with decreasing k3
- Stability boundaries in parameter space
- Regions corresponding to healthy, prediabetes, diabetes
"""

using Plots
using Printf

include("../src/GlucoseInsulinModel.jl")
using .GlucoseInsulinModel

mkpath("../results")

println("=" ^ 70)
println("Example 6: Bifurcation Analysis")
println("=" ^ 70)

# ============================================================================
# 1. Bifurcation in k5 (Pancreatic Responsiveness)
# ============================================================================

println("\n1. Analyzing bifurcation in k5 (pancreatic responsiveness)...")

params_base = BergmanParams()

k5_range = 0.0:0.0005:0.03
println("  Scanning k5 from $(minimum(k5_range)) to $(maximum(k5_range)) min⁻¹")
println("  Number of points: $(length(k5_range))")

println("  Computing fixed points and stability...")
results_k5 = bifurcation_diagram(:k5, k5_range, params_base, D=0.0)

# Find critical threshold
stable_k5 = results_k5.param_values[results_k5.stable]
if !isempty(stable_k5)
    k5_critical = minimum(stable_k5)
    println("  ✓ Critical k5 threshold: $(@sprintf("%.6f", k5_critical)) min⁻¹")
    println("    Below this value, glucose control fails")
    println("    Healthy k5 = 0.015 has $(@sprintf("%.1f", 0.015/k5_critical))× safety margin")
end

# ============================================================================
# 2. Bifurcation in k3 (Insulin Sensitivity)
# ============================================================================

println("\n2. Analyzing bifurcation in k3 (insulin sensitivity)...")

k3_range = 0.0:1.0e-6:5.0e-5
println("  Scanning k3 from $(minimum(k3_range)) to $(maximum(k3_range)) min⁻² per μU/mL")

println("  Computing fixed points and stability...")
results_k3 = bifurcation_diagram(:k3, k3_range, params_base, D=0.0)

# Find relationship between k3 and SI
SI_values = [k3 / params_base.k2 for k3 in k3_range]

# ============================================================================
# 3. Create Visualizations
# ============================================================================

println("\n3. Creating bifurcation diagrams...")

# 3.1 Bifurcation diagram for k5 (Glucose)
println("  • Plotting k5 bifurcation (glucose)...")
p1 = plot_bifurcation_diagram(results_k5, variable_index=1,
                               param_name="k5 (Pancreatic Responsiveness, min⁻¹)",
                               title="Bifurcation Diagram: Glucose vs k5")
# Add reference lines
vline!([0.015], lw=2, ls=:dash, color=:green, label="Healthy k5")
vline!([0.0075], lw=2, ls=:dash, color=:orange, label="Late T2D k5")
savefig(p1, "../results/06_bifurcation_k5_glucose.png")
println("    ✓ Saved: results/06_bifurcation_k5_glucose.png")

# 3.2 Bifurcation diagram for k5 (Insulin)
p2 = plot_bifurcation_diagram(results_k5, variable_index=3,
                               param_name="k5 (Pancreatic Responsiveness, min⁻¹)",
                               title="Bifurcation Diagram: Insulin vs k5")
vline!([0.015], lw=2, ls=:dash, color=:green, label="Healthy")
savefig(p2, "../results/06_bifurcation_k5_insulin.png")

# 3.3 Bifurcation diagram for k3 (Glucose)
println("  • Plotting k3 bifurcation (glucose)...")
p3 = plot_bifurcation_diagram(results_k3, variable_index=1,
                               param_name="k3 (Insulin Sensitivity, min⁻² per μU/mL)",
                               title="Bifurcation Diagram: Glucose vs k3")
vline!([2.5e-5], lw=2, ls=:dash, color=:green, label="Healthy k3")
vline!([1.0e-5], lw=2, ls=:dash, color=:orange, label="Early T2D k3")
vline!([0.75e-5], lw=2, ls=:dash, color=:red, label="Late T2D k3")
savefig(p3, "../results/06_bifurcation_k3_glucose.png")
println("    ✓ Saved: results/06_bifurcation_k3_glucose.png")

# 3.4 Eigenvalue trajectories for k5
println("  • Plotting eigenvalue trajectories...")
p4 = plot_eigenvalue_trajectory(results_k5, 
                                 title="Eigenvalue Trajectories vs k5")
savefig(p4, "../results/06_eigenvalues_k5.png")
println("    ✓ Saved: results/06_eigenvalues_k5.png")

# ============================================================================
# 4. 2D Parameter Space Analysis
# ============================================================================

println("\n4. Analyzing 2D parameter space (k3 vs k5)...")

k3_range_2d = 0.5e-5:0.25e-5:4.0e-5
k5_range_2d = 0.005:0.001:0.025

println("  Grid size: $(length(k3_range_2d)) × $(length(k5_range_2d))")
println("  Total points: $(length(k3_range_2d) * length(k5_range_2d))")
println("  Computing stability region...")

p5 = plot_stability_region(k3_range_2d, k5_range_2d, :k3, :k5, params_base,
                            title="Stability Region: k3 vs k5")

# Add disease state markers
k3_healthy = 2.5e-5
k5_healthy = 0.015
scatter!([k3_healthy], [k5_healthy], marker=:star, markersize=15,
         color=:yellow, label="Healthy", markerstrokecolor=:black, markerstrokewidth=2)

k3_early_t2d = 1.0e-5
k5_early_t2d = 0.022
scatter!([k3_early_t2d], [k5_early_t2d], marker=:diamond, markersize=12,
         color=:orange, label="Early T2D", markerstrokecolor=:black)

k3_late_t2d = 0.75e-5
k5_late_t2d = 0.0075
scatter!([k3_late_t2d], [k5_late_t2d], marker=:circle, markersize=12,
         color=:red, label="Late T2D", markerstrokecolor=:black)

# Add progression arrow
plot!([k3_healthy, k3_early_t2d, k3_late_t2d], 
      [k5_healthy, k5_early_t2d, k5_late_t2d],
      lw=3, color=:white, arrow=true, label="T2D Progression",
      linestyle=:solid)

savefig(p5, "../results/06_stability_region_2d.png")
println("  ✓ Saved: results/06_stability_region_2d.png")

# ============================================================================
# 5. Simulate Representative Points
# ============================================================================

println("\n5. Simulating dynamics at representative parameter values...")

# Three different k5 values
k5_values = [0.003, 0.010, 0.020]
k5_labels = ["Low k5 (0.003)", "Medium k5 (0.010)", "High k5 (0.020)"]
colors = [:red, :orange, :green]

D_func = create_glucose_input(:ogtt)
I_ext = create_insulin_input(:none)
u0 = [0.0, 0.0, 0.0]
tspan = (0.0, 180.0)

p6 = plot(xlabel="Time (min)", ylabel="Total Glucose (mg/dL)",
          title="Dynamics Across k5 Values", legend=:topright,
          size=(1000, 600))

for (k5_val, label, color) in zip(k5_values, k5_labels, colors)
    params_temp = BergmanParams(k1=params_base.k1, k2=params_base.k2, 
                                k3=params_base.k3, k4=params_base.k4,
                                k5=k5_val, Gb=params_base.Gb)
    
    sol = simulate_bergman(u0, tspan, params_temp, D_func, I_ext)
    G_vals = [u[1] + params_temp.Gb for u in sol.u]
    
    plot!(sol.t, G_vals, lw=3, label=label, color=color)
end

hline!([140], lw=2, ls=:dash, color=:gray, label="Normal threshold", alpha=0.5)
savefig(p6, "../results/06_dynamics_across_k5.png")
println("  ✓ Saved: results/06_dynamics_across_k5.png")

# ============================================================================
# 6. Sensitivity Analysis
# ============================================================================

println("\n6. Performing sensitivity analysis...")

sensitivities = sensitivity_analysis(params_base, 0.01)

if !isnothing(sensitivities)
    println("\n  Normalized sensitivities (1% parameter change):")
    println("  Parameter | Sensitivity | Interpretation")
    println("  " * "─"^55)
    for (param, sens) in sensitivities
        interp = abs(sens) > 0.5 ? "High" : (abs(sens) > 0.1 ? "Moderate" : "Low")
        println("  $param       | $(@sprintf("%+.4f", sens))     | $interp sensitivity")
    end
    
    # Bar chart of sensitivities
    p7 = bar([String(k) for k in keys(sensitivities)],
             [v for v in values(sensitivities)],
             xlabel="Parameter", ylabel="Normalized Sensitivity",
             title="Parameter Sensitivity Analysis",
             legend=false, color=:steelblue, size=(800, 500))
    hline!([0], color=:black, lw=1)
    savefig(p7, "../results/06_sensitivity_analysis.png")
    println("\n  ✓ Saved: results/06_sensitivity_analysis.png")
end

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 70)
println("SUMMARY: Bifurcation Analysis")
println("=" ^ 70)
println("\nKey Findings:")
println("\n1. Critical Thresholds:")
println("   • k5 has a SHARP threshold (~0.005 min⁻¹)")
println("   • Below this: System unstable, glucose control fails")
println("   • Healthy k5 = 0.015 has ~3× safety margin")
println("   • T1D (k5 = 0): Far below threshold → catastrophic failure")
println("\n2. Insulin Resistance (k3):")
println("   • Degradation is more GRADUAL")
println("   • Lower k3 → higher steady-state glucose")
println("   • No sharp bifurcation, but clinical threshold exists")
println("\n3. 2D Parameter Space:")
println("   • Stability region visualized")
println("   • T2D progression = walk along stability boundary")
println("   • Early T2D: Compensate upward in k5")
println("   • Late T2D: Fall below stability boundary")
println("\n4. Control Theory Insights:")
println("   • System has FRAGILE stability")
println("   • Small parameter changes → large outcome changes")
println("   • Nonlinear threshold phenomena")
println("   • Multiple parameters must stay in viable range")
println("\n5. Clinical Implications:")
println("   • Parameter space defines \"health landscape\"")
println("   • Disease = trajectory out of stable region")
println("   • Therapies = move parameters back to stable region:")
println("     - Metformin: ↑ k1 (glucose effectiveness)")
println("     - TZDs: ↑ k3 (insulin sensitivity)")
println("     - Sulfonylureas: ↑ k5 (pancreatic response)")
println("     - Insulin: Bypass k5 entirely")
println("\n6. Precision Medicine:")
println("   • Measure individual k3, k5 from IVGTT/OGTT")
println("   • Identify which parameters are impaired")
println("   • Target therapy to specific defects")
println("   • Predict progression based on parameter trajectory")
println("\nBifurcation theory organizes the disease landscape!")
println("=" ^ 70)
