"""
Example 3: Type 2 Diabetes (Insulin Resistance & β-Cell Failure)

Demonstrates the progressive nature of Type 2 diabetes.

Stages simulated:
1. Healthy baseline
2. Early T2D: Insulin resistance with compensation (high insulin)
3. Late T2D: Insulin resistance + β-cell failure

Shows the "march" through parameter space toward metabolic failure.

Expected outcomes:
- Early: Impaired glucose tolerance, hyperinsulinemia
- Late: Hyperglycemia, inadequate insulin
- Progressive deterioration of glucose control
"""

using Plots
using Printf

include("../src/GlucoseInsulinModel.jl")
using .GlucoseInsulinModel

mkpath("../results")

println("=" ^ 70)
println("Example 3: Type 2 Diabetes Progression")
println("=" ^ 70)

# ============================================================================
# 1. Create Parameter Sets for Different Stages
# ============================================================================

println("\n1. Creating parameter sets for disease progression...")

# Stage 0: Healthy
params_healthy = BergmanParams()

# Stage 1: Early T2D (insulin resistance + compensation)
params_early_t2d = BergmanParams(
    k1 = 0.028,
    k2 = 0.025,
    k3 = 1.0e-5,    # 40% of normal (insulin resistance)
    k4 = 0.05,
    k5 = 0.022,     # 150% of normal (compensatory hypersecretion)
    Gb = 81.0
)

# Stage 2: Late T2D (insulin resistance + β-cell failure)
params_late_t2d = BergmanParams(
    k1 = 0.028,
    k2 = 0.025,
    k3 = 0.75e-5,   # 30% of normal (worsening resistance)
    k4 = 0.05,
    k5 = 0.0075,    # 50% of normal (β-cell exhaustion)
    Gb = 81.0
)

println("\n  Parameter comparison:")
println("  " * "─"^60)
println("  Parameter    | Healthy  | Early T2D | Late T2D  |")
println("  " * "─"^60)
println("  k3 (insulin  | 2.50e-5  | 1.00e-5   | 0.75e-5   | (↓ resistance)")
println("  sensitivity) |          | (40%)     | (30%)     |")
println("  " * "─"^60)
println("  k5 (pancreas | 0.015    | 0.022     | 0.0075    | (↑ then ↓)")
println("  response)    |          | (150%)    | (50%)     |")
println("  " * "─"^60)

# Calculate insulin sensitivity indices
SI_healthy = insulin_sensitivity(params_healthy)
SI_early = insulin_sensitivity(params_early_t2d)
SI_late = insulin_sensitivity(params_late_t2d)

println("\n  Insulin Sensitivity Index (SI):")
println("    Healthy:   $(@sprintf("%.2e", SI_healthy)) (μU/mL)⁻¹ min⁻¹")
println("    Early T2D: $(@sprintf("%.2e", SI_early)) ($(@sprintf("%.0f", 100*SI_early/SI_healthy))% of normal)")
println("    Late T2D:  $(@sprintf("%.2e", SI_late)) ($(@sprintf("%.0f", 100*SI_late/SI_healthy))% of normal)")

# ============================================================================
# 2. Simulate All Three Stages
# ============================================================================

println("\n2. Simulating OGTT for all three stages...")

D_func = create_glucose_input(:ogtt)
I_ext = create_insulin_input(:none)

u0 = [0.0, 0.0, 0.0]
tspan = (0.0, 180.0)

# Healthy
println("  • Simulating healthy...")
sol_healthy = simulate_bergman(u0, tspan, params_healthy, D_func, I_ext)
G_healthy = [u[1] + params_healthy.Gb for u in sol_healthy.u]
I_healthy = [u[3] for u in sol_healthy.u]

# Early T2D
println("  • Simulating early T2D...")
sol_early = simulate_bergman(u0, tspan, params_early_t2d, D_func, I_ext)
G_early = [u[1] + params_early_t2d.Gb for u in sol_early.u]
I_early = [u[3] for u in sol_early.u]

# Late T2D
println("  • Simulating late T2D...")
sol_late = simulate_bergman(u0, tspan, params_late_t2d, D_func, I_ext)
G_late = [u[1] + params_late_t2d.Gb for u in sol_late.u]
I_late = [u[3] for u in sol_late.u]

# ============================================================================
# 3. Analyze Results
# ============================================================================

println("\n3. Analyzing glucose tolerance...")

function analyze_ogtt(sol, params, name)
    G_vals = [u[1] + params.Gb for u in sol.u]
    I_vals = [u[3] for u in sol.u]
    
    peak_G = maximum(G_vals)
    peak_G_time = sol.t[argmax(G_vals)]
    
    # 2-hour glucose
    idx_120 = findfirst(t -> t >= 120.0, sol.t)
    G_120 = G_vals[idx_120]
    
    # Diagnosis
    if G_120 < 140.0
        diagnosis = "Normal"
    elseif G_120 < 200.0
        diagnosis = "Prediabetes/IGT"
    else
        diagnosis = "Diabetes"
    end
    
    # Insulin response (AUC approximation)
    insulin_AUC = sum(I_vals) * (sol.t[2] - sol.t[1])
    
    println("\n  $name:")
    println("    Peak glucose:    $(@sprintf("%.1f", peak_G)) mg/dL at $(@sprintf("%.0f", peak_G_time)) min")
    println("    2-hour glucose:  $(@sprintf("%.1f", G_120)) mg/dL")
    println("    Diagnosis:       $diagnosis")
    println("    Insulin AUC:     $(@sprintf("%.0f", insulin_AUC)) (arbitrary units)")
    
    return (peak=peak_G, G_120=G_120, diagnosis=diagnosis, insulin_AUC=insulin_AUC)
end

results_healthy = analyze_ogtt(sol_healthy, params_healthy, "Healthy")
results_early = analyze_ogtt(sol_early, params_early_t2d, "Early T2D")
results_late = analyze_ogtt(sol_late, params_late_t2d, "Late T2D")

# ============================================================================
# 4. Visualizations
# ============================================================================

println("\n4. Creating visualizations...")

# 4.1 Three-way glucose comparison
p1 = plot(sol_healthy.t, G_healthy, lw=3, label="Healthy", color=:blue,
          xlabel="Time (min)", ylabel="Total Glucose (mg/dL)",
          title="Type 2 Diabetes Progression: Glucose Response",
          legend=:topright, size=(1000, 600))
plot!(sol_early.t, G_early, lw=3, label="Early T2D (compensated)", 
      color=:orange, ls=:dash)
plot!(sol_late.t, G_late, lw=3, label="Late T2D (decompensated)", 
      color=:red, ls=:dashdot)
hline!([140], lw=2, ls=:dot, color=:green, label="Normal threshold", alpha=0.5)
hline!([200], lw=2, ls=:dot, color=:red, label="Diabetes threshold", alpha=0.5)
savefig(p1, "../results/03_t2d_progression_glucose.png")
println("    ✓ Saved: results/03_t2d_progression_glucose.png")

# 4.2 Three-way insulin comparison
p2 = plot(sol_healthy.t, I_healthy, lw=3, label="Healthy", color=:blue,
          xlabel="Time (min)", ylabel="Insulin above basal (μU/mL)",
          title="Type 2 Diabetes Progression: Insulin Response",
          legend=:topright, size=(1000, 600))
plot!(sol_early.t, I_early, lw=3, label="Early T2D (hyperinsulinemia)", 
      color=:orange, ls=:dash)
plot!(sol_late.t, I_late, lw=3, label="Late T2D (inadequate)", 
      color=:red, ls=:dashdot)
annotate!(60, maximum(I_early)*0.9, text("Compensatory\nhypersecretion", :orange, 10))
savefig(p2, "../results/03_t2d_progression_insulin.png")
println("    ✓ Saved: results/03_t2d_progression_insulin.png")

# 4.3 E/I ratio over time
EI_ratio_healthy = G_healthy ./ max.(I_healthy, 0.1)
EI_ratio_early = G_early ./ max.(I_early, 0.1)
EI_ratio_late = G_late ./ max.(I_late, 0.1)

p3 = plot(sol_healthy.t, EI_ratio_healthy, lw=3, label="Healthy", color=:blue,
          xlabel="Time (min)", ylabel="Glucose/Insulin Ratio",
          title="Type 2 Diabetes: Glucose/Insulin Ratio",
          legend=:topright, size=(1000, 600))
plot!(sol_early.t, EI_ratio_early, lw=3, label="Early T2D", color=:orange, ls=:dash)
plot!(sol_late.t, EI_ratio_late, lw=3, label="Late T2D", color=:red, ls=:dashdot)
savefig(p3, "../results/03_t2d_ei_ratio.png")
println("    ✓ Saved: results/03_t2d_ei_ratio.png")

# 4.4 Combined view (2x2 grid)
p4 = plot(
    plot(sol_healthy.t, G_healthy, lw=2, label="", color=:blue, 
         title="Healthy", ylabel="Glucose (mg/dL)"),
    plot(sol_healthy.t, I_healthy, lw=2, label="", color=:blue,
         ylabel="Insulin (μU/mL)"),
    plot(sol_early.t, G_early, lw=2, label="", color=:orange,
         title="Early T2D", ylabel="Glucose (mg/dL)"),
    plot(sol_early.t, I_early, lw=2, label="", color=:orange,
         ylabel="Insulin (μU/mL)"),
    plot(sol_late.t, G_late, lw=2, label="", color=:red,
         title="Late T2D", xlabel="Time (min)", ylabel="Glucose (mg/dL)"),
    plot(sol_late.t, I_late, lw=2, label="", color=:red,
         xlabel="Time (min)", ylabel="Insulin (μU/mL)"),
    layout=(3, 2), size=(1200, 900)
)
savefig(p4, "../results/03_t2d_progression_grid.png")
println("    ✓ Saved: results/03_t2d_progression_grid.png")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 70)
println("SUMMARY: Type 2 Diabetes Progression")
println("=" ^ 70)
println("\nStage 1: Healthy")
println("  • Normal insulin sensitivity (SI = $(@sprintf("%.2e", SI_healthy)))")
println("  • 2-hour glucose: $(@sprintf("%.0f", results_healthy.G_120)) mg/dL")
println("  • Diagnosis: $(results_healthy.diagnosis)")
println("\nStage 2: Early T2D (Compensated)")
println("  • Insulin resistance (SI = $(@sprintf("%.2e", SI_early)), $(@sprintf("%.0f", 100*SI_early/SI_healthy))%)")
println("  • Compensatory hyperinsulinemia (insulin ↑ $(@sprintf("%.0f", 100*results_early.insulin_AUC/results_healthy.insulin_AUC))%)")
println("  • 2-hour glucose: $(@sprintf("%.0f", results_early.G_120)) mg/dL")
println("  • Diagnosis: $(results_early.diagnosis)")
println("  • Key: β-cells working overtime to compensate")
println("\nStage 3: Late T2D (Decompensated)")
println("  • Severe insulin resistance (SI = $(@sprintf("%.2e", SI_late)), $(@sprintf("%.0f", 100*SI_late/SI_healthy))%)")
println("  • β-cell failure (insulin response inadequate)")
println("  • 2-hour glucose: $(@sprintf("%.0f", results_late.G_120)) mg/dL")
println("  • Diagnosis: $(results_late.diagnosis)")
println("  • Key: β-cells exhausted, cannot compensate")
println("\nControl Theory Perspective:")
println("  • Early T2D: ACTUATOR FAILURE (insulin resistance)")
println("    → System compensates by increasing GAIN (k5 ↑)")
println("    → Maintains control, but at high metabolic cost")
println("  • Late T2D: ACTUATOR + SENSOR FAILURE")
println("    → Actuator broken (k3 ↓) AND gain saturated (k5 ↓)")
println("    → Control lost, open-loop dynamics")
println("\nPathophysiology Insight:")
println("  T2D is a PROGRESSIVE disease:")
println("  1. Insulin resistance develops (obesity, genetics, aging)")
println("  2. Pancreas compensates by overproducing insulin")
println("  3. Over years/decades, β-cells burn out")
println("  4. Insulin production falls despite high demand")
println("  5. Result: Hyperglycemia + relative insulin deficiency")
println("\nClinical Implications:")
println("  • EARLY INTERVENTION is critical")
println("  • Lifestyle changes can preserve β-cell function")
println("  • Once β-cells fail, insulin may be required")
println("  • Preventing progression is KEY: \"β-cells don't regenerate\"")
println("=" ^ 70)
