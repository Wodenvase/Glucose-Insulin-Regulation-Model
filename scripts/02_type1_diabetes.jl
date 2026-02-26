"""
Example 2: Type 1 Diabetes (Insulin Deficiency)

Demonstrates complete insulin deficiency due to β-cell destruction.

Simulates:
- Same OGTT protocol as Example 1
- No endogenous insulin secretion (k5 = 0)
- Comparison with normal response
- Insulin therapy intervention

Expected outcomes:
- Severe hyperglycemia (>300 mg/dL)
- No insulin response
- No glucose clearance
- System failure without treatment
"""

using Plots
using Printf

include("../src/GlucoseInsulinModel.jl")
using .GlucoseInsulinModel

mkpath("../results")

println("=" ^ 70)
println("Example 2: Type 1 Diabetes (Insulin Deficiency)")
println("=" ^ 70)

# ============================================================================
# 1. Setup Type 1 Diabetes Parameters
# ============================================================================

println("\n1. Creating Type 1 diabetes parameters...")

params_t1d = BergmanParams(
    k1 = 0.028,
    k2 = 0.025,
    k3 = 2.5e-5,
    k4 = 0.05,
    k5 = 0.0,    # NO insulin secretion (β-cell destruction)
    Gb = 81.0
)

# Normal parameters for comparison
params_normal = BergmanParams()

println("  Key parameter change:")
println("    k5 (pancreatic responsiveness): 0.0 min⁻¹ (was 0.015)")
println("    → Complete absence of insulin secretion")

# ============================================================================
# 2. Simulate T1D Without Treatment
# ============================================================================

println("\n2. Simulating Type 1 diabetes WITHOUT insulin therapy...")

D_func = create_glucose_input(:ogtt)
I_ext_none = create_insulin_input(:none)

u0 = [0.0, 0.0, 0.0]
tspan = (0.0, 180.0)

sol_t1d = simulate_bergman(u0, tspan, params_t1d, D_func, I_ext_none)

G_t1d = [u[1] + params_t1d.Gb for u in sol_t1d.u]
I_t1d = [u[3] for u in sol_t1d.u]

peak_G_t1d = maximum(G_t1d)
final_G_t1d = G_t1d[end]

println("  Results:")
println("    Peak glucose: $(@sprintf("%.1f", peak_G_t1d)) mg/dL")
println("    Final glucose (180 min): $(@sprintf("%.1f", final_G_t1d)) mg/dL")
println("    ✗ SEVERE HYPERGLYCEMIA - Life-threatening!")

# ============================================================================
# 3. Simulate Normal for Comparison
# ============================================================================

println("\n3. Simulating normal response for comparison...")

sol_normal = simulate_bergman(u0, tspan, params_normal, D_func, I_ext_none)
G_normal = [u[1] + params_normal.Gb for u in sol_normal.u]

# ============================================================================
# 4. Simulate T1D WITH Insulin Therapy
# ============================================================================

println("\n4. Simulating Type 1 diabetes WITH insulin therapy...")

# Basal-bolus insulin regimen
basal_insulin = 10.0  # μU/mL constant background
bolus_insulin = 50.0  # μU/mL at meal time

I_ext_therapy = create_insulin_input(:basal_bolus, basal_insulin, [(0.0, bolus_insulin)])

sol_t1d_therapy = simulate_bergman(u0, tspan, params_t1d, D_func, I_ext_therapy)
G_t1d_therapy = [u[1] + params_t1d.Gb for u in sol_t1d_therapy.u]
I_t1d_therapy = [u[3] for u in sol_t1d_therapy.u]

peak_G_therapy = maximum(G_t1d_therapy)
final_G_therapy = G_t1d_therapy[end]

println("  Insulin regimen:")
println("    Basal: $(basal_insulin) μU/mL (continuous)")
println("    Bolus: $(bolus_insulin) μU/mL (at meal)")
println("\n  Results with therapy:")
println("    Peak glucose: $(@sprintf("%.1f", peak_G_therapy)) mg/dL")
println("    Final glucose: $(@sprintf("%.1f", final_G_therapy)) mg/dL")
println("    ✓ Much improved control!")

# ============================================================================
# 5. Visualizations
# ============================================================================

println("\n5. Creating visualizations...")

# 5.1 Comparison: Normal vs T1D
p1 = plot(sol_normal.t, G_normal, lw=3, label="Normal", color=:blue,
          xlabel="Time (min)", ylabel="Total Glucose (mg/dL)",
          title="Comparison: Normal vs Type 1 Diabetes (Untreated)",
          legend=:topleft, size=(1000, 600))
plot!(sol_t1d.t, G_t1d, lw=3, label="Type 1 Diabetes (no insulin)", 
      color=:red, ls=:dash)
hline!([140], lw=2, ls=:dot, color=:green, label="Normal threshold", alpha=0.5)
hline!([250], lw=2, ls=:dot, color=:red, label="Severe hyperglycemia", alpha=0.5)
savefig(p1, "../results/02_comparison_normal_t1d.png")
println("    ✓ Saved: results/02_comparison_normal_t1d.png")

# 5.2 T1D: Glucose, Insulin, and Insulin Action
p2 = plot(sol_t1d.t, G_t1d, lw=3, label="Glucose", color=:blue,
          xlabel="", ylabel="Glucose (mg/dL)",
          title="Type 1 Diabetes: Complete Regulatory Failure",
          legend=:topright, size=(1000, 800), layout=(3,1))

plot!(subplot=2, sol_t1d.t, I_t1d, lw=3, label="Insulin", color=:red,
      xlabel="", ylabel="Insulin (μU/mL)", legend=:topright)
annotate!(subplot=2, 90, 0.5, text("No insulin secretion!", :red, 12))

X_t1d = [u[2] for u in sol_t1d.u]
plot!(subplot=3, sol_t1d.t, X_t1d, lw=3, label="Insulin Action", color=:green,
      xlabel="Time (min)", ylabel="Insulin Action (min⁻¹)", legend=:topright)
annotate!(subplot=3, 90, 0.0001, text("No insulin action!", :green, 12))

savefig(p2, "../results/02_t1d_no_insulin.png")
println("    ✓ Saved: results/02_t1d_no_insulin.png")

# 5.3 Three-way comparison: Normal, T1D untreated, T1D treated
p3 = plot(sol_normal.t, G_normal, lw=3, label="Normal", color=:blue,
          xlabel="Time (min)", ylabel="Total Glucose (mg/dL)",
          title="Type 1 Diabetes: Effect of Insulin Therapy",
          legend=:topright, size=(1000, 600))
plot!(sol_t1d.t, G_t1d, lw=3, label="T1D (no treatment)", color=:red, ls=:dash)
plot!(sol_t1d_therapy.t, G_t1d_therapy, lw=3, label="T1D (with insulin)", 
      color=:green, ls=:dashdot)
hline!([140], lw=2, ls=:dot, color=:gray, label="Normal threshold", alpha=0.5)
savefig(p3, "../results/02_t1d_comparison_three_way.png")
println("    ✓ Saved: results/02_t1d_comparison_three_way.png")

# 5.4 Insulin profile with therapy
p4 = plot(sol_t1d_therapy.t, I_t1d_therapy, lw=3, color=:red,
          xlabel="Time (min)", ylabel="Insulin (μU/mL)",
          title="Type 1 Diabetes: Exogenous Insulin Profile",
          label="Insulin (basal + bolus)", legend=:topright,
          size=(1000, 500))
annotate!(30, maximum(I_t1d_therapy)*0.9, text("Bolus\n(meal-time)", :red, 10))
annotate!(120, 5, text("← Basal (background)", :red, 10))
savefig(p4, "../results/02_t1d_insulin_therapy.png")
println("    ✓ Saved: results/02_t1d_insulin_therapy.png")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 70)
println("SUMMARY: Type 1 Diabetes")
println("=" ^ 70)
println("\nPathophysiology:")
println("  ✗ Autoimmune destruction of pancreatic β-cells")
println("  ✗ k5 = 0 → No endogenous insulin secretion")
println("  ✗ Open-loop system: no feedback control")
println("\nWithout treatment:")
println("  ✗ Peak glucose: $(@sprintf("%.0f", peak_G_t1d)) mg/dL (normal: ~140)")
println("  ✗ No glucose clearance")
println("  ✗ Life-threatening hyperglycemia")
println("  ✗ Risk of diabetic ketoacidosis (DKA)")
println("\nWith insulin therapy:")
println("  ✓ Peak glucose: $(@sprintf("%.0f", peak_G_therapy)) mg/dL")
println("  ✓ Glucose control restored (near-normal)")
println("  ✓ Exogenous insulin replaces missing feedback signal")
println("\nControl theory perspective:")
println("  • Type 1 diabetes = SENSOR FAILURE")
println("  • Pancreas can no longer detect/respond to glucose")
println("  • Treatment = Manual/automated open-loop control")
println("  • Artificial pancreas = Closing the loop with technology")
println("\nClinical implications:")
println("  • Insulin therapy is MANDATORY for survival")
println("  • Requires careful dosing (basal + bolus)")
println("  • Risk: Hypoglycemia from excessive insulin")
println("  • Goal: Time-in-range 70-180 mg/dL >70% of the time")
println("=" ^ 70)
