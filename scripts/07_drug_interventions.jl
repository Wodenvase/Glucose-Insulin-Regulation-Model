"""
Example 7: Drug Interventions (Pharmacological Therapy)

Compares different drug classes for Type 2 diabetes treatment.

Drugs simulated:
1. Metformin (↑ glucose effectiveness)
2. Sulfonylurea (↑ insulin secretion)
3. TZD/Thiazolidinedione (↑ insulin sensitivity)
4. Combination therapy

Shows how each drug modulates specific model parameters.

Expected outcomes:
- Metformin: Modest glucose reduction, safe
- Sulfonylurea: Larger effect, hypoglycemia risk
- TZD: Addresses root cause (insulin resistance)
- Combination: Additive effects, best control
"""

using Plots
using Printf

include("../src/GlucoseInsulinModel.jl")
using .GlucoseInsulinModel

mkpath("../results")

println("=" ^ 70)
println("Example 7: Drug Interventions for Type 2 Diabetes")
println("=" ^ 70)

# ============================================================================
# 1. Setup Baseline Type 2 Diabetes
# ============================================================================

println("\n1. Creating baseline Type 2 diabetes parameters...")

params_t2d_baseline = BergmanParams(
    k1 = 0.028,
    k2 = 0.025,
    k3 = 1.0e-5,    # 40% of normal (insulin resistance)
    k4 = 0.05,
    k5 = 0.010,     # 67% of normal (β-cell dysfunction)
    Gb = 81.0
)

println("  Baseline T2D parameters:")
println("    k1 (glucose effectiveness):     $(params_t2d_baseline.k1) min⁻¹")
println("    k3 (insulin sensitivity):       $(params_t2d_baseline.k3) min⁻² per μU/mL (40% of normal)")
println("    k5 (pancreatic responsiveness): $(params_t2d_baseline.k5) min⁻¹ (67% of normal)")

# ============================================================================
# 2. Define Drug Effects
# ============================================================================

println("\n2. Defining drug effects on parameters...")

# Drug 1: Metformin (increases glucose effectiveness)
println("\n  Drug 1: METFORMIN")
println("    Mechanism: ↑ Glucose effectiveness, ↓ hepatic glucose production")
println("    Model effect: k1 → 1.5 × k1")

params_metformin = BergmanParams(
    k1 = 1.5 * params_t2d_baseline.k1,  # Increase glucose effectiveness
    k2 = params_t2d_baseline.k2,
    k3 = params_t2d_baseline.k3,
    k4 = params_t2d_baseline.k4,
    k5 = params_t2d_baseline.k5,
    Gb = params_t2d_baseline.Gb
)

# Drug 2: Sulfonylurea (increases insulin secretion)
println("\n  Drug 2: SULFONYLUREA (e.g., Glyburide)")
println("    Mechanism: Stimulates β-cell insulin secretion")
println("    Model effect: k5 → 1.5 × k5")

params_sulfonylurea = BergmanParams(
    k1 = params_t2d_baseline.k1,
    k2 = params_t2d_baseline.k2,
    k3 = params_t2d_baseline.k3,
    k4 = params_t2d_baseline.k4,
    k5 = 1.5 * params_t2d_baseline.k5,  # Increase insulin secretion
    Gb = params_t2d_baseline.Gb
)

# Drug 3: TZD/Pioglitazone (increases insulin sensitivity)
println("\n  Drug 3: TZD (e.g., Pioglitazone)")
println("    Mechanism: PPAR-γ agonist, ↑ insulin sensitivity")
println("    Model effect: k3 → 1.8 × k3")

params_tzd = BergmanParams(
    k1 = params_t2d_baseline.k1,
    k2 = params_t2d_baseline.k2,
    k3 = 1.8 * params_t2d_baseline.k3,  # Increase insulin sensitivity
    k4 = params_t2d_baseline.k4,
    k5 = params_t2d_baseline.k5,
    Gb = params_t2d_baseline.Gb
)

# Drug 4: Combination (Metformin + Sulfonylurea)
println("\n  Drug 4: COMBINATION (Metformin + Sulfonylurea)")
println("    Model effect: k1 → 1.5 × k1 AND k5 → 1.5 × k5")

params_combination = BergmanParams(
    k1 = 1.5 * params_t2d_baseline.k1,
    k2 = params_t2d_baseline.k2,
    k3 = params_t2d_baseline.k3,
    k4 = params_t2d_baseline.k4,
    k5 = 1.5 * params_t2d_baseline.k5,
    Gb = params_t2d_baseline.Gb
)

# ============================================================================
# 3. Simulate OGTT for Each Treatment
# ============================================================================

println("\n3. Simulating OGTT for each treatment...")

D_func = create_glucose_input(:ogtt)
I_ext = create_insulin_input(:none)
u0 = [0.0, 0.0, 0.0]
tspan = (0.0, 180.0)

# Baseline (no treatment)
println("  • Simulating baseline (no treatment)...")
sol_baseline = simulate_bergman(u0, tspan, params_t2d_baseline, D_func, I_ext)
G_baseline = [u[1] + params_t2d_baseline.Gb for u in sol_baseline.u]

# Metformin
println("  • Simulating with Metformin...")
sol_metformin = simulate_bergman(u0, tspan, params_metformin, D_func, I_ext)
G_metformin = [u[1] + params_metformin.Gb for u in sol_metformin.u]

# Sulfonylurea
println("  • Simulating with Sulfonylurea...")
sol_sulfonylurea = simulate_bergman(u0, tspan, params_sulfonylurea, D_func, I_ext)
G_sulfonylurea = [u[1] + params_sulfonylurea.Gb for u in sol_sulfonylurea.u]

# TZD
println("  • Simulating with TZD...")
sol_tzd = simulate_bergman(u0, tspan, params_tzd, D_func, I_ext)
G_tzd = [u[1] + params_tzd.Gb for u in sol_tzd.u]

# Combination
println("  • Simulating with Combination...")
sol_combination = simulate_bergman(u0, tspan, params_combination, D_func, I_ext)
G_combination = [u[1] + params_combination.Gb for u in sol_combination.u]

# ============================================================================
# 4. Analyze Drug Efficacy
# ============================================================================

println("\n4. Analyzing drug efficacy...")

function analyze_efficacy(sol, params, name)
    G_vals = [u[1] + params.Gb for u in sol.u]
    
    peak_G = maximum(G_vals)
    
    idx_120 = findfirst(t -> t >= 120.0, sol.t)
    G_120 = G_vals[idx_120]
    
    final_G = G_vals[end]
    
    # Calculate AUC (glucose exposure)
    AUC = sum(G_vals) * (sol.t[2] - sol.t[1])
    
    return (peak=peak_G, G_120=G_120, final=final_G, AUC=AUC)
end

results_baseline = analyze_efficacy(sol_baseline, params_t2d_baseline, "Baseline")
results_metformin = analyze_efficacy(sol_metformin, params_metformin, "Metformin")
results_sulfonylurea = analyze_efficacy(sol_sulfonylurea, params_sulfonylurea, "Sulfonylurea")
results_tzd = analyze_efficacy(sol_tzd, params_tzd, "TZD")
results_combination = analyze_efficacy(sol_combination, params_combination, "Combination")

println("\n  Efficacy Summary:")
println("  " * "─"^70)
println("  Treatment      | Peak (mg/dL) | 2-hr (mg/dL) | AUC (rel) |")
println("  " * "─"^70)

function print_result(name, result, baseline_auc)
    auc_reduction = (1 - result.AUC/baseline_auc) * 100
    println("  $(rpad(name, 14)) | $(lpad(@sprintf("%.0f", result.peak), 12)) | $(lpad(@sprintf("%.0f", result.G_120), 12)) | $(lpad(@sprintf("%+.0f%%", auc_reduction), 9)) |")
end

print_result("Baseline", results_baseline, results_baseline.AUC)
print_result("Metformin", results_metformin, results_baseline.AUC)
print_result("Sulfonylurea", results_sulfonylurea, results_baseline.AUC)
print_result("TZD", results_tzd, results_baseline.AUC)
print_result("Combination", results_combination, results_baseline.AUC)
println("  " * "─"^70)

# ============================================================================
# 5. Visualizations
# ============================================================================

println("\n5. Creating visualizations...")

# 5.1 All treatments comparison
p1 = plot(sol_baseline.t, G_baseline, lw=3, label="Baseline (no Rx)", color=:black,
          xlabel="Time (min)", ylabel="Total Glucose (mg/dL)",
          title="Drug Interventions: Glucose Response Comparison",
          legend=:topright, size=(1000, 600), linestyle=:dash)
plot!(sol_metformin.t, G_metformin, lw=3, label="Metformin", color=:blue)
plot!(sol_sulfonylurea.t, G_sulfonylurea, lw=3, label="Sulfonylurea", color=:orange)
plot!(sol_tzd.t, G_tzd, lw=3, label="TZD", color=:green)
plot!(sol_combination.t, G_combination, lw=3, label="Combination", color=:purple)
hline!([140], lw=2, ls=:dot, color=:gray, label="Normal threshold", alpha=0.5)
savefig(p1, "../results/07_drug_comparison_all.png")
println("    ✓ Saved: results/07_drug_comparison_all.png")

# 5.2 Individual drug effects (2x2 grid)
p2 = plot(
    plot(sol_baseline.t, G_baseline, lw=2, color=:black, ls=:dash,
         label="Baseline", title="Metformin", ylabel="Glucose (mg/dL)"),
    plot(sol_baseline.t, G_baseline, lw=2, color=:black, ls=:dash,
         label="Baseline", title="Sulfonylurea"),
    plot(sol_baseline.t, G_baseline, lw=2, color=:black, ls=:dash,
         label="Baseline", title="TZD", xlabel="Time (min)", ylabel="Glucose (mg/dL)"),
    plot(sol_baseline.t, G_baseline, lw=2, color=:black, ls=:dash,
         label="Baseline", title="Combination", xlabel="Time (min)"),
    layout=(2,2), size=(1200, 800)
)
plot!(p2[1], sol_metformin.t, G_metformin, lw=3, color=:blue, label="With Drug")
plot!(p2[2], sol_sulfonylurea.t, G_sulfonylurea, lw=3, color=:orange, label="With Drug")
plot!(p2[3], sol_tzd.t, G_tzd, lw=3, color=:green, label="With Drug")
plot!(p2[4], sol_combination.t, G_combination, lw=3, color=:purple, label="With Drug")

savefig(p2, "../results/07_drug_individual_effects.png")
println("    ✓ Saved: results/07_drug_individual_effects.png")

# 5.3 Efficacy bar charts
treatments = ["Baseline", "Metformin", "Sulfonylurea", "TZD", "Combination"]
peak_values = [results_baseline.peak, results_metformin.peak, 
               results_sulfonylurea.peak, results_tzd.peak, results_combination.peak]
G_120_values = [results_baseline.G_120, results_metformin.G_120,
                results_sulfonylurea.G_120, results_tzd.G_120, results_combination.G_120]

data_matrix = hcat(peak_values, G_120_values)
p3 = bar(treatments, data_matrix,
         bar_position=:dodge,
         label=["Peak Glucose", "2-hr Glucose"],
         xlabel="Treatment", ylabel="Glucose (mg/dL)",
         title="Drug Efficacy Comparison",
         xticks=(1:5, treatments), xrotation=20,
         legend=:topright, size=(1000, 600),
         color=[:steelblue :coral])
hline!([140], lw=2, ls=:dash, color=:green, label="Normal threshold")
savefig(p3, "../results/07_drug_efficacy_bars.png")
println("    ✓ Saved: results/07_drug_efficacy_bars.png")

# 5.4 AUC reduction
auc_reductions = [(results_baseline.AUC - r.AUC)/results_baseline.AUC * 100 
                  for r in [results_baseline, results_metformin, results_sulfonylurea, 
                           results_tzd, results_combination]]

p4 = bar(treatments, auc_reductions, 
         xlabel="Treatment", ylabel="AUC Reduction (%)",
         title="Glucose Exposure Reduction",
         legend=false, color=:forestgreen,
         xticks=(1:5, treatments), xrotation=20,
         size=(900, 600))
savefig(p4, "../results/07_auc_reduction.png")
println("    ✓ Saved: results/07_auc_reduction.png")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 70)
println("SUMMARY: Drug Interventions")
println("=" ^ 70)
println("\n1. METFORMIN (First-line therapy)")
println("   Mechanism: ↑ Glucose effectiveness (k1 × 1.5)")
println("   Effect:")
println("     • Peak glucose: $(@sprintf("%.0f", results_metformin.peak)) mg/dL (vs $(@sprintf("%.0f", results_baseline.peak)) baseline)")
println("     • AUC reduction: $(@sprintf("%.0f", auc_reductions[2]))%")
println("   Pros: Safe, no hypoglycemia, CV benefits")
println("   Cons: Modest efficacy, GI side effects")
println("\n2. SULFONYLUREA (β-cell stimulant)")
println("   Mechanism: ↑ Insulin secretion (k5 × 1.5)")
println("   Effect:")
println("     • Peak glucose: $(@sprintf("%.0f", results_sulfonylurea.peak)) mg/dL")
println("     • AUC reduction: $(@sprintf("%.0f", auc_reductions[3]))%")
println("   Pros: Effective, inexpensive")
println("   Cons: Hypoglycemia risk, weight gain, β-cell stress")
println("\n3. TZD/Pioglitazone (insulin sensitizer)")
println("   Mechanism: ↑ Insulin sensitivity (k3 × 1.8)")
println("   Effect:")
println("     • Peak glucose: $(@sprintf("%.0f", results_tzd.peak)) mg/dL")
println("     • AUC reduction: $(@sprintf("%.0f", auc_reductions[4]))%")
println("   Pros: Addresses root cause (insulin resistance)")
println("   Cons: Weight gain, fracture risk, slow onset")
println("\n4. COMBINATION (Metformin + Sulfonylurea)")
println("   Mechanism: ↑ k1 AND ↑ k5")
println("   Effect:")
println("     • Peak glucose: $(@sprintf("%.0f", results_combination.peak)) mg/dL")
println("     • AUC reduction: $(@sprintf("%.0f", auc_reductions[5]))%")
println("   Pros: Additive/synergistic effects, best glucose control")
println("   Cons: Multiple side effects, polypharmacy")
println("\nRanking (by efficacy):")
println("  1. Combination therapy (AUC ↓ $(@sprintf("%.0f", auc_reductions[5]))%)")
println("  2. TZD (AUC ↓ $(@sprintf("%.0f", auc_reductions[4]))%)")
println("  3. Sulfonylurea (AUC ↓ $(@sprintf("%.0f", auc_reductions[3]))%)")
println("  4. Metformin (AUC ↓ $(@sprintf("%.0f", auc_reductions[2]))%)")
println("\nRanking (by safety):")
println("  1. Metformin (no hypoglycemia)")
println("  2. TZD (no hypoglycemia)")
println("  3. Combination (some hypoglycemia risk)")
println("  4. Sulfonylurea (highest hypoglycemia risk)")
println("\nClinical Decision-Making:")
println("  • First-line: Metformin (safety + modest efficacy)")
println("  • Add-on therapy depends on phenotype:")
println("    - Predominant insulin resistance → TZD or GLP-1")
println("    - β-cell dysfunction → Sulfonylurea or insulin")
println("  • Combination therapy for inadequate monotherapy response")
println("  • Model-based medicine: Estimate k3, k5 → tailor therapy")
println("\nPrecision Endocrinology:")
println("  • Measure patient-specific k3 (insulin sensitivity) from IVGTT")
println("  • Measure k5 (β-cell function) from OGTT")
println("  • Simulate drug effects IN SILICO before prescribing")
println("  • Optimize regimen for individual parameter deficits")
println("  • This is the future: personalized, quantitative medicine")
println("=" ^ 70)
