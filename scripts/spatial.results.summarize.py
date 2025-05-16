import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import scikit_posthocs as sp

# Load the CSV file
file_path = "results/cell2location/merged_all.csv"  # Replace with your actual path
df = pd.read_csv(file_path)

# 1. Summary Statistics
summary_stats = {
    "Total Rows (Spots)": df.shape[0],
    "Unique patient_id": df["patient_id"].nunique() if "patient_id" in df.columns else None,
    "Unique sample_id": df["sample_id"].nunique() if "sample_id" in df.columns else None,
    "Unique Condition": df["Condition"].nunique() if "Condition" in df.columns else None,
    "Unique TissueStage": df["TissueStage"].nunique() if "TissueStage" in df.columns else None,
    "Unique UICC": df["UICC"].nunique() if "UICC" in df.columns else None,
}
print("Summary Statistics:")
for k, v in summary_stats.items():
    print(f"{k}: {v}")

# 2. Grouping and Averaging
tam_columns = ["RTM_TAMs", "Prolif_TAMs", "Angio_TAMs", "LA_TAMs", 
               "IFN_TAMs", "Inflam_TAMs", "Reg_TAMs"]
grouped_means = df.groupby(["Condition", "TissueStage"])[tam_columns].mean().reset_index()
print("\nGrouped Averages:\n", grouped_means)

# 3. Prepare data for plotting
melted_df = df.melt(id_vars=["Condition", "TissueStage"], 
                    value_vars=tam_columns, 
                    var_name="TAM_Type", 
                    value_name="Proportion")

# Ensure order of TissueStage for consistent plot coloring
melted_df["TissueStage"] = pd.Categorical(
    melted_df["TissueStage"],
    categories=["Healthy", "Background", "Tumour"],
    ordered=True
)

# Run Dunn's post-hoc test with Bonferroni correction for each TAM type and collect significant pairs
significance_marks = []

for tam in tam_columns:
    posthoc = sp.posthoc_dunn(
        df[df["TissueStage"].isin(["Healthy", "Background", "Tumour"])][["TissueStage", tam]].dropna(),
        val_col=tam,
        group_col="TissueStage",
        p_adjust="bonferroni"
    )

    # Mark significant comparisons with asterisks
    for i, group1 in enumerate(posthoc.index):
        for j, group2 in enumerate(posthoc.columns):
            if j <= i:
                continue
            p = posthoc.iloc[i, j]
            if p < 0.001:
                mark = "***"
            elif p < 0.01:
                mark = "**"
            elif p < 0.05:
                mark = "*"
            else:
                continue
            significance_marks.append((tam, group1, group2, mark))

# Plot with significance annotations
plt.figure(figsize=(14, 8))
ax = sns.boxplot(data=melted_df, x="TAM_Type", y="Proportion", hue="TissueStage")

# Add asterisks for significant differences
for i, (tam, group1, group2, mark) in enumerate(significance_marks):
    xpos = tam_columns.index(tam)
    ypos = df[df["TissueStage"].isin([group1, group2])][tam].max() + 0.1 + i * 0.02
    x1 = xpos - 0.2 + 0.2 * ["Healthy", "Background", "Tumour"].index(group1)
    x2 = xpos - 0.2 + 0.2 * ["Healthy", "Background", "Tumour"].index(group2)
    ax.plot([x1, x1, x2, x2], [ypos - 0.01, ypos, ypos, ypos - 0.01], lw=1.2, c='k')
    ax.text((x1 + x2) / 2, ypos + 0.01, mark, ha='center', va='bottom', color='k', fontsize=10)

# Final plot adjustments
ax.set_title("TAM Subtype Proportions with Statistical Significance", fontsize=16)
ax.set_xlabel("TAM Subtype", fontsize=12)
ax.set_ylabel("Proportion", fontsize=12)
plt.xticks(rotation=45)
plt.legend(title="TissueStage", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# 5. Save the plot as PDF
pdf_path = "figures/tam_boxplot_with_tissue_stage.pdf"
plt.savefig(pdf_path, format='pdf')
plt.show()
