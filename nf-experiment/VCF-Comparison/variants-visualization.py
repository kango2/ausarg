# Variants Visualization Script
# This script provides a comprehensive visualization of variant comparison metrics
# between SNPs and INDELs. The visualizations include:
# - Violin plots comparing the distribution of metrics like TP, FP, UNK, Precision, Recall, and F1 Score for SNPs and INDELs.
# - Pair plots providing scatter plots and histograms to understand bivariate relationships and individual metric distributions for SNPs and INDELs.
# - Heatmap plots showing the correlation between the metrics for SNPs and INDELs.
# - Joint plots illustrating the relationships between specific metric pairs for SNPs and INDELs.

# Authors: Kosar Hooshmand

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

# Load data
# "pass" dataset
snp_pass_data = pd.read_csv("/g/data/te53/kh3349/varaints-comparison/hap.py_output/SNP.roc.Locations.SNP.PASS.csv.gz", compression='gzip')
indel_pass_data = pd.read_csv("/g/data/te53/kh3349/varaints-comparison/hap.py_output/INDEL.roc.Locations.INDEL.PASS.csv.gz", compression='gzip')
snp_pass_data['Type'] = 'SNP'
indel_pass_data['Type'] = 'INDEL'
combined_pass_data = pd.concat([snp_pass_data, indel_pass_data])

# "all" dataset
snp_all_data = pd.read_csv("/g/data/te53/kh3349/varaints-comparison/hap.py_output/SNP.roc.Locations.SNP.csv.gz", compression='gzip')
indel_all_data = pd.read_csv("/g/data/te53/kh3349/varaints-comparison/hap.py_output/INDEL.roc.Locations.INDEL.csv.gz", compression='gzip')
snp_all_data['Type'] = 'SNP'
indel_all_data['Type'] = 'INDEL'
combined_all_data = pd.concat([snp_all_data, indel_all_data])

# Output directory
output_directory = "/g/data/te53/kh3349/varaints-comparison"

# Simplified names for metrics
simplified_names = {
    "QUERY.TP": "TP",
    "QUERY.FP": "FP",
    "QUERY.UNK": "UNK",
    "METRIC.Precision": "Precision",
    "METRIC.Recall": "Recall",
    "METRIC.F1_Score": "F1"
}

# Metrics to plot
metrics_to_plot = ["QUERY.TP", "QUERY.FP", "QUERY.UNK", "METRIC.Precision", "METRIC.Recall", "METRIC.F1_Score"]

# Metric pairs for joint plots
metric_pairs = [("METRIC.Precision", "METRIC.Recall"), ("QUERY.TP", "QUERY.FP")]

def generate_visualizations(data, prefix):
    # Melted data for violin plots
    melted_data = data.melt(id_vars=['Type'], value_vars=metrics_to_plot, var_name='Metric', value_name='Value')

# ---- Heatmap Plots ----
for dtype in ['SNP', 'INDEL']:
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    data_pass = snp_pass_data if dtype == 'SNP' else indel_pass_data
    data_all = snp_all_data if dtype == 'SNP' else indel_all_data

    sns.heatmap(data_pass[metrics_to_plot].corr(), annot=True, cmap="coolwarm", 
                xticklabels=[simplified_names[m] for m in metrics_to_plot], 
                yticklabels=[simplified_names[m] for m in metrics_to_plot], ax=axes[0])
    axes[0].set_title(f'Correlation Heatmap for {dtype} (Pass)')
    
    sns.heatmap(data_all[metrics_to_plot].corr(), annot=True, cmap="coolwarm", 
                xticklabels=[simplified_names[m] for m in metrics_to_plot], 
                yticklabels=[simplified_names[m] for m in metrics_to_plot], ax=axes[1])
    axes[1].set_title(f'Correlation Heatmap for {dtype} (All)')

    plt.tight_layout()
    plt.savefig(f"{output_directory}/heatmap_{dtype}_combined.png", dpi=300)
    plt.close()

# ---- Joint Plots ----   

# Plotting individual joint plots
def plot_individual_jointplots(data, label, variant_type):
    cmap_choice = 'Blues' if variant_type == 'SNPs' else 'Reds'
    for metric_x, metric_y in metric_pairs:
        sns.jointplot(data=data, x=metric_x, y=metric_y, kind='kde', fill=True, cmap=cmap_choice).fig.suptitle(f'{variant_type} ({label}): {metric_x.split(".")[-1]} vs {metric_y.split(".")[-1]}')
        plt.tight_layout()
        plt.savefig(f"{output_directory}/jointplot_{variant_type}_{label}_{metric_x.split('.')[-1]}_vs_{metric_y.split('.')[-1]}.png", dpi=300)
        plt.close()

# Generate individual joint plots for SNPs
plot_individual_jointplots(snp_all_data, "All", "SNPs")
plot_individual_jointplots(snp_pass_data, "Pass", "SNPs")

# Generate individual joint plots for INDELs
plot_individual_jointplots(indel_all_data, "All", "INDELs")
plot_individual_jointplots(indel_pass_data, "Pass", "INDELs")

# ---- Pair Plots ----

# Generate individual pair plots for SNPs
sns.pairplot(data=snp_all_data, vars=metrics_to_plot, diag_kind="kde", markers="o", corner=True, plot_kws={'alpha':0.5}).fig.suptitle("SNPs (All)")
plt.tight_layout()
plt.savefig(f"{output_directory}/pairplot_SNPs_All.png", dpi=300)
plt.close()

sns.pairplot(data=snp_pass_data, vars=metrics_to_plot, diag_kind="kde", markers="o", corner=True, plot_kws={'alpha':0.5}).fig.suptitle("SNPs (Pass)")
plt.tight_layout()
plt.savefig(f"{output_directory}/pairplot_SNPs_Pass.png", dpi=300)
plt.close()

# Generate individual pair plots for INDELs
sns.pairplot(data=indel_all_data, vars=metrics_to_plot, diag_kind="kde", markers="o", corner=True, plot_kws={'alpha':0.5}).fig.suptitle("INDELs (All)")
plt.tight_layout()
plt.savefig(f"{output_directory}/pairplot_INDELs_All.png", dpi=300)
plt.close()

sns.pairplot(data=indel_pass_data, vars=metrics_to_plot, diag_kind="kde", markers="o", corner=True, plot_kws={'alpha':0.5}).fig.suptitle("INDELs (Pass)")
plt.tight_layout()
plt.savefig(f"{output_directory}/pairplot_INDELs_Pass.png", dpi=300)
plt.close()

# ---- Violin Plots ----

# Metrics to plot
metrics_to_plot = ["QUERY.TP", "QUERY.FP", "QUERY.UNK", "METRIC.Precision", "METRIC.Recall", "METRIC.F1_Score"]

# Prepare the melted data for SNPs
melted_all_snp_data = snp_all_data.melt(value_vars=metrics_to_plot, var_name='Metric', value_name='Value')
melted_all_snp_data['Dataset'] = 'All'
melted_pass_snp_data = snp_pass_data.melt(value_vars=metrics_to_plot, var_name='Metric', value_name='Value')
melted_pass_snp_data['Dataset'] = 'Pass'
melted_combined_snp_data = pd.concat([melted_all_snp_data, melted_pass_snp_data]).reset_index(drop=True)

# Prepare the melted data for INDELs
melted_all_indel_data = indel_all_data.melt(value_vars=metrics_to_plot, var_name='Metric', value_name='Value')
melted_all_indel_data['Dataset'] = 'All'
melted_pass_indel_data = indel_pass_data.melt(value_vars=metrics_to_plot, var_name='Metric', value_name='Value')
melted_pass_indel_data['Dataset'] = 'Pass'
melted_combined_indel_data = pd.concat([melted_all_indel_data, melted_pass_indel_data]).reset_index(drop=True)

# Plot for SNPs
plt.figure(figsize=(20, 10))
for idx, metric in enumerate(metrics_to_plot, 1):
    plt.subplot(2, 3, idx)
    sns.violinplot(x='Metric', y='Value', hue='Dataset', data=melted_combined_snp_data[melted_combined_snp_data['Metric'] == metric], palette='muted', split=True)
    plt.title(f"SNPs: {simplified_names[metric]}")
plt.tight_layout()
plt.savefig(f"{output_directory}/SNPs_violin_plots_combined.png", dpi=300)
plt.close()

# Plot for INDELs
plt.figure(figsize=(20, 10))
for idx, metric in enumerate(metrics_to_plot, 1):
    plt.subplot(2, 3, idx)
    sns.violinplot(x='Metric', y='Value', hue='Dataset', data=melted_combined_indel_data[melted_combined_indel_data['Metric'] == metric], palette='muted', split=True)
    plt.title(f"INDELs: {simplified_names[metric]}")
plt.tight_layout()
plt.savefig(f"{output_directory}/INDELs_violin_plots_combined.png", dpi=300)
plt.close()
