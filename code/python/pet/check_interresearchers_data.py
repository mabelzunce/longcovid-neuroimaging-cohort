import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set the data path
repo_path = "/home/martin/data/UNSAM/CovidProject/longcovid-neuroimaging-cohort/data/processed_pet/"
output_fig_path = os.path.join(repo_path, "suvr_vs_pet_scatter_plots")
# Create the output directory if it doesn't exist
os.makedirs(output_fig_path, exist_ok=True)
# Read the CSV files
pet_file = os.path.join(repo_path, "PET_CERMEP_CEUNIM_merged.csv")
suvr_file = os.path.join(repo_path, "merged_summary_csv", "summary_suvr_mean_values.csv")

pet_df = pd.read_csv(pet_file)
suvr_df = pd.read_csv(suvr_file)

# Keep only rows with matching IDs in both dataframes, in the same order
common_ids = set(pet_df['ID']).intersection(set(suvr_df['ID']))
pet_df = pet_df[pet_df['ID'].isin(common_ids)].set_index('ID').loc[sorted(common_ids)].reset_index()
suvr_df = suvr_df[suvr_df['ID'].isin(common_ids)].set_index('ID').loc[sorted(common_ids)].reset_index()
# Identify numeric columns with 'suvr'
metadata_cols = {'ID', 'group', 'Group', 'site', 'SITE', 'age', 'sex', 'sexo', 'edad', 'weight', 'anonimizado'}
# Keep only the numeric regional uptake measures; metadata columns are not valid predictors.
suvr_columns = [col for col in suvr_df.columns if col not in metadata_cols]
feature_df = suvr_df[suvr_columns].apply(pd.to_numeric, errors='coerce')
# print columnd names and shapes
print("SUVR DataFrame columns:", suvr_columns)
suvr_values = feature_df.to_numpy(dtype=float)
# Repeat for pet_df to identify numeric columns
numeric_cols_pet = [col for col in pet_df.columns if col not in metadata_cols]
# print columnd names and shapes
print("PET DataFrame columns:", numeric_cols_pet)
feature_df = pet_df[numeric_cols_pet].apply(pd.to_numeric, errors='coerce')
pet_values = feature_df.to_numpy(dtype=float)

# print shapes
print("PET DataFrame shape:", pet_df.shape)
print("SUVR DataFrame shape:", suvr_df.shape)

# Simple plot, convert matrix into a vector and plot one against the other
plt.figure(figsize=(10, 6))
# Flatten the matrices to 1D arrays for plotting
suvr_values_flat = suvr_values.flatten()
pet_values_flat = pet_values.flatten()
# Create a scatter plot
sns.scatterplot(x=suvr_values_flat, y=pet_values_flat, alpha=0.5)
plt.title('Scatter Plot of SUVR Values vs PET Values')
plt.xlabel('SUVR Values')
plt.ylabel('PET Values')
plt.grid(True)
# Save the plot
output_plot_path = os.path.join(output_fig_path, "suvr_vs_pet_scatter_plot.png")
plt.savefig(output_plot_path)

# Now go through each column in suvr_columns and plot against the corresponding column in pet_df, if it exists
for col in suvr_columns:
    if col in pet_df.columns:
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x=suvr_df[col], y=pet_df[col], alpha=0.5)
        plt.title(f'Scatter Plot of {col} SUVR vs PET Values')
        plt.xlabel(f'{col} SUVR Values')
        plt.ylabel(f'{col} PET Values')
        plt.grid(True)
        # Save the plot
        output_plot_path = os.path.join(output_fig_path, f"{col}_suvr_vs_pet_scatter_plot.png")
        plt.savefig(output_plot_path)
        plt.close()

# Now do the same per rows
for idx in range(pet_df.shape[0]):
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=suvr_values[idx, :], y=pet_values[idx, :], alpha=0.5)
    plt.title(f'Scatter Plot of SUVR vs PET Values for ID {pet_df.loc[idx, "ID"]}')
    plt.xlabel('SUVR Values')
    plt.ylabel('PET Values')
    plt.grid(True)
    # Save the plot
    output_plot_path = os.path.join(output_fig_path, f"{pet_df.loc[idx, 'ID']}_suvr_vs_pet_scatter_plot.png")
    plt.savefig(output_plot_path)
    plt.close()