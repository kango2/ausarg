# Redefine the column names for the first 12 columns
column_names = [
    "Query sequence name", "Query sequence length", "Query start", "Query end", "Relative strand",
    "Target sequence name", "Target sequence length", "Target start", "Target end", "Number of residue matches",
    "Alignment block length", "Mapping quality"
]

# Reading the file with more control over parsing
data = []
with open(file_path, 'r') as file:
    for line in file:
        # Split the line by tab
        parts = line.strip().split('\t')
        data.append(parts)

# Convert the list of lists to a DataFrame
data_df = pd.DataFrame(data)

# Assign column names to the first 12 columns
data_df.columns = column_names + [f'Var_col_{i}' for i in range(1, len(data_df.columns) - 12 + 1)]

# Extract remaining columns for tags
remaining_columns = data_df.iloc[:, 12:]

# Parse remaining columns into separate tag columns
tags = remaining_columns.apply(lambda x: x.str.split(':', expand=True).stack().reset_index(level=1, drop=True)).reset_index()
tags.columns = ['Row', 'Tag', 'Type', 'Value']

# Pivot tags DataFrame to get tags as columns
tags_pivot = tags.pivot(index='Row', columns='Tag', values='Value').reset_index(drop=True)

# Combine fixed columns and tag columns
fixed_columns = data_df.iloc[:, :12]
fixed_columns.columns = column_names
result = pd.concat([fixed_columns, tags_pivot], axis=1)

# Filter for 18S and 28S sequences
result['Query type'] = result['Query sequence name'].apply(lambda x: '18S' if '_18S' in x else '28S')

# Convert target start, target end and divergence to numeric types
result['Target start'] = result['Target start'].astype(int)
result['Target end'] = result['Target end'].astype(int)
result['dv'] = result['dv'].astype(float)

# Identify the least divergent alignments for each target location
grouped = result.groupby(['Target sequence name', 'Target start', 'Target end', 'Query type'])
least_divergent = grouped.apply(lambda x: x.loc[x['dv'].idxmin()]).reset_index(drop=True)

# Extract relevant columns for output
output = least_divergent[['Target sequence name', 'Target start', 'Target end', 'Query type', 'dv']]

import ace_tools as tools; tools.display_dataframe_to_user(name="rDNA Alignments", dataframe=output)

output.head()
