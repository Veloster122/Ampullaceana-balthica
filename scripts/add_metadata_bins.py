#!/usr/bin/env python3
import pandas as pd
import sys

def main():
    input_file = "00-Ampullaceana_balthica_raw_data/Ampullaceana_balthica_metadata.tsv"
    output_file = "00-Ampullaceana_balthica_raw_data/Ampullaceana_balthica_metadata_binned.tsv"
    
    # Read all lines
    with open(input_file, 'r') as f:
        lines = f.readlines()
        
    headers = lines[0].strip().split('\t')
    types = lines[1].strip().split('\t')
    
    # Load into dataframe skipping the type row, forcing 'ID' to be read as string
    df = pd.read_csv(input_file, sep='\t', skiprows=[1], dtype={'ID': str})
    
    new_columns = {}
    new_types = []
    
    # Find numeric columns
    for i, (col, qtype) in enumerate(zip(headers, types)):
        if qtype == 'numeric' and col != 'ID' and not col.startswith('S'): # Skip ID columns if any
            try:
                # Convert to numeric, forcing errors to NaN for safety
                numeric_series = pd.to_numeric(df[col], errors='coerce')
                
                # Attempt to bin into 3 quantiles. If too many duplicate values (e.g., lots of zeros),
                # duplicates='drop' will reduce the number of bins. Returns interval objects.
                binned = pd.qcut(numeric_series, q=3, duplicates='drop')
                
                # Convert intervals to string so it's categorical
                new_col_name = f"{col}_bin"
                # Rename the interval categories to something readable
                # (Low, Med, High) logic is tricky if bins drop.
                # So we just use the interval string (e.g. "(0.0, 15.2]")
                df[new_col_name] = binned.astype(str)
                # Replace 'nan' string with empty or 'Missing'
                df.loc[df[new_col_name] == 'nan', new_col_name] = 'Missing'
                
                new_columns[new_col_name] = 'categorical'
                print(f"Binned: {col} -> {new_col_name}")
            except Exception as e:
                print(f"Skipping {col} due to error: {e}")
                
    # Add new column types
    final_headers = headers + list(new_columns.keys())
    final_types = types + list(new_columns.values())
    
    # Reconstruct lines
    with open(output_file, 'w') as f:
        f.write('\t'.join(final_headers) + '\n')
        f.write('\t'.join(final_types) + '\n')
        
        for index, row in df.iterrows():
            row_data = [str(row[col]) if pd.notna(row[col]) else "" for col in final_headers]
            f.write('\t'.join(row_data) + '\n')
            
    print(f"\nSuccess! Binned metadata saved to: {output_file}")

if __name__ == "__main__":
    main()
