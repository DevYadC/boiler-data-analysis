import pandas as pd
from data_cleaning_utils import clean_cell_data, extract_units, check_threshold

# 1. Read CSV
df = pd.read_csv('CEC_Data_2021.csv')

# 2. Parse timestamp and set as index
df['Timestamp'] = pd.to_datetime(
    df['Timestamp'].str.replace(' Rel', '', regex=False),
    utc=True
)
df.set_index('Timestamp', inplace=True)

# 3. Rename columns to include units (peek at first non‑null cell)
unit_map = {}
for col in df.columns:
    raw = df[col].dropna().iloc[0]
    unit = extract_units(raw)
    unit_map[col] = f"{col} ({unit})"
df.rename(columns=unit_map, inplace=True)

# 4. Strip unit text from every cell, leaving only numbers
df = df.applymap(clean_cell_data)

# 5. **Convert all Celsius→Kelvin** now, before filtering
c_cols = [col for col in df.columns if '(°C)' in col]
# add 273.15 and then bulk‑rename
for col in c_cols:
    df[col] = df[col].astype(float) + 273.15

df.rename(
    columns={col: col.replace('(°C)', '(K)') for col in c_cols},
    inplace=True
)

# 6. Forward‑fill up to 2 gaps, then drop any remaining NaNs
df = df.fillna(method='ffill', limit=2).dropna()

# 7. Apply threshold‑based outlier removal **per column**
for col in df.columns:
    df[col] = df[col].astype(float).apply(lambda x: check_threshold(x, col))
df = df.dropna()

# 8. Print total gas flow for Boiler B‑2
total_gas = df["Campus Energy Centre Boiler B-2 Gas Flow Rate (m³/h)"].sum()
print("Total Gas Flow Rate:", total_gas)

# 9. Select B‑2 columns + Temp(K) & Humidity, then save
boiler_2_cols = [
    c for c in df.columns
    if "B-2" in c or c in ["UBC Temp (K)", "UBC Humidity (%RH)"]
]
boiler_2_df = df[boiler_2_cols].dropna()
boiler_2_df.to_csv('B2_Cleaned_CEC_Data_2021.csv', index=True)
