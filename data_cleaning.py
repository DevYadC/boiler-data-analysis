import pandas as pd
import numpy as np
from data_cleaning_utils import clean_cell_data, extract_units, check_threshold

# read  CSV file
df = pd.read_csv('CEC_Data_2021.csv')

# clean and convert timestamp; remove ' Rel' and set as datetime index
df['Timestamp'] = pd.to_datetime(df['Timestamp'].str.replace(' Rel', ''), utc=True)
df.set_index('Timestamp', inplace=True)

# rename columns to include units, from cells
for column in df.columns:
    first_value = df[column].dropna().iloc[0]
    unit = extract_units(first_value)
    new_column_name = f"{column} ({unit})"
    df.rename(columns={column: new_column_name}, inplace=True)

# remove units from cells, leave only numerical data
df = df.applymap(clean_cell_data)


# forward-fill missing data with a limit of 2 and drop remaining NaN rows
df = df.fillna(method='ffill', limit=2)
df = df.dropna()

# apple check threshold function to each cell, removing outliers
for column in df.columns:
    df[column] = df.apply(lambda row: check_threshold(row[column], column), axis=1)

# drop rows that have NaN values
df = df.dropna()

# convert temperature columns from Celsius to Kelvin
temperature_columns = [col for col in df.columns if '°C' in col]
for col in temperature_columns:
    df[col] = df[col] + 273.15
    new_col_name = col.replace('°C', 'K')
    df.rename(columns={col: new_col_name}, inplace=True)

# save data for boiler B-2, as well as temperature and humidity outside plant
boiler_2_columns = [col for col in df.columns if "B-2" in col or col in ['UBC Temp (K)', 'UBC Humidity (%RH)']]
boiler_2_df = df[boiler_2_columns].dropna()
boiler_2_df.to_csv('Clean_CEC_Data_B2_2021.csv')
