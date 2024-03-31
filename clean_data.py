#Author: Yadpreet Cheema
#script to clean data set from UBC energy centre

import pandas as pd
import numpy as np

##Helper Functions
#checks if a character is a superscript
def is_superscript(char):
    superscript_chars = '⁰¹²³⁴⁵⁶⁷⁸⁹'
    return char in superscript_chars

#returns data only containing digits, '.', '-', or 'E'
def clean_cell_data(data):
    data_str = str(data)
    # Ensure that only digits, '.', '-', 'E' are included
    cleaned_data = ''.join([i for i in data_str if i.isdigit() and not is_superscript(i) or i in ('.', '-', 'E')])
    # Convert the cleaned string to numeric, handling errors by coercing to NaN
    return pd.to_numeric(cleaned_data, errors='coerce')

#returns string of units for data in given column
def extract_units(value):
    return ''.join([i for i in value  if not(i.isdigit() and not is_superscript(i) or i=='.' or i=='-' or i=="E" or i=="Â")]).strip()



# read the CSV file
df = pd.read_csv('CEC_Data_2021.csv')

# remove the ' Rel' from the timestamp strings
df['Timestamp'] = df['Timestamp'].str.replace(' Rel', '')

# convert to datetime
df['Timestamp'] = pd.to_datetime(df['Timestamp'], utc=True)

#takes the 'Timestamp' column of the DataFrame df and makes it the new index (row labels) of the DataFrame
df.set_index('Timestamp', inplace=True)

# rename columns to include units
for column in df.columns:
    # extract the first non-null value in the column to use it to determine units for data in column
    first_value = df[column].dropna().iloc[0]
    
    #store units for column in variable
    unit=extract_units(first_value)

    # create a new column name that includes the unit string 
    new_column_name = f"{column} ({unit})"
    
    # Rename the column in the dataframe with units included
    df.rename(columns={column: new_column_name}, inplace=True)


# Apply the clean_cell_data function to every element in the DataFrame to remove units from data cells
df = df.applymap(clean_cell_data)

# Save the dataframe with remove units to be used to determine which dataframe is best
#df.to_csv('Removed_Units_CEC_Data_2021.csv')
df.to_csv('CEC_Data_2021_cleaned.csv')
#Forward fill with a limit of 2
df = df.fillna(method='ffill', limit=2)

#delete rows that still have nan (ie. columns where 3 or more consecutive missing data)
df = df.dropna()


def check_threshold(cell, column):
    if 'UBC Temp' in column: #threshold for ambient ubc temperature
        lower_bound = -40
        upper_bound = 40
    if '°C' in column and 'UBC Temp' not in column: #thresholds for temps part of process and not ambient ubc temperature
        lower_bound = 0
        upper_bound = 300
    elif "%" in column:
        lower_bound = 0
        upper_bound = 100
    elif 'ppm' in column:
        lower_bound = 0
        upper_bound = 500
    elif 'MWh' in column:
        lower_bound = 0
        upper_bound = 100
    elif 'L/s' in column:
        lower_bound = 0
        upper_bound = 500
    elif 'm³/h' in column:
        lower_bound = 0
        upper_bound = 2000
    elif 'MW' in column:
        lower_bound = 0
        upper_bound = 100
    elif 'kPa' in column:
        lower_bound = 0
        upper_bound = 300
    else:
        return cell
    if cell < lower_bound or cell > upper_bound:
        return np.nan
    else:
        return cell 

# Apply check_threshold element-wise
for column in df.columns:
    df[column] = df.apply(lambda row: check_threshold(row[column], column), axis=1)

# Now, drop any rows that have NaN values (indicating a threshold violation)
df = df.dropna()

total_gas_flow_rate = df["Campus Energy Centre Boiler B-2 Gas Flow Rate (m³/h)"].sum()

print("Total Gas Flow Rate:", total_gas_flow_rate)

print(df)
# Save the cleaned dataframe to a new CSV file
#df.to_csv('Cleaned_CEC_Data_2021.csv')

boiler_2_columns = [column for column in df.columns if "B-2" in column or column == 'UBC Temp (°C)' or column == 'UBC Humidity (%RH)']
boiler_2_df = df[boiler_2_columns].dropna()
#boiler_2_df.to_csv('B2_Cleaned_CEC_Data_2021.csv')