import pandas as pd
import numpy as np


file_path = 'Removed_Units_CEC_Data_2021.csv' 
data = pd.read_csv(file_path)

# Define outlier checking function for data
def check_threshold(cell, column):
    # Define thresholds based on column types
    if 'UBC Temp' in column:
        lower_bound, upper_bound = -40, 40
    elif '°C' in column and 'UBC Temp' not in column:
        lower_bound, upper_bound = 0, 300
    elif "%" in column:
        lower_bound, upper_bound = 0, 100
    elif 'ppm' in column:
        lower_bound, upper_bound = 0, 500
    elif 'MWh' in column:
        lower_bound, upper_bound = 0, 100
    elif 'L/s' in column:
        lower_bound, upper_bound = 0, 500
    elif 'm³/h' in column:
        lower_bound, upper_bound = 0, 2000
    elif 'MW' in column:
        lower_bound, upper_bound = 0, 100
    elif 'kPa' in column:
        lower_bound, upper_bound = 0, 300
    else:
        return False  
    return cell < lower_bound or cell > upper_bound


boiler_columns = {
    'Boiler 1': [],
    'Boiler 2': [],
    'Boiler 3': []
}


for col_name in data.columns:
    if 'Boiler B-1' in col_name:
        boiler_columns['Boiler 1'].append(col_name)
    elif 'Boiler B-2' in col_name:
        boiler_columns['Boiler 2'].append(col_name)
    elif 'Boiler B-3' in col_name:
        boiler_columns['Boiler 3'].append(col_name)


missing_values = {}
outliers_counts = {}

for boiler, cols in boiler_columns.items():
    missing_count = 0
    outlier_count = 0
    
    for col in cols:
        # Count missing values
        missing_count += data[col].isna().sum()
        
        # Count outliers by applying the check_threshold function row by row
        for value in data[col]:
            if pd.isna(value):  # Skip missing values
                continue
            if check_threshold(value, col):
                outlier_count += 1
                
    missing_values[boiler] = missing_count
    outliers_counts[boiler] = outlier_count

print("Missing Values:", missing_values)
print("Outliers Counts:", outliers_counts)
