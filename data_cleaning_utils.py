import pandas as pd
import numpy as np

def is_superscript(char):
    superscript_chars = '⁰¹²³⁴⁵⁶⁷⁸⁹'
    return char in superscript_chars

def clean_cell_data(data):
    data_str = str(data)
    # Include only digits (ignoring superscripts) and specific characters
    cleaned_data = ''.join([i for i in data_str if (i.isdigit() and not is_superscript(i)) or i in ('.', '-', 'E')])
    return pd.to_numeric(cleaned_data, errors='coerce')

def extract_units(value):
    # Remove digits (except superscripts) and specific characters to extract unit string
    return ''.join([i for i in value if not((i.isdigit() and not is_superscript(i)) or i in ('.', '-', 'E', 'Â'))]).strip()

def check_threshold(cell, column):
    thresholds = {
        '°C': (-30, 300),
        '%': (0, 100),
        'ppm': (0, 500),
        'MWh': (0, 100),
        'L/s': (0, 500),
        'm³/h': (0, 1500),
        'MW': (0, 100),
        'kPa': (0, 300),
    }

    bounds = None
    for k, v in thresholds.items():
        if k in column:
            bounds = v
            break

    if bounds:
        lower, upper = bounds
        return cell if lower <= cell <= upper else pd.NA
    return cell
