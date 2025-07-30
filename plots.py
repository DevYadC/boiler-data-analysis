
import pandas as pd
import matplotlib.pyplot as plt

# Load all precomputed data
# Assumes 'Timestamp' column exists and is parseable
df = pd.read_csv(
    'B2_Cleaned_2021_MassBalance_Efficiency.csv',
    index_col='Timestamp',
    parse_dates=True
)

# Simple IQR filter
def iqr_filter(s, k=1.5):
    q1, q3 = s.quantile([0.25, 0.75])
    return s[s.between(q1 - k*(q3 - q1), q3 + k*(q3 - q1))]

# 1. Flue gas flow & UBC Temp
flow = iqr_filter(df['flue gas molar flowrate (mol/h)'])
plt.figure()
plt.plot(flow.index, flow, label='Flue Gas Flow (mol/h)')
plt.plot(df.index, df['UBC Temp (K)'], label='UBC Temp (K)')
plt.legend()
plt.show()

# 2. Flue gas composition
plt.figure()
for col in [
    'mole fraction O2 in flue gas',
    'mole fraction N2 in flue gas',
    'mole fraction H2O in flue gas',
    'mole fraction CO2 in flue gas'
]:
    plt.plot(df.index, df[col].clip(lower=0), label=col)
plt.legend()
plt.show()

# 3. Boiler efficiency comparison
plt.figure()
plt.plot(df.index, df['boiler_efficiency_percent'], label='Calculated η')
plt.plot(df.index, df['Campus Energy Centre Boiler B-2 Efficiency (%)'], label='Sensor η')
plt.legend()
plt.show()

# 4. Boiler energy output (MWh)
plt.figure()
plt.plot(df.index, df['Campus Energy Centre Boiler B-2 Energy (MWh)'],
         label='Energy (MWh)')
plt.legend()
plt.show()

# 5. Exhaust NOx (ppm)
plt.figure()
plt.scatter(df.index,
            df['Campus Energy Centre Boiler B-2 Exhaust NOx (ppm)'],
            s=5)
plt.title('Exhaust NOx (ppm)')
plt.show()
