import pandas as pd
import numpy as np
from scipy.integrate import quad

# ----- Constants -----
CP_COEFF_A_H2O = 7.243e1      # J/mol·K
CP_COEFF_B_H2O = 1.039e-2
CP_COEFF_C_H2O = -1.497e-6
CP_COEFF_D_H2O = 0.0

ENTHALPY_FORMATION_CO2 = -393.5e3   # J/mol
ENTHALPY_FORMATION_H2O = -248.1e3   
ENTHALPY_FORMATION_CH4 =  -74.6e3   
ENTHALPY_FORMATION_C2H6 =  -83.75e3
ENTHALPY_FORMATION_O2  =    0.0     

FUEL_FRACTION_CH4   = 0.95
FUEL_FRACTION_C2H6  = 0.05

REFERENCE_TEMPERATURE_K = 298.15     # K

# ----- Thermo Calculations -----

def integrate_cp_water(T_k: np.ndarray, T_ref: float = REFERENCE_TEMPERATURE_K) -> np.ndarray:
    """
    Integrate Cp(T) = A + B·T + C·T^2 + D·T^3
    from reference temperature to each final temperature.
    Returns delta enthalpy (J/mol).
    """
    def cp_function(T, a, b, c, d):
        return a + b*T + c*T**2 + d*T**3

    return np.array([
        quad(cp_function, T_ref, T, 
             args=(CP_COEFF_A_H2O, CP_COEFF_B_H2O, CP_COEFF_C_H2O, CP_COEFF_D_H2O))[0]
        for T in T_k
    ])

def heat_of_combustion_per_mol() -> float:
    """
    Compute heat of combustion of the CH4/C2H6 mixture (J/mol).
    """
    # CH4 + 2 O2 -> CO2 + 2 H2O
    heat_ch4  = ((ENTHALPY_FORMATION_CO2 + 2*ENTHALPY_FORMATION_H2O)
                 - (ENTHALPY_FORMATION_CH4 + 2*ENTHALPY_FORMATION_O2))
    # 2 C2H6 + 7 O2 -> 4 CO2 + 6 H2O
    heat_c2h6 = ((4*ENTHALPY_FORMATION_CO2 + 6*ENTHALPY_FORMATION_H2O)
                 - (2*ENTHALPY_FORMATION_C2H6 + 7*ENTHALPY_FORMATION_O2))

    return abs(FUEL_FRACTION_CH4 * heat_ch4 + FUEL_FRACTION_C2H6 * heat_c2h6)

# ----- Efficiency Calculation -----
def compute_boiler_efficiency(df: pd.DataFrame) -> pd.Series:
    """
    Calculate boiler efficiency (%) for each row:
      efficiency = Q_water / Q_fuel * 100
    where:
      Q_water = molar_flow_water * (enthalpy_hot - enthalpy_cold)
      Q_fuel  = fuel_molar_flow_mol_s * heat_of_combustion
    """
    # 1) water molar flow (mol/s)
    water_molar_flow = df['Campus Energy Centre Boiler B-2 Water Flow Rate (L/s)']*1000.0 / 18.0

    # 2) enthalpy change of water streams (J/mol)
    inlet_K  = df['Campus Energy Centre Boiler B-2 Entering Water Temp (K)'] 
    outlet_K = df['Campus Energy Centre Boiler B-2 Leaving Water Temp (K)']  
    enthalpy_cold = integrate_cp_water(inlet_K.values)
    enthalpy_hot  = integrate_cp_water(outlet_K.values)

    # 3) fuel molar flow in mol/s (input CSV gives mol/h)
    fuel_mol_s = df['fuel molar flowrate (mol/h)']/3600.0
    

    # 4) heat of combustion (J/mol)
    heating_value = heat_of_combustion_per_mol()

    # 5) compute efficiency
    efficiency = (
        water_molar_flow * (enthalpy_hot - enthalpy_cold)
        / (fuel_mol_s * heating_value)
    ) * 100.0

    return efficiency


def main():
    # load boiler data
    df = pd.read_csv('B2_Cleaned_Data_2021_with_mass_balance.csv')

    # compute and append efficiency as column
    df['boiler_efficiency_percent'] = compute_boiler_efficiency(df)

    # save results to new CSV
    df.to_csv(
        'B2_Cleaned_2021_MassBalance_Efficiency.csv',
        index=False
    )
    


main()
