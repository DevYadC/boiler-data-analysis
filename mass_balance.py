import pandas as pd



def molar_flowrate_fuel(V_L_h: pd.Series, P_kPa: pd.Series, T_K: pd.Series) -> pd.Series:
    """
    calculate fuel molar flowrate in mol/h using ideal gas law:
      n = P*V / (R*T)
    where V is in L/h , P in kPa, T in K.
    """
    R_L_kPa_per_mol_K = 8.314        # L·kPa/(mol·K)
    fuel_n_mol_h = (P_kPa * V_L_h) / (R_L_kPa_per_mol_K * T_K)
    return fuel_n_mol_h.rename('fuel molar flowrate (mol/h)')
                   



def air_mole_fractions(T_K: pd.Series, RH_frac: pd.Series) -> pd.DataFrame:
    """
    computes inlet air mole fractions of H2O, O2, and N2, using temperature and relative humidity
    returns a DataFrame with columns y_H2O, y_O2, y_N2.
    """
    ANTOINE_A = 5.40221 # Antoine equation constants for water
    ANTOINE_B = 1838.675 
    ANTOINE_C = -31.737 
    P_ATM_Pa = 101_325 # atmospheric pressure in Pa

    # Saturation pressure of water (Pa) via Antoine
    P_sat_Pa = (10 ** (ANTOINE_A - ANTOINE_B / (ANTOINE_C + T_K))) * 1e5
    P_H2O = RH_frac * P_sat_Pa
    P_O2 = 0.21 * P_ATM_Pa
    P_N2 = 0.79 * P_ATM_Pa
    P_total = P_H2O + P_O2 + P_N2

    return pd.DataFrame({
        'mole fraction H2O in air': P_H2O / P_total,
        'mole fraction O2 in air':  P_O2 / P_total,
        'mole fraction N2 in air':  P_N2 / P_total,
    })




def flue_gas_balance(fuel_n_mol_h: pd.Series,
                     frac_CH4: float,
                     frac_C2H6: float,
                     y_air: pd.DataFrame,
                     fraction_O2_stack: pd.Series) -> pd.DataFrame:
    """
    Perform elemental balances to get:
      • A_in (air molar flow in mol/h)
      • G    (flue gas molar flow  mol/h)
      • y_i  (flue gas mole fractions)

    Reaction equations:
    Methane: CH4 + 2O2 -> CO2 + 2H20 
    Ethane: 2C2H6 + 7O2 -> 4CO2 + 6H2O
    """
    # split fuel stream
    n_CH4  = frac_CH4  * fuel_n_mol_h
    n_C2H6 = frac_C2H6 * fuel_n_mol_h

    # stoichiometry 
    O2_needed   = 2.0 * n_CH4 + 3.5 * n_C2H6
    CO2_prod    = 1.0 * n_CH4 + 2.0 * n_C2H6
    H2O_prod    = 2.0 * n_CH4 + 6.0 * n_C2H6

    # inlet air fractions
    y_H2O = y_air['mole fraction H2O in air']
    y_O2  = y_air['mole fraction O2 in air']
    y_N2  = y_air['mole fraction N2 in air']

    #derived by expanding fraction_O2_stack= molar_flow_rate_O2/ molar_flow_rate_flue_gas, and isolating for molar air flow rate in (A_in)
    A_in = (
        fraction_O2_stack * (O2_needed - CO2_prod - H2O_prod)
        - O2_needed
    ) / (fraction_O2_stack * (y_H2O + y_O2 + y_N2) - y_O2)

    # total flue gas flow G (mol/h)
    G = (
        (A_in * y_O2 - O2_needed)
      + (A_in * y_N2)
      + (A_in * y_H2O + H2O_prod)
      + CO2_prod
    )

    # individual component flows
    G_O2   = A_in * y_O2   - O2_needed
    G_N2   = A_in * y_N2
    G_H2O  = A_in * y_H2O  + H2O_prod
    G_CO2  = CO2_prod

    return pd.DataFrame({
        'air inlet molar flowrate (mol/h)':         A_in,
        'flue gas molar flowrate (mol/h)':            G,
        'mole fraction O2 in flue gas':          G_O2 / G,
        'mole fraction N2 in flue gas':          G_N2 / G,
        'mole fraction H2O in flue gas':         G_H2O / G,
        'mole fraction CO2 in flue gas':         G_CO2 / G,
    }, index=fuel_n_mol_h.index)





def main():
    df = pd.read_csv('B2_Cleaned_CEC_Data_2021.csv')

    # note: df['UBC Temp (°C)'] is already in Kelvin
    T_K = df['UBC Temp (K)']

    fuel_n = molar_flowrate_fuel(
        V_L_h = df['Campus Energy Centre Boiler B-2 Gas Flow Rate (m³/h)']*1e3,
        P_kPa  = df['Campus Energy Centre Boiler B-2 Gas Pressure (kPa)'],
        T_K    = T_K,
    )

    air_y = air_mole_fractions(
        T_K    = T_K,
        RH_frac= df['UBC Humidity (%RH)'] / 100,
    )

    flue = flue_gas_balance(
        fuel_n_mol_h       = fuel_n,
        frac_CH4           = 0.95,
        frac_C2H6          = 0.05,
        y_air              = air_y,
        fraction_O2_stack  = df['Campus Energy Centre Boiler B-2 Exhaust O2 (%)'] / 100,
    )

    # — join everything back onto the original df —
    results = pd.concat([fuel_n, air_y, flue], axis=1)
    df_with_results = pd.concat([df, results], axis=1)

    # — save the enriched DataFrame —
    df_with_results.to_csv(
        'B2_Cleaned_Data_2021_with_mass_balance.csv',
        index=False
    )

main()
