import pandas as pd
import numpy as np
from scipy.integrate import quad

# ───── Constants ─────
# Water Cp polynomial coefficients (J/mol·K)
A_H2O, B_H2O, C_H2O, D_H2O = 7.243e1, 1.039e-2, -1.497e-6, 0.0
# Standard heats of formation (J/mol)
ΔHf_CO2  = -393.5e3
ΔHf_H2O  = -248.1e3
ΔHf_CH4  =  -74.6e3
ΔHf_C2H6 =  -83.75e3
ΔHf_O2   =      0.0

# Assumed fuel composition
FRAC_CH4, FRAC_C2H6 = 0.95, 0.05

# Reference temperature for ΔH integration (K)
T_REF = 298.15


# ───── Functions ─────
def read_data(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


def molar_flow_water(df: pd.DataFrame) -> pd.Series:
    """
    Convert water flow rate (L/s) → mol/s
    """
    W_L_s = df['Campus Energy Centre Boiler B-2 Water Flow Rate (L/s)']
    density_H2O = 1000.0   # g/L
    MW_H2O     =   18.0    # g/mol
    return W_L_s * density_H2O / MW_H2O


def integrate_cp(T_K: np.ndarray) -> np.ndarray:
    """
    Integrate Cp(T) = A + B·T + C·T² + D·T³
    from T_REF to each T_K to get ΔH (J/mol).
    """
    def cp(T, A, B, C, D):
        return A + B*T + C*T**2 + D*T**3

    # vectorized quad
    return np.array([
        quad(cp, T_REF, T, args=(A_H2O, B_H2O, C_H2O, D_H2O))[0]
        for T in T_K
    ])


def heat_of_combustion() -> float:
    """
    Weighted average heat of combustion (J/mol) of the CH4/C2H6 mix.
    """
    Hc_CH4  =  (ΔHf_CO2 + 2*ΔHf_H2O) - ΔHf_CH4    # CH4 → CO2 + 2H2O
    Hc_C2H6 = (4*ΔHf_CO2 + 6*ΔHf_H2O) - 2*ΔHf_C2H6 # 2C2H6 → 4CO2 + 6H2O
    return FRAC_CH4*Hc_CH4 + FRAC_C2H6*Hc_C2H6


def compute_efficiency(df: pd.DataFrame) -> pd.Series:
    """
    Compute boiler efficiency (%) row by row:
      η = Q_water_out / Q_fuel_in
        = [W·(ΔH_hw - ΔH_cw)] / [ṅ_fuel (mol/s) · H_combustion] × 100
    """
    # 1) water molar flow (mol/s)
    W = molar_flow_water(df)

    # 2) water enthalpy rises (J/mol)
    T_in_K  = df['Campus Energy Centre Boiler B-2 Entering Water Temp (°C)'] + 273.15
    T_out_K = df['Campus Energy Centre Boiler B-2 Leaving Water Temp (°C)']  + 273.15
    ΔH_cw = integrate_cp(T_in_K.values)
    ΔH_hw = integrate_cp(T_out_K.values)

    # 3) fuel molar flow (mol/s) — convert from the mass‑balance output (mol/h)
    n_fuel_h = df['fuel_molar_flow_mol_h']
    n_fuel_s = n_fuel_h / 3600.0

    # 4) fuel heating value (J/mol)
    H_comb = heat_of_combustion()

    # 5) efficiency
    η = W * (ΔH_hw - ΔH_cw) / (n_fuel_s * H_comb) * 100

    return η


def main():
    # a) load the mass‑balance‑enriched data
    df = read_data('B2_Cleaned_Data_2021_with_mass_balance.csv')

    # b) compute and append efficiency
    df['boiler_efficiency_%'] = compute_efficiency(df)

    # c) save to a new file
    df.to_csv(
      'B2_Cleaned_Data_2021_with_mass_balance_and_efficiency.csv',
      index=False
    )
    print("Wrote new file with boiler_efficiency_%")

if __name__ == '__main__':
    main()
