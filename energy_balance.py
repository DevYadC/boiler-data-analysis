import sympy as sp
import pandas as pd
# create dataframe with boiler 2 data
df = pd.read_csv('B2_Cleaned_CEC_Data_2021.csv')
b2_row=df.iloc[0] #first row of data for boiler 2

#variables
T_ref=298.15 #reference temperature in kelvin
water_flowrate=b2_row['Campus Energy Centre Boiler B-2 Water Flow Rate (L/s)'] #water flow rate (L/s)
cw_T_in=b2_row['Campus Energy Centre Boiler B-2 Entering Water Temp (°C)'] + 273.15 #cold water inlet temperature (K)
hw_T_out=b2_row['Campus Energy Centre Boiler B-2 Leaving Water Temp (°C)'] + 273.15 #hot water outlet temperature (K)



#return enthalpy in J/mol
def calculate_enthalpy(a,b,c,d,T_f, T_reference):
    T = sp.symbols('x')
    f = a + b*T + c*T**2 + d*T**3
    enthalpy = sp.integrate(f, (T, T_reference, T_f))
    return enthalpy

#enthalpy (J/mol) , molar_flow (mol/s)
#returns q in (J/s)
def calculate_heat_transfer(enthalpy, molar_flow):
    return enthalpy*molar_flow


#Determine heat transfer from cold water (q_cw) and hot water (q_hw)
##heat capacity coefficients water
a_H2O=7.243E1
b_H2O=1.039E-2
c_H2O=-1.497E-6
d_H2O=0
##convert water flow rate from L/s -> mol/s
density_H2O= 1000 # g/L
MW_H2O= 18 # molar weight water = 18g/mol
W = water_flowrate*density_H2O*MW_H2O #molar flow rate water (mol/s)

##calculate enthalpy cold water
H_cw=calculate_enthalpy(a_H2O, b_H2O, c_H2O, d_H2O, cw_T_in, T_ref)
print(f'cold water inlet stream heat capacity J/mol: {H_cw}')

##calculate enthalpy hot water
H_hw=calculate_enthalpy(a_H2O, b_H2O, c_H2O, d_H2O, hw_T_out, T_ref)
print(f'hot water inlet stream heat capacity J/mol: {H_hw}')

##calculate heat transfer with eqn: q_cw= W * Cp_cw * (T_cw – T_o) = W * H_cw
q_cw= W * H_cw
q_hw= W * H_hw
print(f'Heat transfer from cold water inlet stream(J/s): {q_cw}')
print(f'Heat transfer from hot water outlet stream(J/s): {q_hw}')






#Determine heat transfer from fuel stream (q_f) , fuel composition: 95% CH4, 5% C2H6

##heat capacity coefficients
a_CH4=19.25
b_CH4=0.05213
c_CH4=1.197E-5
d_CH4=-1.132E-8

a_C2H6=5.409
b_C2H6=0.1781
c_C2H6=-6.938E-5
d_C2H6=8.713E-9

##mole fractions in fuel
y_CH4=0.95
y_C2H6=0.05

#assume fuel temperature is ambient temperature
T_fuel= b2_row['UBC Temp (°C)'] +273.15 #Kelvins

#calculate enthalpy of fuel with equation: H_fuel = y_CH4 * H_CH4 + y_C2H6 * H_C2H6
H_CH4 = calculate_enthalpy(a_CH4, b_CH4, c_CH4, d_CH4, T_fuel, T_ref)
H_C2H6 = calculate_enthalpy(a_C2H6, b_C2H6, c_C2H6, d_C2H6, T_fuel, T_ref)

H_fuel = y_CH4*H_CH4 + y_C2H6*H_C2H6
print(f'Enthalpy of fuel mixture J/mol: {H_fuel}')






#Determine heat of combustion of fuel
deltaHf_CO2 = -393.5E3 # value for ΔHf° of CO2 J/mol
deltaHf_H2O = -248.1E3 # value for ΔHf° of H2O vapor J/mol
deltaHf_CH4 = -74.6E3# value for ΔHf° of CH4 J/mol
deltaHf_C2H6 = -83.75# value for ΔHf° of C2H6 J/mol
deltaHf_O2 = 0 # ΔHf° of O2 is zero since it's an element in its standard state 

# derived from combustion reactions of CH4 and C2H6
def heat_of_combustion_CH4():
    return deltaHf_CO2 + 2 * deltaHf_H2O - deltaHf_CH4 - 2 * deltaHf_O2

def heat_of_combustion_C2H6():
    return 4 * deltaHf_CO2 + 6 * deltaHf_H2O - 2 * deltaHf_C2H6 - 7 * deltaHf_O2

H_combustion_fuel= y_CH4 * heat_of_combustion_CH4() + y_C2H6 * heat_of_combustion_C2H6()

print(f'Heat of combustion of fuel J/mol: {H_combustion_fuel}')





#Determine heat transfer of fuel, use equation: q_F=F*{Cp_F*(T_F – T_o) + H_cF} = F * {H_F + H_cF }

#calculate fuel molar flow rate (mol/s)
fuel_volumetric_flowrate= b2_row['Campus Energy Centre Boiler B-2 Gas Flow Rate (m³/h)'] * (1000/1) * 3600 # L/s
fuel_pressure = b2_row['Campus Energy Centre Boiler B-2 Gas Pressure (kPa)'] 
R= 8.314 #L*kPa/mol*K

#use P*V=n*R*T to calculate molar flowrate
def molar_flowrate(P,V,R,T):
    n=(P*V)/(R*T)
    return n

F = molar_flowrate(fuel_pressure, fuel_volumetric_flowrate, R, T_fuel) #fuel molar flow rate (mol/s)

qF = F * (H_fuel + H_combustion_fuel)
print(f'Heat transfer of fuel (Js/): {qF}')





#Determine enthalpy of dry air, using eqn: H_AD = y_O2 * H_O2 + y_N2 * H_N2

y_O2=0.21 #mole fraction N2 in dry air
y_N2=0.79 #mole fraction O2 in dry air


# heat capacity coefficients for O2 and N2
a_O2 = 28.11
b_O2 = -3.7E-6
c_O2 = 1.746E-5
d_O2 = 1.065E-8

a_N2 = 31.15
b_N2 = -0.01357
c_N2 = 2.68E-6
d_N2 = -1.168E-8

T_air= b2_row['UBC Temp (°C)'] +273.15 #air temperature (K)

H_O2 = calculate_enthalpy(a_O2, b_O2, c_O2, d_O2, T_air, T_ref)

H_N2 = calculate_enthalpy(a_N2, b_N2, c_N2, d_N2, T_air, T_ref)

H_AD = y_O2 * H_O2 + y_N2 * H_N2


#Determine enthalpy of water vapour in air

a_H2O_vapor=32.24
b_H2O_vapor=0.001924
c_H2O_vapor=1.055E-5
d_H2O_vapor=-3.596E-9

H_H2O_gas = calculate_enthalpy(a_H2O_vapor, b_H2O_vapor, c_H2O_vapor, d_H2O_vapor, T_air, T_ref)