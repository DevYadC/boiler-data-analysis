import pandas as pd
import numpy as np
import sympy as sp
# create dataframe with boiler 2 data
df = pd.read_csv('B2_Cleaned_CEC_Data_2021.csv')


b2_row=df.iloc[0] #first row of data for boiler 2



#*********Mass Balance********

#***Determine the molar flow rate of fuel at operating conditions***
#variables 
fuel_volumetric_flowrate=b2_row['Campus Energy Centre Boiler B-2 Gas Flow Rate (m³/h)'] #m³/h
fuel_pressure=b2_row['Campus Energy Centre Boiler B-2 Gas Pressure (kPa)'] #kPa
R= 8.314 #L*kPa/mol*K
T=b2_row['UBC Temp (°C)'] #°C
#use P*V=n*R*T to calculate molar flowrate
def molar_flowrate(P,V,R,T):
    n=(P*V)/(R*T)
    return n
#unit conversions
fuel_volumetric_flowrate=fuel_volumetric_flowrate*(1000/1) #m³/h -> L/h
T = T + 273.15 # °C -> K
#calculate molar flow rate (moles/h)
print(f'Molar flowrate fuel in (moles/h): {molar_flowrate(fuel_pressure, fuel_volumetric_flowrate, R, T)}')



#***Determine the mole fractions of wet air***
#variables 
rel_humidity=b2_row['UBC Humidity (%RH)']/100 #relative humidity as fraction
T=b2_row['UBC Temp (°C)'] #°C
P_atmosphere=101325# Pascals
A= 5.40221 #antoine A constant water (	Bridgeman and Aldrich, 1964	)
B= 1838.675 #antoine B constant water (	Bridgeman and Aldrich, 1964	)
C= -31.737 #antoine C constant water (	Bridgeman and Aldrich, 1964	)

# Convert T to Kelvin for consistency with P_total in Pa
T_K = T + 273.15

# Calculate saturation pressure using the Antoine equation
P_sat = 10 ** (A - (B / (C + T_K))) # P_sat in bar
P_sat_Pa = P_sat * 100000 # Convert bar to Pa

# Calculate partial pressure of water vapor using relative humidity
P_H2O = rel_humidity * P_sat_Pa

#Calculate partial pressures of dry air (79% N2 and 21% O2)
P_O2= 0.21*P_atmosphere
P_N2=0.79*P_atmosphere

#Calculate P total
P_total=P_H2O+P_O2+P_N2

# Calculate mole fractions
y_H2O_air = P_H2O / P_total
y_O2_air= P_O2/P_total
y_N2_air=P_N2/P_total

print(f'Air in mole fraction water: {y_H2O_air}, mole fraction oxygen: {y_O2_air}, mole fraction nitrogen: {y_N2_air}')




#***Determine the mole fractions of flue gas using atomic carbon, nitrogen, oxygen, and hydrogen balances***
# Given values
fuel_molar_flowrate=molar_flowrate(fuel_pressure, fuel_volumetric_flowrate, R, T) #(mol/h)
percent_CH4 = 0.95
percent_C2H6 = 0.05

# Calculate the molar flow rates of CH4 and C2H6
flow_rate_CH4 = percent_CH4 * fuel_molar_flowrate  # mol/h
flow_rate_C2H6 = percent_C2H6 * fuel_molar_flowrate  # mol/h

#Methane: CH4 + 2O2 -> CO2 + 2H20 
#Ethane: 2C2H6 + 7O2 -> 4CO2 + 6H2O

# Combustion equations coefficients for CH4 and C2H6
O2_coeff_CH4 = 2
O2_coeff_C2H6 = 3.5
CO2_coeff_CH4 = 1
CO2_coeff_C2H6 = 2
H2O_coeff_CH4 = 2
H2O_coeff_C2H6 = 6

# O2 consumed, CO2 and H2O produced
O2_consumed = O2_coeff_CH4 * flow_rate_CH4 + O2_coeff_C2H6 * flow_rate_C2H6
CO2_produced = CO2_coeff_CH4 * flow_rate_CH4 + CO2_coeff_C2H6 * flow_rate_C2H6
H2O_produced = H2O_coeff_CH4 * flow_rate_CH4 + H2O_coeff_C2H6 * flow_rate_C2H6

# Define the percentage of O2 in the stack
fraction_O2_stack = b2_row['Campus Energy Centre Boiler B-2 Exhaust O2 (%)']/100 #from data set for boiler 2

#Note: y_H2O , y_N2, _O2 are from previous calculations and represent mole fractions of inlet air stream

#derived by expanding fraction_O2_stack= molar_flow_rate_O2/ molar_flow_rate_flue_gas, and isolating for molar air flow rate in (A_in)
A_in = (fraction_O2_stack * (O2_consumed - H2O_produced - CO2_produced) - O2_consumed) / (fraction_O2_stack * (y_H2O_air + y_N2_air + y_O2_air) - y_O2_air)
print(f'Molar flow rate Air in (mol/h): {A_in}')

#Determine molar flow rate of flue gas (N_flue (mol/h))
N_flue= (A_in*y_O2_air-O2_consumed)+(A_in*y_N2_air)+(A_in*y_H2O_air+H2O_produced)+CO2_produced
print(f'Molar flowrate flue gas out(mol/h)')
#determine molar flow rates of components in flue gas (mol/h)
N_O2_flue= A_in*y_O2_air-O2_consumed
N_N2_flue=A_in*y_N2_air
N_H2O_flue=A_in*y_H2O_air+H2O_produced
N_CO2_flue=CO2_produced

#determine mole fractions of components in flue glass
y_O2_flue= N_O2_flue/N_flue
y_N2_flue= N_N2_flue/N_flue
y_H2O_flue= N_H2O_flue/N_flue
y_CO2_flue= N_CO2_flue/N_flue

print(f'Flue gas mole fraction O2, N2, H2O, CO2: {y_O2_flue, y_N2_flue, y_H2O_flue,y_CO2_flue }')








#**********ENERGY BALANCE**********

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
print(f'Heat transfer of fuel (J/s): {qF}')





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


#Determine enthalpy of water vapour in air stream

a_H2O_vapor=32.24
b_H2O_vapor=0.001924
c_H2O_vapor=1.055E-5
d_H2O_vapor=-3.596E-9

H_AH2O = calculate_enthalpy(a_H2O_vapor, b_H2O_vapor, c_H2O_vapor, d_H2O_vapor, T_air, T_ref) #enthalpy of water vapor in inlet air stream

#Determine heat transfer of dry air using eqn: q_AD = A_D * H_AD
# A_D = molar flow rate of dry air (mol/s)
# H_AD = enthalpy dry air (J/mol)
A_D = y_O2 * A_in + y_N2 * A_in
q_AD = A_in * H_AD


#Determine heat transfer of water vapour in air stream using eqn: q_AH2O = A_H2O * ( H_AH2O + H_vap )
# A_H2O = molar flow rate of water vapour in air stream (mol/s)
# H_AH2O = #enthalpy of water vapor in inlet air stream (J/mol)
H_vap = 40656  #Heat of vaporization of water (J/mol) at reference temperature
A_H2O = y_H2O_air * A_in
q_AH2O = A_H2O * ( H_AH2O + H_vap )

#Determine heat transfer of Air stream using eqn: q_A = q_AD + q_AH2O
q_A = q_AD + q_AH2O
print(f'heat transfer of Air stream (J/mol): {q_A}')







#Determine enthalpy of dry flue gas (ie excluding water vapor)
#use equation: H_flue_dry = y_CO2_flue * H_CO2 + y_N2_flue * H_N2 + y_O2_flue * H_O2

#Heat capacity coefficients 
a_CO2 = 19.8 
b_CO2 = 0.07344
c_CO2 = -5.602E-5
d_CO2 = 1.715E-8

#mole fractions of components excluding H2O (note: this is not same as overall mole fractions of flue gas) 
y_O2_flue_dry = y_O2_flue/ ( y_O2_flue + y_N2_flue + y_CO2_flue )
y_N2_flue_dry = y_N2_flue/ ( y_O2_flue + y_N2_flue + y_CO2_flue )
y_CO2_flue_dry = y_CO2_flue/ ( y_O2_flue + y_N2_flue + y_CO2_flue )

#mole fraction of water vapor in flue gas


T_flue = b2_row['Campus Energy Centre Boiler B-2 Exhaust Temp (°C)'] + 273.15 #flue gas temp (K)
H_CO2 = calculate_enthalpy(a_CO2, b_CO2, c_CO2, d_CO2, T_flue, T_ref)
H_flue_dry = y_CO2_flue_dry * H_CO2 + y_N2_flue_dry * H_N2 + y_O2_flue_dry * H_O2


#Determine enthalpy of wet flue gas (ie water vapor in flue gas)
H_flue_H2O = calculate_enthalpy(a_H2O_vapor, b_H2O_vapor, c_H2O_vapor, d_H2O_vapor, T_flue, T_ref)

#Determine enthalpy of flue gas
H_flue = ( y_O2_flue + y_N2_flue + y_CO2_flue ) * H_flue_dry  + y_H2O_flue * (H_flue_H2O)
print(f'Enthalpy of flue gas (J/mol): {H_flue}')