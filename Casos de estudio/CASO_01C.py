# -*- coding: utf-8 -*-
"""
FTES - Casos de analisis TFM
CASO 01C
@author: gumer freire
"""

import FTES
import numpy as np

def windload(z, years=50):
    """
    Function to get the wind load for urban area zone IV depending on the height (h)
    of a given point. Returns wind load per m2.
    Pressure coefficient cp given for facade zone D.
    Input parameters:
        z - height above ground of the point subjected to wind load
        years - return period for probability of wind value. Default = 50 years
    This calculation is introduced as Cprob coefficient following Eurocode EN 1991-1-4
    """
    
    #Probability coefficient for period of time considered. According to EN 1991-1-4
    K = 0.2 #Recommended value National Annex EN 1991-1-4
    n = 0.5 #Recommended value National Annex EN 1991-1-4
    p = 1/years
    if p > 0.99: p = 0.99 
    Cprob = ((1 - K * (-np.log(1 - p))) / (1 - K * (-np.log(0.98))))**n
    
    vb =  29 #basic speed of wind, according to CTE 29 m/s for zone C
    d_density = 1.25 #kg/m3
    qb = 0.5 * d_density * ((vb*Cprob)**2)
    
    #CTE parameters, Table D.2 for urban zone IV
    k = 0.22
    L = 0.3
    Z = 5.0
    if z > Z: maxz = z
    else: maxz = Z
    F = k * np.log(maxz/L)
    ce = F * (F+7*k) #exposure coefficient
    cp = 0.8
    
    return qb*ce*cp # wind load in N/m2
    

# #%% ELU - Ultimate limit state

# load_area = 4.8077
# glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
# G_coef = 1.35 # coeficiente de acciones permanentes
# Q_coef = 1.5 # coeficiente de acciones variables



# C01C_ELU = FTES.FTES_Structure()
# C01C_ELU.import_DxfStructure('CASO_01C/CASO_01C.dxf')

# #Set mechanical properties
# C01C_ELU.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
# C01C_ELU.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
# C01C_ELU.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
# C01C_ELU.set_Cable(Memberfilter='Allmembers')
# C01C_ELU.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# # Set restraints
# C01C_ELU.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
# #Set glazing load
# C01C_ELU.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

# #Set wind load
# C01C_ELU.addForce(load_area * windload(1.5625 * 1) * Q_coef, ['y'], Nodefilter='Windload_01')
# C01C_ELU.addForce(load_area * windload(1.5625 * 2) * Q_coef, ['y'], Nodefilter='Windload_02')
# C01C_ELU.addForce(load_area * windload(1.5625 * 3) * Q_coef, ['y'], Nodefilter='Windload_03')
# C01C_ELU.addForce(load_area * windload(1.5625 * 4) * Q_coef, ['y'], Nodefilter='Windload_04')
# C01C_ELU.addForce(load_area * windload(1.5625 * 5) * Q_coef, ['y'], Nodefilter='Windload_05')
# C01C_ELU.addForce(load_area * windload(1.5625 * 6) * Q_coef, ['y'], Nodefilter='Windload_06')
# C01C_ELU.addForce(load_area * windload(1.5625 * 7) * Q_coef, ['y'], Nodefilter='Windload_07')
# C01C_ELU.addForce(load_area * windload(1.5625 * 8) * Q_coef, ['y'], Nodefilter='Windload_08')
# C01C_ELU.addForce(load_area * windload(1.5625 * 9) * Q_coef, ['y'], Nodefilter='Windload_09')
# C01C_ELU.addForce(load_area * windload(1.5625 * 10) * Q_coef, ['y'], Nodefilter='Windload_10')
# C01C_ELU.addForce(load_area * windload(1.5625 * 11) * Q_coef, ['y'], Nodefilter='Windload_11')
# C01C_ELU.addForce(load_area * windload(1.5625 * 12) * Q_coef, ['y'], Nodefilter='Windload_12')
# C01C_ELU.addForce(load_area * windload(1.5625 * 13) * Q_coef, ['y'], Nodefilter='Windload_13')
# C01C_ELU.addForce(load_area * windload(1.5625 * 14) * Q_coef, ['y'], Nodefilter='Windload_14')
# C01C_ELU.addForce(load_area * windload(1.5625 * 15) * Q_coef, ['y'], Nodefilter='Windload_15')

# #Structure analysis
# C01C_ELU.set_solverParameters(4000, 50, 5000)
# C01C_ELU.calculate_structure()

# #Plot and export data
# C01C_ELU.plotStructure()
# C01C_ELU.export_Data('CASO_01C/CASO_01C_ELU.xlsx')
# C01C_ELU.export_DxfStructure('CASO_01C/CASO_01C_ELU.dxf')

# #%% ELS - Serviceability limit state

# load_area = 4.8077
# glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
# G_coef = 1.0 # coeficiente de acciones permanentes ELS
# Q_coef = 1.0 # coeficiente de acciones variables ELS

# C01C_ELS = FTES.FTES_Structure()
# C01C_ELS.import_DxfStructure('CASO_01C/CASO_01C.dxf')

# #Set mechanical properties
# C01C_ELS.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
# C01C_ELS.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
# C01C_ELS.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
# C01C_ELS.set_Cable(Memberfilter='Allmembers')
# C01C_ELS.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# # Set restraints
# C01C_ELS.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
# #Set glazing load
# C01C_ELS.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

# #Set wind load
# C01C_ELS.addForce(load_area * windload(1.5625 * 1) * Q_coef, ['y'], Nodefilter='Windload_01')
# C01C_ELS.addForce(load_area * windload(1.5625 * 2) * Q_coef, ['y'], Nodefilter='Windload_02')
# C01C_ELS.addForce(load_area * windload(1.5625 * 3) * Q_coef, ['y'], Nodefilter='Windload_03')
# C01C_ELS.addForce(load_area * windload(1.5625 * 4) * Q_coef, ['y'], Nodefilter='Windload_04')
# C01C_ELS.addForce(load_area * windload(1.5625 * 5) * Q_coef, ['y'], Nodefilter='Windload_05')
# C01C_ELS.addForce(load_area * windload(1.5625 * 6) * Q_coef, ['y'], Nodefilter='Windload_06')
# C01C_ELS.addForce(load_area * windload(1.5625 * 7) * Q_coef, ['y'], Nodefilter='Windload_07')
# C01C_ELS.addForce(load_area * windload(1.5625 * 8) * Q_coef, ['y'], Nodefilter='Windload_08')
# C01C_ELS.addForce(load_area * windload(1.5625 * 9) * Q_coef, ['y'], Nodefilter='Windload_09')
# C01C_ELS.addForce(load_area * windload(1.5625 * 10) * Q_coef, ['y'], Nodefilter='Windload_10')
# C01C_ELS.addForce(load_area * windload(1.5625 * 11) * Q_coef, ['y'], Nodefilter='Windload_11')
# C01C_ELS.addForce(load_area * windload(1.5625 * 12) * Q_coef, ['y'], Nodefilter='Windload_12')
# C01C_ELS.addForce(load_area * windload(1.5625 * 13) * Q_coef, ['y'], Nodefilter='Windload_13')
# C01C_ELS.addForce(load_area * windload(1.5625 * 14) * Q_coef, ['y'], Nodefilter='Windload_14')
# C01C_ELS.addForce(load_area * windload(1.5625 * 15) * Q_coef, ['y'], Nodefilter='Windload_15')

# #Structure analysis
# C01C_ELS.set_solverParameters(2000, 50, 5000)
# C01C_ELS.calculate_structure()

# #Plot and export data
# C01C_ELS.plotStructure()
# C01C_ELS.export_Data('CASO_01C/CASO_01C_ELS.xlsx')
# C01C_ELS.export_DxfStructure('CASO_01C/CASO_01C_ELS.dxf')

#%% ELS - Serviceability limit state WIND 1 YEAR

load_area = 4.8077
glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
G_coef = 1.0 # coeficiente de acciones permanentes ELS
Q_coef = 1.0 # coeficiente de acciones variables ELS

C01C_ELS_1YEAR = FTES.FTES_Structure()
C01C_ELS_1YEAR.import_DxfStructure('CASO_01C/CASO_01C.dxf')

#Set mechanical properties
C01C_ELS_1YEAR.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
C01C_ELS_1YEAR.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
C01C_ELS_1YEAR.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
C01C_ELS_1YEAR.set_Cable(Memberfilter='Allmembers')
C01C_ELS_1YEAR.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# Set restraints
C01C_ELS_1YEAR.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
#Set glazing load
C01C_ELS_1YEAR.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

#Set wind load
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 1, years=1) * Q_coef, ['y'], Nodefilter='Windload_01')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 2, years=1) * Q_coef, ['y'], Nodefilter='Windload_02')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 3, years=1) * Q_coef, ['y'], Nodefilter='Windload_03')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 4, years=1) * Q_coef, ['y'], Nodefilter='Windload_04')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 5, years=1) * Q_coef, ['y'], Nodefilter='Windload_05')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 6, years=1) * Q_coef, ['y'], Nodefilter='Windload_06')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 7, years=1) * Q_coef, ['y'], Nodefilter='Windload_07')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 8, years=1) * Q_coef, ['y'], Nodefilter='Windload_08')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 9, years=1) * Q_coef, ['y'], Nodefilter='Windload_09')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 10, years=1) * Q_coef, ['y'], Nodefilter='Windload_10')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 11, years=1) * Q_coef, ['y'], Nodefilter='Windload_11')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 12, years=1) * Q_coef, ['y'], Nodefilter='Windload_12')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 13, years=1) * Q_coef, ['y'], Nodefilter='Windload_13')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 14, years=1) * Q_coef, ['y'], Nodefilter='Windload_14')
C01C_ELS_1YEAR.addForce(load_area * windload(1.5625 * 15, years=1) * Q_coef, ['y'], Nodefilter='Windload_15')

#Structure analysis
C01C_ELS_1YEAR.set_solverParameters(2000, 50, 5000)
C01C_ELS_1YEAR.calculate_structure()

#Plot and export data
# C01C_ELS_1YEAR.plotStructure()
C01C_ELS_1YEAR.export_Data('CASO_01C/CASO_01C_ELS_1YEAR.xlsx')
C01C_ELS_1YEAR.export_DxfStructure('CASO_01C/CASO_01C_ELS_1YEAR.dxf')

#%% ELS - Serviceability limit state WIND 2 YEARS

load_area = 4.8077
glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
G_coef = 1.0 # coeficiente de acciones permanentes ELS
Q_coef = 1.0 # coeficiente de acciones variables ELS

C01C_ELS_2YEAR = FTES.FTES_Structure()
C01C_ELS_2YEAR.import_DxfStructure('CASO_01C/CASO_01C.dxf')

#Set mechanical properties
C01C_ELS_2YEAR.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
C01C_ELS_2YEAR.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
C01C_ELS_2YEAR.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
C01C_ELS_2YEAR.set_Cable(Memberfilter='Allmembers')
C01C_ELS_2YEAR.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# Set restraints
C01C_ELS_2YEAR.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
#Set glazing load
C01C_ELS_2YEAR.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

#Set wind load
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 1, years=2) * Q_coef, ['y'], Nodefilter='Windload_01')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 2, years=2) * Q_coef, ['y'], Nodefilter='Windload_02')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 3, years=2) * Q_coef, ['y'], Nodefilter='Windload_03')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 4, years=2) * Q_coef, ['y'], Nodefilter='Windload_04')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 5, years=2) * Q_coef, ['y'], Nodefilter='Windload_05')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 6, years=2) * Q_coef, ['y'], Nodefilter='Windload_06')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 7, years=2) * Q_coef, ['y'], Nodefilter='Windload_07')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 8, years=2) * Q_coef, ['y'], Nodefilter='Windload_08')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 9, years=2) * Q_coef, ['y'], Nodefilter='Windload_09')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 10, years=2) * Q_coef, ['y'], Nodefilter='Windload_10')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 11, years=2) * Q_coef, ['y'], Nodefilter='Windload_11')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 12, years=2) * Q_coef, ['y'], Nodefilter='Windload_12')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 13, years=2) * Q_coef, ['y'], Nodefilter='Windload_13')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 14, years=2) * Q_coef, ['y'], Nodefilter='Windload_14')
C01C_ELS_2YEAR.addForce(load_area * windload(1.5625 * 15, years=2) * Q_coef, ['y'], Nodefilter='Windload_15')

#Structure analysis
C01C_ELS_2YEAR.set_solverParameters(2000, 50, 5000)
C01C_ELS_2YEAR.calculate_structure()

#Plot and export data
C01C_ELS_2YEAR.plotStructure(savePdf='CASO_01C/CASO_01C.pdf')
C01C_ELS_2YEAR.export_Data('CASO_01C/CASO_01C_ELS_2YEAR.xlsx')
C01C_ELS_2YEAR.export_DxfStructure('CASO_01C/CASO_01C_ELS_2YEAR.dxf')
