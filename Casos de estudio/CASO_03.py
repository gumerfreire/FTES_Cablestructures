# -*- coding: utf-8 -*-
"""
FTES - Casos de analisis TFM
CASO 03
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

# load_area = 3.8462
# glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
# G_coef = 1.35 # coeficiente de acciones permanentes
# Q_coef = 1.5 # coeficiente de acciones variables



# C03_ELU = FTES.FTES_Structure()
# C03_ELU.import_DxfStructure('CASO_03/CASO_03.dxf')

# #Set mechanical properties
# C03_ELU.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
# C03_ELU.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
# C03_ELU.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
# C03_ELU.set_Cable(Memberfilter='Allmembers')
# C03_ELU.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# # Set restraints
# C03_ELU.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
# #Set glazing load
# C03_ELU.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

# #Set wind load
# C03_ELU.addForce(load_area * windload(1.9231 * 1) * Q_coef, ['y'], Nodefilter='Windload_01')
# C03_ELU.addForce(load_area * windload(1.9231 * 2) * Q_coef, ['y'], Nodefilter='Windload_02')
# C03_ELU.addForce(load_area * windload(1.9231 * 3) * Q_coef, ['y'], Nodefilter='Windload_03')
# C03_ELU.addForce(load_area * windload(1.9231 * 4) * Q_coef, ['y'], Nodefilter='Windload_04')
# C03_ELU.addForce(load_area * windload(1.9231 * 5) * Q_coef, ['y'], Nodefilter='Windload_05')
# C03_ELU.addForce(load_area * windload(1.9231 * 6) * Q_coef, ['y'], Nodefilter='Windload_06')
# C03_ELU.addForce(load_area * windload(1.9231 * 7) * Q_coef, ['y'], Nodefilter='Windload_07')
# C03_ELU.addForce(load_area * windload(1.9231 * 8) * Q_coef, ['y'], Nodefilter='Windload_08')
# C03_ELU.addForce(load_area * windload(1.9231 * 9) * Q_coef, ['y'], Nodefilter='Windload_09')
# C03_ELU.addForce(load_area * windload(1.9231 * 10) * Q_coef, ['y'], Nodefilter='Windload_10')
# C03_ELU.addForce(load_area * windload(1.9231 * 11) * Q_coef, ['y'], Nodefilter='Windload_11')
# C03_ELU.addForce(load_area * windload(1.9231 * 12) * Q_coef, ['y'], Nodefilter='Windload_12')

# #Structure analysis
# C03_ELU.set_solverParameters(2000, 50, 5000)
# C03_ELU.calculate_structure()

# #Plot and export data
# C03_ELU.plotStructure()
# C03_ELU.export_Data('CASO_03/CASO_03_ELU.xlsx')
# C03_ELU.export_DxfStructure('CASO_03/CASO_03_ELU.dxf')

# #%% ELS - Serviceability limit state

# load_area = 3.8462
# glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
# G_coef = 1.0 # coeficiente de acciones permanentes
# Q_coef = 1.0 # coeficiente de acciones variables



# C03_ELS = FTES.FTES_Structure()
# C03_ELS.import_DxfStructure('CASO_03/CASO_03.dxf')

# #Set mechanical properties
# C03_ELS.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
# C03_ELS.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
# C03_ELS.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
# C03_ELS.set_Cable(Memberfilter='Allmembers')
# C03_ELS.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# # Set restraints
# C03_ELS.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
# #Set glazing load
# C03_ELS.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

# #Set wind load
# C03_ELS.addForce(load_area * windload(1.9231 * 1) * Q_coef, ['y'], Nodefilter='Windload_01')
# C03_ELS.addForce(load_area * windload(1.9231 * 2) * Q_coef, ['y'], Nodefilter='Windload_02')
# C03_ELS.addForce(load_area * windload(1.9231 * 3) * Q_coef, ['y'], Nodefilter='Windload_03')
# C03_ELS.addForce(load_area * windload(1.9231 * 4) * Q_coef, ['y'], Nodefilter='Windload_04')
# C03_ELS.addForce(load_area * windload(1.9231 * 5) * Q_coef, ['y'], Nodefilter='Windload_05')
# C03_ELS.addForce(load_area * windload(1.9231 * 6) * Q_coef, ['y'], Nodefilter='Windload_06')
# C03_ELS.addForce(load_area * windload(1.9231 * 7) * Q_coef, ['y'], Nodefilter='Windload_07')
# C03_ELS.addForce(load_area * windload(1.9231 * 8) * Q_coef, ['y'], Nodefilter='Windload_08')
# C03_ELS.addForce(load_area * windload(1.9231 * 9) * Q_coef, ['y'], Nodefilter='Windload_09')
# C03_ELS.addForce(load_area * windload(1.9231 * 10) * Q_coef, ['y'], Nodefilter='Windload_10')
# C03_ELS.addForce(load_area * windload(1.9231 * 11) * Q_coef, ['y'], Nodefilter='Windload_11')
# C03_ELS.addForce(load_area * windload(1.9231 * 12) * Q_coef, ['y'], Nodefilter='Windload_12')

# #Structure analysis
# C03_ELS.set_solverParameters(2000, 50, 5000)
# C03_ELS.calculate_structure()

# #Plot and export data
# C03_ELS.plotStructure()
# C03_ELS.export_Data('CASO_03/CASO_03_ELS.xlsx')
# C03_ELS.export_DxfStructure('CASO_03/CASO_03_ELS.dxf')

#%% ELS - Serviceability limit state WIND 1 YEAR


load_area = 3.8462
glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
G_coef = 1.0 # coeficiente de acciones permanentes
Q_coef = 1.0 # coeficiente de acciones variables



C03_ELS_1YEAR = FTES.FTES_Structure()
C03_ELS_1YEAR.import_DxfStructure('CASO_03/CASO_03.dxf')

#Set mechanical properties
C03_ELS_1YEAR.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
C03_ELS_1YEAR.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
C03_ELS_1YEAR.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
C03_ELS_1YEAR.set_Cable(Memberfilter='Allmembers')
C03_ELS_1YEAR.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# Set restraints
C03_ELS_1YEAR.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
#Set glazing load
C03_ELS_1YEAR.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

#Set wind load
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 1, years=1) * Q_coef, ['y'], Nodefilter='Windload_01')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 2, years=1) * Q_coef, ['y'], Nodefilter='Windload_02')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 3, years=1) * Q_coef, ['y'], Nodefilter='Windload_03')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 4, years=1) * Q_coef, ['y'], Nodefilter='Windload_04')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 5, years=1) * Q_coef, ['y'], Nodefilter='Windload_05')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 6, years=1) * Q_coef, ['y'], Nodefilter='Windload_06')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 7, years=1) * Q_coef, ['y'], Nodefilter='Windload_07')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 8, years=1) * Q_coef, ['y'], Nodefilter='Windload_08')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 9, years=1) * Q_coef, ['y'], Nodefilter='Windload_09')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 10, years=1) * Q_coef, ['y'], Nodefilter='Windload_10')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 11, years=1) * Q_coef, ['y'], Nodefilter='Windload_11')
C03_ELS_1YEAR.addForce(load_area * windload(1.9231 * 12, years=1) * Q_coef, ['y'], Nodefilter='Windload_12')

#Structure analysis
C03_ELS_1YEAR.set_solverParameters(2000, 50, 5000)
C03_ELS_1YEAR.calculate_structure()

#Plot and export data
# C03_ELS_1YEAR.plotStructure()
C03_ELS_1YEAR.export_Data('CASO_03/CASO_03_ELS_1YEAR.xlsx')
C03_ELS_1YEAR.export_DxfStructure('CASO_03/CASO_03_ELS_1YEAR.dxf')

#%% ELS - Serviceability limit state WIND 2 YEARS

load_area = 3.8462
glassload = 55 * 9.81 #glass load N/m2 (for 55 kg/m2)
G_coef = 1.0 # coeficiente de acciones permanentes
Q_coef = 1.0 # coeficiente de acciones variables



C03_ELS_2YEAR = FTES.FTES_Structure()
C03_ELS_2YEAR.import_DxfStructure('CASO_03/CASO_03.dxf')

#Set mechanical properties
C03_ELS_2YEAR.set_E(160e9, Memberfilter='Allmembers') #galvanized steel cable
C03_ELS_2YEAR.set_A(634/1e6, Memberfilter='Allmembers') #set area for cable 31mm
C03_ELS_2YEAR.set_sw(5.3*9.81, Memberfilter='Allmembers') #set self weight for cable 31 mm
C03_ELS_2YEAR.set_Cable(Memberfilter='Allmembers')
C03_ELS_2YEAR.set_P(75000, Memberfilter='Allmembers') #set the same pretension for all cables

# Set restraints
C03_ELS_2YEAR.addRestraint(['x','y','z'], Nodefilter='RestraintsXYZ')
#Set glazing load
C03_ELS_2YEAR.addForce(-load_area * glassload * G_coef, ['z'], Nodefilter='Glassload')

#Set wind load
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 1, years=2) * Q_coef, ['y'], Nodefilter='Windload_01')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 2, years=2) * Q_coef, ['y'], Nodefilter='Windload_02')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 3, years=2) * Q_coef, ['y'], Nodefilter='Windload_03')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 4, years=2) * Q_coef, ['y'], Nodefilter='Windload_04')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 5, years=2) * Q_coef, ['y'], Nodefilter='Windload_05')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 6, years=2) * Q_coef, ['y'], Nodefilter='Windload_06')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 7, years=2) * Q_coef, ['y'], Nodefilter='Windload_07')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 8, years=2) * Q_coef, ['y'], Nodefilter='Windload_08')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 9, years=2) * Q_coef, ['y'], Nodefilter='Windload_09')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 10, years=2) * Q_coef, ['y'], Nodefilter='Windload_10')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 11, years=2) * Q_coef, ['y'], Nodefilter='Windload_11')
C03_ELS_2YEAR.addForce(load_area * windload(1.9231 * 12, years=2) * Q_coef, ['y'], Nodefilter='Windload_12')

#Structure analysis
C03_ELS_2YEAR.set_solverParameters(2000, 50, 5000)
C03_ELS_2YEAR.calculate_structure()

#Plot and export data
C03_ELS_2YEAR.plotStructure(savePdf='CASO_03/CASO_03.pdf')
C03_ELS_2YEAR.export_Data('CASO_03/CASO_03_ELS_2YEAR.xlsx')
C03_ELS_2YEAR.export_DxfStructure('CASO_03/CASO_03_ELS_2YEAR.dxf')