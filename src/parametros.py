from unidades import *

# Definición de propiedades de materiales y secciones
# fc = 210 * kgf / cm**2
# E = 15100 * (fc / (kgf / cm**2))**0.5 * kgf / cm**2
# G = 0.5 * E / (1 + 0.2)

# Propiedades del acero
fy = 4200 * kgf / cm**2  # Límite de fluencia del acero
Es = 2.1 * 10**6 * kgf / cm**2  # Módulo de elasticidad del acero

# Propiedades del concreto
fc = 210  # Resistencia a la compresión del concreto (kg/cm^2)
E = 15100 * fc**0.5 * kgf / cm**2  # Módulo de elasticidad del concreto
G = 0.5 * E / (1 + 0.2)  # Módulo de corte del concreto
fc = fc * kgf / cm**2  # Convertir fc a las unidades adecuadas
cover = 4 * cm  # Recubrimiento de vigas y columnas

# Parámetros del concreto no lineal
fc1 = -fc  # Resistencia a la compresión (negativa para simular compresión)
Ec1 = E  # Módulo de elasticidad del concreto
nuc1 = 0.2  # Coeficiente de Poisson del concreto
Gc1 = Ec1 / (2 * (1 + nuc1))  # Módulo de corte del concreto

# Propiedades del concreto confinado
Kfc = 1.0  # Ratio de resistencia del concreto confinado respecto al no confinado
Kres = 0.2  # Ratio de resistencia residual respecto a la máxima
fpc1 = Kfc * fc1  # Resistencia máxima del concreto confinado
epsc01 = 2 * fpc1 / Ec1  # Deformación unitaria correspondiente a fpc1
fpcu1 = Kres * fpc1  # Resistencia residual del concreto confinado
epsU1 = 5 * epsc01  # Deformación unitaria máxima del concreto confinado
lambda1 = 0.1  # Factor para definir el descenso de la curva de esfuerzo-deformación

# Propiedades del concreto no confinado
fpc2 = fc1  # Resistencia máxima del concreto no confinado
epsc02 = -0.003  # Deformación unitaria correspondiente a fpc2
fpcu2 = Kres * fpc2  # Resistencia residual del concreto no confinado
epsU2 = -0.006  # Deformación unitaria máxima del concreto no confinado

# Propiedades de resistencia a la tracción del concreto
ft1 = -0.14 * fpc1  # Resistencia a la tracción para concreto confinado
ft2 = -0.14 * fpc2  # Resistencia a la tracción para concreto no confinado
Ets = ft2 / 0.002  # Módulo de elasticidad en la zona de tracción

# Densidad del concreto
rho = 2400 * kg / m**3  # Densidad del concreto