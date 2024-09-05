from modelado.ModeloRegular import ModeloRegular
from modelado.SeccionesGeometricas import SeccionRectangular
from analisis.AnalisisEstatico import AnalisisEstatico
from analisis.ModalAnalysis import ModalAnalysis
from unidades import m, kg, cm, kgf, concrete_density, Pa
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import opsvis as opsv
import numpy as np
from math import sqrt

# Definir resistencia del concreto (fc) y propiedades elásticas
fc = 210 * kg / cm**2
E = 151 * fc**0.5 * kgf / cm**2  # Módulo de elasticidad
G = 0.5 * E / (1 + 0.2)  # Módulo de rigidez

# Secciones de elementos
viga = SeccionRectangular(0.30 * m, 0.60 * m)  # Sección de viga (30 cm x 60 cm)
columna = SeccionRectangular(0.60 * m, 0.60 * m)  # Sección de columna (60 cm x 60 cm)

# Propiedades geométricas de la viga
Av = viga.calcular_area()  # Área (m²)
Iz_v, Iy_v = viga.calcular_momentos_inercia()  # Momentos de inercia (m⁴)
Jxx_v = viga.calcular_momento_polar()  # Momento polar de inercia (m⁴)

# Propiedades geométricas de la columna
Ac = columna.calcular_area()  # Área (m²)
Iz_c, Iy_c = columna.calcular_momentos_inercia()  # Momentos de inercia (m⁴)
Jxx_c = columna.calcular_momento_polar()  # Momento polar de inercia (m⁴)

# Parámetros del modelo estructural
RigidDiaphragm = 'ON'
dx, dy, dz = 4.0 * m, 4.0 * m, 3.0 * m  # Dimensiones en metros
nx, ny, nz = 5, 4, 8 # 2, 3, 2 #5, 4, 8  # Número de divisiones en X, Y, Z

# Crear el objeto del modelo estructural
portico = ModeloRegular(dx, dy, dz, nx, ny, nz)

# Generar nodos, elementos y diafragmas
Nodes = portico.generar_nodos()
Elems = portico.generar_elementos_vigas_columnas()
Diap = portico.generar_diaphragmas()

# Inicializar el modelo en OpenSees
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

# Crear los nodos en OpenSees
for Ni in Nodes:
    ops.node(int(Ni[0]), *Ni[1:4])

# Definir diafragmas rígidos
if RigidDiaphragm == 'ON':
    dirDia = 3  # Eje perpendicular al plano del diafragma
    for Nd in Diap:
        ops.node(int(Nd[0]), *Nd[1:4])
        ops.fix(int(Nd[0]), *[0, 0, 1, 1, 1, 0])
        NodesDi = [int(Ni[0]) for Ni in Nodes if Ni[3] == Nd[3]]
        ops.rigidDiaphragm(dirDia, int(Nd[0]), *NodesDi)

# Definir restricciones en Z=0
ops.fixZ(0.0, *[1, 1, 1, 1, 1, 1], '-tol', 1e-6)

# Definir transformaciones geométricas
ops.geomTransf('PDelta', 1, *[1, 0, 0])
ops.geomTransf('Linear', 2, *[1, -1, 0])

# Crear elementos en OpenSees
for Ele in Elems:
    if int(Ele[3]) == 1:  # Columna
        ops.element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Ac, E, G, Jxx_c, Iy_c, Iz_c, int(Ele[3]), '-mass', concrete_density * Ac)
    else:  # Viga
        ops.element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Av, E, G, Jxx_v, Iy_v, Iz_v, int(Ele[3]), '-mass', concrete_density * Av)

# # Visualizar el modelo
# opsv.plot_model(fig_wi_he=(30., 30.), az_el=(-20, 20), node_labels=0, element_labels=1, local_axes=False)
# plt.show()

# Definir masa por unidad de área
wLive = 250 * kg / m**2
wLosa = 300 * kg / m**2
wAcab = 100 * kg / m**2
wTabi = 150 * kg / m**2
wTotal = 1.0 * (wLosa + wAcab + wTabi) + 0.25 * wLive

# Aplicar masa
Carga = wTotal * dx * dy * m**2
for Ni in Nodes:
    ops.mass(int(Ni[0]), Ni[4] * Carga, Ni[4] * Carga, 0.0)
Nmodes = 3 * nz

# define some analysis settings
ops.wipeAnalysis()
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormUnbalance", 0.0001, 10)
ops.algorithm("Linear")
ops.integrator("LoadControl", 0.0)
ops.analysis("Static")

# run the eigenvalue analysis with 7 modes
# and obtain the eigenvalues
eigs = ops.eigen("-genBandArpack", Nmodes)

# compute the modal properties
ops.modalProperties("-file", "ModalReport.txt", "-unorm")

vals = ops.eigen(Nmodes)
Tmodes = np.zeros(len(vals))
for i in range(Nmodes):
    Tmodes[i] = 2*np.pi/vals[i]**0.5
    print("T[%i]: %.5f"%(i+1,Tmodes[i]))

# Realizamos un análisis para obtener la matriz de Masas
ops.wipeAnalysis()
ops.constraints('Transformation') 
ops.numberer("Plain")
ops.system("FullGeneral")
ops.algorithm('Linear')
ops.integrator('GimmeMCK',1.0,0.0,0.0)
ops.analysis('Transient')
ops.analyze(1,0.0) 

# Obtenemos la matriz de Masas
N = ops.systemSize()         # Número de Grados de Libertad
Mmatrix = ops.printA('-ret')
Mmatrix = np.array(Mmatrix).reshape((N,N))
MF = Mmatrix[-3*nz:,-3*nz:]

np.set_printoptions(precision=3,linewidth=300,suppress=True)
H = np.arange(1,nz+1)*dz
P = sum(MF[0::3,0::3])*9.80665 # Peso por nivel
Ro = 8.

# Instanciar la clase y realizar el análisis
SisEst = AnalisisEstatico(Tmodes, Z=0.45, U=1.5, S=1.0, Tp=0.4, Tl=2.5, R=Ro)

# Calcular el espectro y las cargas
zucs_r = SisEst.calc_espectro()
F_lat, k = SisEst.calc_cargas(zucs_r[0], P, H)

# Ejecutar el análisis estático en X
SisEstXX = SisEst.aplicar_y_analizar(F_lat, Diap, nx, ny, nz, dx, dy, dz, Ro, dir='x')

# Ejecutar el análisis estático en Y
SisEstYY = SisEst.aplicar_y_analizar(F_lat, Diap, nx, ny, nz, dx, dy, dz, Ro, dir='y')
