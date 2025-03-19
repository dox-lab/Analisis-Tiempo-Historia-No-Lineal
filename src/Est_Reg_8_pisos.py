from modelado.ModeloRegular import ModeloRegular
from modelado.SeccionesGeometricas import SeccionRectangular
from analisis.AnalisisEstatico import AnalisisEstatico
from analisis.CalculoMasasEfectivas import CalculoMasasEfectivas
from analisis.AnalisisModalEspectral import AnalisisModalEspectral
from unidades import m, kg, cm, kgf, concrete_density
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import opsvis as opsv
import numpy as np

# Propiedades del material
fc = 210 * kg / cm**2  # Resistencia del concreto
E = 151 * fc**0.5 * kgf / cm**2  # Módulo de elasticidad del concreto
G = 0.5 * E / (1 + 0.2)  # Módulo de rigidez del concreto

# Definición de secciones geométricas
viga = SeccionRectangular(0.30 * m, 0.60 * m)  # Sección de viga
columna = SeccionRectangular(0.60 * m, 0.60 * m)  # Sección de columna

# Cálculo de propiedades geométricas
Av, Iz_v, Iy_v, Jxx_v = viga.calcular_area(), *viga.calcular_momentos_inercia(), viga.calcular_momento_polar()
Ac, Iz_c, Iy_c, Jxx_c = columna.calcular_area(), *columna.calcular_momentos_inercia(), columna.calcular_momento_polar()

# Parámetros del modelo estructural
dx, dy, dz = 4.0 * m, 4.0 * m, 3.0 * m  # Dimensiones en metros
nx, ny, nz = 5, 5, 8  # Número de divisiones en X, Y, Z
RigidDiaphragm = 'ON'  # Activación de diafragmas rígidos
Ro = 8.0  # Factor de sobrerresistencia

# Creación del modelo estructural
portico = ModeloRegular(dx, dy, dz, nx, ny, nz)

# Generación de nodos, elementos y diafragmas
Nodes = portico.generar_nodos()
Elems = portico.generar_elementos_vigas_columnas()
Diap = portico.generar_diaphragmas()

# Inicialización de OpenSees
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

# Creación de nodos en OpenSees
for Ni in Nodes:
    ops.node(int(Ni[0]), *Ni[1:4])

# Definición de diafragmas rígidos
if RigidDiaphragm == 'ON':
    dirDia = 3  # Eje perpendicular al plano del diafragma
    for Nd in Diap:
        ops.node(int(Nd[0]), *Nd[1:4])
        ops.fix(int(Nd[0]), *[0, 0, 1, 1, 1, 0])
        NodesDi = [int(Ni[0]) for Ni in Nodes if Ni[3] == Nd[3]]
        ops.rigidDiaphragm(dirDia, int(Nd[0]), *NodesDi)

# Restricciones en la base (Z=0)
ops.fixZ(0.0, *[1, 1, 1, 1, 1, 1], '-tol', 1e-6)

# Definición de transformaciones geométricas
ops.geomTransf('PDelta', 1, *[1, 0, 0])
ops.geomTransf('Linear', 2, *[1, -1, 0])

# Creación de elementos en OpenSees (vigas y columnas)
for Ele in Elems:
    props = (Ac, E, G, Jxx_c, Iy_c, Iz_c) if int(Ele[3]) == 1 else (Av, E, G, Jxx_v, Iy_v, Iz_v)
    ops.element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), *props, int(Ele[3]), '-mass', concrete_density * props[0])

# Visualización del modelo
opsv.plot_model(fig_wi_he=(30., 30.), az_el=(-20, 20), node_labels=0, element_labels=0, local_axes=False)
plt.show()

# Definición de cargas
wLive = 250 * kg / m**2  # Carga viva
wLosa = 300 * kg / m**2  # Carga de losa
wAcab = 100 * kg / m**2  # Carga de acabados
wTabi = 150 * kg / m**2  # Carga de tabiques
wTotal = 1.0 * (wLosa + wAcab + wTabi) + 0.25 * wLive  # Carga total
Carga = wTotal * dx * dy * m**2  # Carga aplicada

# Aplicación de masa a los nodos
for Ni in Nodes:
    ops.mass(int(Ni[0]), Ni[4] * Carga, Ni[4] * Carga, 0.0)

Nmodes = 3 * nz  # Número de modos

# Configuración del análisis estático
ops.wipeAnalysis()
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormUnbalance", 0.0001, 10)
ops.algorithm("Linear")
ops.integrator("LoadControl", 0.0)
ops.analysis("Static")

# Análisis de valores propios (modos de vibración)
eigs = ops.eigen("-genBandArpack", Nmodes)
ops.modalProperties("-file", "resultados/Est_Reg_8_pisos.txt", "-unorm")

# Cálculo de modos y periodos
vals = ops.eigen(Nmodes)
Tmodes = np.array([2 * np.pi / val**0.5 for val in vals])
for i, T in enumerate(Tmodes):
    print(f"T[{i+1}]: {T:.5f}")

# Obtención de la matriz de masas del sistema
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer("Plain")
ops.system("FullGeneral")
ops.algorithm('Linear')
ops.integrator('GimmeMCK', 1.0, 0.0, 0.0)
ops.analysis('Transient')
ops.analyze(1, 0.0)

# Procesamiento de la matriz de masas
Mmatrix = np.array(ops.printA('-ret')).reshape((ops.systemSize(), ops.systemSize()))
Mf = Mmatrix[-3*nz:, -3*nz:]  # Matriz de masas filtrada por grados de libertad

# Configuración del análisis estático
H = np.arange(1, nz + 1) * dz  # Altura de cada nivel
P = np.sum(Mf[0::3, 0::3]) * 9.80665  # Peso por nivel

# Realización del análisis estático
SisEst = AnalisisEstatico(Tmodes, Z=0.45, U=1.0, S=1.0, Tp=0.4, Tl=2.5, R=Ro)
zucs_r = SisEst.calc_espectro()  # Cálculo del espectro
F_lat, k = SisEst.calc_cargas(zucs_r[0], P, H)  # Cálculo de cargas sísmicas

# Ejecución del análisis estático en X e Y
corte_acum_x, _ = SisEst.aplicar_y_analizar(F_lat, Diap, nx, ny, nz, dx, dy, dz, Ro, dir='x')
corte_acum_y, _ = SisEst.aplicar_y_analizar(F_lat, Diap, nx, ny, nz, dx, dy, dz, Ro, dir='y')

# Obtención de modos de vibración
Tags = ops.getNodeTags()
modo = np.zeros((Nmodes, 3 * nz))
for j in range(1, Nmodes + 1):
    ind = 0
    for i in Tags[-nz:]:
        temp = ops.nodeEigenvector(i, j)
        modo[j - 1, [ind, ind + 1, ind + 2]] = temp[0], temp[1], temp[-1]
        ind += 3

# Cálculo de masas efectivas
calculo_masas = CalculoMasasEfectivas(Nmodes, nz, Mf, Tmodes)
masas_efectivas, num_modos_min = calculo_masas.calcular_masas(Tags)
num_modos_min = 24

# Análisis modal espectral
analisis_modal = AnalisisModalEspectral(zucs_r, Mf, modo, Tmodes, nz, num_modos_min, Ro, corte_acum_x, corte_acum_y, dz)

# Combinación modal
disp_comb_x, drift_comb_x, shear_comb_x, disp_comb_y, drift_comb_y, shear_comb_y, a_din_sin_escalar = analisis_modal.combinacion_modal()

# Escalamiento de resultados
scale_factor_x, scale_factor_y, message_x, message_y = analisis_modal.escalamiento(shear_comb_x, shear_comb_y)
print(message_x)
print(message_y)

# Aplicación de factores de escala y obtención de resultados finales
ADME = analisis_modal.aplicar_escala(scale_factor_x, scale_factor_y, disp_comb_x, drift_comb_x, shear_comb_x, disp_comb_y, drift_comb_y, shear_comb_y)

# Graficación de distorsiones
analisis_modal.plot_derivas(ADME)

print(Mf)