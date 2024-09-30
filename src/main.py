from modelado.DXFProcessor import DXFProcessor
from modelado.SectionCreator import StructureModel
from modelado.SeccionesGeometricas import SeccionRectangular, Muro
from analisis.ReadSMDFile import ReadSMDFile
from analisis.ReadRecord import ReadRecord
import openseespy.opensees as ops
import math
import matplotlib.pyplot as plt
import numpy as np
import vfo.vfo as vfo
from parametros import *


# Procesar datos de la estructura desde el archivo DXF
data_estructura = DXFProcessor("src/modelado/Estructura_Irregular_8_Pisos.dxf")
nodes = data_estructura['Nodes']
elements = data_estructura['Elems']
sec_vigas = data_estructura['sec_vig']
sec_columnas = data_estructura['sec_col']
num_vigas = data_estructura['nsec_vig']
num_columnas = data_estructura['nsec_col']
diafragmas = data_estructura['Diap']
aPlanta = data_estructura['AreaPlanta']
rest_node = data_estructura['Rest_node']
levels = data_estructura['Levels']
nz = len(data_estructura["Levels"]) - 1
dz = levels[1] - levels[0]

# Inicializar el modelo en OpenSees
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)
rigid_diaphragm = True  # Usar diafragmas rígidos

# Definición de materiales
ops.uniaxialMaterial('Elastic', 3, G)  # Concreto Cortante
ops.uniaxialMaterial('Concrete02', 4, fpc1, epsc01, fpcu1, epsU1, lambda1, ft1, Ets)  # Concreto confinado
ops.uniaxialMaterial('Concrete02', 5, fpc2, epsc02, fpcu2, epsU2, lambda1, ft2, Ets)  # Concreto no confinado
ops.uniaxialMaterial('Steel02', 6, fy, Es, 0.01, 18, 0.925, 0.15)  # Acero de refuerzo

# Crear las secciones para columnas y vigas
columnas = [SeccionRectangular(col[0], col[1]) for col in sec_columnas]
vigas = [SeccionRectangular(vig[0], vig[1]) for vig in sec_vigas]

# Crear el modelo de la estructura
structure_model = StructureModel(num_columnas, num_vigas, columnas, vigas, G, cover)
structure_model.create_structure()

# Crear nodos en OpenSees
for nodo in nodes:
    ops.node(int(nodo[0]), *nodo[1:4])

# Definir diafragmas rígidos si está activado
if rigid_diaphragm:
    dir_diaphragm = 3  # Eje perpendicular al plano del diafragma
    for d in diafragmas:
        ops.node(int(d[0]), *d[1:4])
        ops.fix(int(d[0]), *[0, 0, 1, 1, 1, 0])
        nodes_at_level = [int(nodo[0]) for nodo in nodes if nodo[3] == d[3]]
        ops.rigidDiaphragm(dir_diaphragm, int(d[0]), *nodes_at_level)

# Fijar los nodos en la base (Z=0)
ops.fixZ(0.0, *[1, 1, 1, 1, 1, 1], '-tol', 1e-6)

# Definir transformaciones geométricas
ops.geomTransf('PDelta', 1, *[1, 0, 0])
ops.geomTransf('Linear', 2, *[1, -1, 0])

# Integración de elementos usando la regla de Lobatto
num_int_points = 7
for i in range(num_columnas):
    ops.beamIntegration('Lobatto', i + 1, i + 1, num_int_points)
for i in range(num_vigas):
    ops.beamIntegration('Lobatto', i + 1 + num_columnas, i + 1 + num_columnas, num_int_points)

# Crear elementos en OpenSees
for elem in elements:
    elem_type = 'forceBeamColumn'
    if int(elem[5]) == 1:  # Columna
        area_columna = elem[3] * elem[4]
        ops.element(elem_type, int(elem[0]), int(elem[1]), int(elem[2]), 1, int(elem[6]), '-mass', rho * area_columna * m**2)
    elif int(elem[5]) == 2:  # Viga
        area_viga = elem[3] * elem[4]
        ops.element(elem_type, int(elem[0]), int(elem[1]), int(elem[2]), 2, int(elem[6]), '-mass', rho * area_viga * m**2)
    elif int(elem[5]) == 3:  # Muro en x
        muro_x = Muro(abs(nodes[int(elem[1])][1] - nodes[int(elem[2])][1]), elem[6])
        nf_x, thk_x, w_list_x, rho_x, conc_x, steel_x = muro_x.prop_muros()
        ops.element('MVLEM_3D', int(elem[0]), int(elem[1]), int(elem[2]), int(elem[3]), int(elem[4]), nf_x,
                    '-thick', *thk_x, '-width', *w_list_x, '-rho', *rho_x, '-matConcrete', *conc_x, '-matSteel', *steel_x, '-matShear', 3,
                    '-Poisson', 0.2, '-Density', rho)
    elif int(elem[5]) == 4:  # Muro en y
        muro_y = Muro(abs(nodes[int(elem[1])][2] - nodes[int(elem[2])][2]), elem[6])
        nf_y, thk_y, w_list_y, rho_y, conc_y, steel_y = muro_y.prop_muros()
        ops.element('MVLEM_3D', int(elem[0]), int(elem[1]), int(elem[2]), int(elem[3]), int(elem[4]), nf_y,
                    '-thick', *thk_y, '-width', *w_list_y, '-rho', *rho_y, '-matConcrete', *conc_y, '-matSteel', *steel_y, '-matShear', 3,
                    '-Poisson', 0.2, '-Density', rho)

# # Visualizar el modelo
# vfo.plot_model(show_nodes="yes", line_width=5)
# plt.show()

# ASIGNACIÓN DE MASAS Y MODOS DE VIBRACIÓN
wLive = 250*kg/m**2
wLosa = 300*kg/m**2
wAcab = 100*kg/m**2
wTabi = 150*kg/m**2
wTotal = 1.0*(wLosa+wAcab+wTabi)+0.25*wLive

Carga = wTotal*aPlanta*m**2
for ni in nodes:
    ops.mass(int(ni[0]), ni[4] * Carga/((len(nodes)-len(rest_node))/nz), ni[4] *Carga/((len(nodes)-len(rest_node))/nz), 0.0)

Nmodes = 3*nz

vals = ops.eigen(Nmodes)
Tmodes = np.zeros(len(vals))
for i in range(Nmodes):
    Tmodes[i] = 2*np.pi/vals[i]**0.5
    print("T[%i]: %.5f" % (i+1, Tmodes[i]))

# vfo.plot_modeshape(modenumber=1, scale=2500, line_width=3)
# vfo.plot_modeshape(modenumber=2, scale=2500, line_width=3)
# vfo.plot_modeshape(modenumber=3, scale=2500, line_width=3)
# plt.show()

vfo.createODB('src/modelado/3D_Building', "Gravity", Nmodes=3)

ops.timeSeries('Linear', 1)  # tag
ops.pattern('Plain', 1, 1)  # tag, timeSeries_tag

# Peso propio de los elementos
for ele in elements:
    if int(ele[5]) == 1:  # 1 Columna
        Ac = ele[3]*ele[4]
        ops.eleLoad('-ele', int(ele[0]), '-type',
                '-beamUniform', 0.0, 0.0, -rho*Ac*g)  # Wy, Wz, Wx
    elif int(ele[5]) == 2:  # 2 viga
        Av = ele[3]*ele[4]
        ops.eleLoad('-ele', int(ele[0]), '-type',
                '-beamUniform', -rho*Av*g, 0.0, 0.0)
    elif int(ele[5]) == 3:  # 3 Muro en x
        Lmx = abs(nodes[int(ele[1])][1]-nodes[int(ele[2])][1])
        t = ele[6]
        ops.load(int(ele[1]), 0., 0., -0.25*dz*Lmx*t*rho*g, 0., 0., 0.)
        ops.load(int(ele[2]), 0., 0., -0.25*dz*Lmx*t*rho*g, 0., 0., 0.)
        ops.load(int(ele[3]), 0., 0., -0.25*dz*Lmx*t*rho*g, 0., 0., 0.)
        ops.load(int(ele[4]), 0., 0., -0.25*dz*Lmx*t*rho*g, 0., 0., 0.)
    elif int(ele[5]) == 4:  # 4 Muro en y
        Lmy = abs(nodes[int(ele[1])][2]-nodes[int(ele[2])][2])
        t = ele[6]
        ops.load(int(ele[1]), 0., 0., -0.25*dz*Lmy*t*rho*g, 0., 0., 0.)
        ops.load(int(ele[2]), 0., 0., -0.25*dz*Lmy*t*rho*g, 0., 0., 0.)
        ops.load(int(ele[3]), 0., 0., -0.25*dz*Lmy*t*rho*g, 0., 0., 0.)
        ops.load(int(ele[4]), 0., 0., -0.25*dz*Lmy*t*rho*g, 0., 0., 0.)
    else:
        print("Error!. Es un elemento no definido!")

# Carga viva y muerta sobre los nodos

for ni in nodes:
    if int(ni[3]) != 0:
        ops.load(int(ni[0]), 0.0, 0.0, -ni[4]*Carga*g/((len(nodes)-len(rest_node))/nz), 0.0, 0.0, 0.0)

# GRAVITY-ANALYSIS PARAMETERS -------------------------------------------------------------
Tol = 1.0e-8  # convergence tolerance for test
constraintsTypeGravity = "Plain"
if rigid_diaphragm:
    constraintsTypeGravity = "Transformation"
ops.constraints(constraintsTypeGravity)  # how it handles boundary conditions
ops.numberer('RCM')  # renumber dof's to minimize band-width (optimization)
ops.system('BandGeneral')  # system of equations
ops.test('EnergyIncr', Tol, 6)  # test convergence
ops.algorithm('Newton')  # Newton's solution algorithm
NstepGravity = 10  # apply gravity in 10 steps
DGravity = 1.0 / NstepGravity  # first load increment
ops.integrator('LoadControl', DGravity)
ops.analysis('Static')  # define type of analysis (static or transient)
ops.analyze(NstepGravity)  # apply gravity

# Maintain constant gravity loads and reset time to zero
ops.loadConst('-time', 0.0)

ops.wipeAnalysis()
vfo.createODB('src/modelado/3D_Building', "Dynamic_GM1")

# --------------------------------------------------------------------------------------------------
# Example 8. 3D Uniform Earthquake Excitation
# Execute this after the model is built and gravity is applied

# Set up ground-motion parameters
GMdirection = 1  # ground-motion direction
GMfile = "elCentro"  # ground-motion filenames
GMdir = "data/sismos"  # ground-motion directory
GMfact = 1  # ground-motion scaling factor

# Set up ground-motion-analysis parameters
DtAnalysis = 0.02  # time-step for lateral analysis
TmaxAnalysis = 30.0  # maximum duration of ground-motion analysis

# Configurar parámetros de análisis
# Configurar parámetros de análisis ---------------------------------------------
# MANEJADOR DE RESTRICCIONES -- Determina cómo se aplican las ecuaciones de restricción en el análisis
constraintsTypeDynamic = "Transformation"  # Definir el tipo de restricciones
ops.constraints(constraintsTypeDynamic)  # Aplicar las restricciones

# NUMERADOR DE GRADOS DE LIBERTAD: determina la asignación entre los números de ecuaciones y los grados de libertad
numbererTypeDynamic = "RCM"  # Tipo de numerador de DOF
ops.numberer(numbererTypeDynamic)  # Aplicar numeración de grados de libertad

# SISTEMA -- Solucionadores de ecuaciones lineales (cómo almacenar y resolver el sistema de ecuaciones)
systemTypeDynamic = "BandGeneral"  # Probar UmfPack para problemas grandes
ops.system(systemTypeDynamic)  # Aplicar el sistema

# TEST: prueba de convergencia 
# Convergence TEST -- determina si se ha logrado la convergencia al final de un paso de iteración
TolDynamic = 1e-8  # Tolerancia para la prueba de convergencia
maxNumIterDynamic = 10  # Número máximo de iteraciones antes de fallar
printFlagDynamic = 0  # Indicador para imprimir información de convergencia (opcional)
testTypeDynamic = "EnergyIncr"  # Tipo de prueba de convergencia
ops.test(testTypeDynamic, TolDynamic, maxNumIterDynamic, printFlagDynamic)  # Aplicar la prueba de convergencia

# Para procedimientos de mejora de convergencia:
maxNumIterConvergeDynamic = 2000
printFlagConvergeDynamic = 0

# ALGORITMO DE SOLUCIÓN: -- Iterar desde el último paso de tiempo al actual
algorithmTypeDynamic = "ModifiedNewton"  # Tipo de algoritmo de solución
ops.algorithm(algorithmTypeDynamic)  # Aplicar el algoritmo

# INTEGRADOR ESTÁTICO: -- Determina el siguiente paso de tiempo para un análisis
# INTEGRADOR TRANSITORIO: -- Determina el siguiente paso de tiempo para un análisis incluyendo efectos inerciales
NewmarkGamma = 0.5  # Parámetro gamma del integrador de Newmark (también HHT)
NewmarkBeta = 0.25  # Parámetro beta del integrador de Newmark
integratorTypeDynamic = "Newmark"  # Tipo de integrador
ops.integrator(integratorTypeDynamic, NewmarkGamma, NewmarkBeta)  # Aplicar el integrador

# ANÁLISIS -- Define qué tipo de análisis se va a realizar
analysisTypeDynamic = "Transient"  # Tipo de análisis
ops.analysis(analysisTypeDynamic)  # Ejecutar el análisis

# ----------- Set up analysis parameters
# Here, you'd need to implement or source the parameters defined in LibAnalysisDynamicParameters
# Assuming you've set up similar commands with constraints, DOFnumberer, system, convergenceTest, etc.

# ------------ Define & apply damping (RAYLEIGH damping)
xDamp = 0.02  # damping ratio
MpropSwitch = 1.0
KcurrSwitch = 0.0
KcommSwitch = 1.0
KinitSwitch = 0.0
nEigenI = 1  # mode 1
nEigenJ = 3  # mode 3

# Eigenvalue analysis for nEigenJ modes
lambdaN = ops.eigen(nEigenJ)
lambdaI = lambdaN[nEigenI - 1]  # eigenvalue mode i
lambdaJ = lambdaN[nEigenJ - 1]  # eigenvalue mode j
omegaI = math.sqrt(lambdaI)
omegaJ = math.sqrt(lambdaJ)

# Calculate Rayleigh damping coefficients
alphaM = MpropSwitch * xDamp * (2 * omegaI * omegaJ) / (omegaI + omegaJ)
betaKcurr = KcurrSwitch * 2.0 * xDamp / (omegaI + omegaJ)
betaKcomm = KcommSwitch * 2.0 * xDamp / (omegaI + omegaJ)
betaKinit = KinitSwitch * 2.0 * xDamp / (omegaI + omegaJ)

# Apply Rayleigh damping
ops.rayleigh(alphaM, betaKcurr, betaKinit, betaKcomm)

#  --------------------------------- Perform Dynamic Ground-Motion Analysis
IDloadTag = 2  # for uniformSupport excitation
inFile = f"{GMdir}/{GMfile}.at2"
outFile = f"{GMdir}/{GMfile}.g3"  # Output file

# Read the ground-motion file (This would be a custom function to handle file reading and formatting)
dt, npts= ReadRecord(inFile, outFile)  # Ensure this function is implemented in Python
# Read the file line by line and extract all numerical values into a list

# signal_data = []

# with open(outFile, 'r') as file:
#     for line in file:
#         # Split the line into individual numbers and convert them to float
#         numbers = line.split()
#         signal_data.extend([float(num) for num in numbers])

# # Convert the list to a NumPy array for further processing
# signal = np.array(signal_data)

# # Create a time vector assuming a dt of 0.005 seconds
# time = np.arange(0, len(signal) * dt, dt)

# # Plot the seismic signal
# plt.figure(figsize=(10, 6))
# plt.plot(time, signal, label="Señal sísmica")
# plt.title("Señal Sísmica - Aceleración vs Tiempo")
# plt.xlabel("Tiempo (s)")
# plt.ylabel("Aceleración (g)")
# plt.grid(True)
# plt.legend()
# plt.show()

GMfatt = g * GMfact  # Datos del archivo de entrada en unidades de g -- aceleración
ops.timeSeries('Path', IDloadTag, '-dt', dt, '-filePath', outFile, '-factor', GMfatt) # define acceleration vector from file (dt=0.005 is associated with the input file gm)
ops.pattern("UniformExcitation", IDloadTag, GMdirection, "-accel", IDloadTag)  # Crear excitación uniforme

Nsteps = int(TmaxAnalysis / DtAnalysis)
ok = ops.analyze(Nsteps, DtAnalysis)  # Realizar análisis, devuelve ok=0 si tuvo éxito

if ok != 0:  # El análisis no fue exitoso
    # Cambiar algunos parámetros de análisis para lograr la convergencia
    # El rendimiento es más lento dentro de este bucle
    ok = 0
    controlTime = ops.getTime()
    while controlTime < TmaxAnalysis and ok == 0:
        controlTime = ops.getTime()
        ok = ops.analyze(1, DtAnalysis)
        if ok != 0:
            print("Trying Newton with Initial Tangent ..")
            ops.test("NormDispIncr", Tol, 1000, 0)
            ops.algorithm("Newton", "-initial")
            ok = ops.analyze(1, DtAnalysis)
            ops.test(testTypeDynamic, TolDynamic, maxNumIterDynamic, 0)
            ops.algorithm(algorithmTypeDynamic)
        if ok != 0:
            print("Trying Broyden ..")
            ops.algorithm("Broyden", 8)
            ok = ops.analyze(1, DtAnalysis)
            ops.algorithm(algorithmTypeDynamic)
        if ok != 0:
            print("Trying NewtonWithLineSearch ..")
            ops.algorithm("NewtonLineSearch", 0.8)
            ok = ops.analyze(1, DtAnalysis)
            ops.algorithm(algorithmTypeDynamic)

print(f"Ground Motion Done. End Time: {ops.getTime()}")