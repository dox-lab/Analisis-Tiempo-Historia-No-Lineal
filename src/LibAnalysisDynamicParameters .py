import openseespy.opensees as ops
# --------------------------------------------------------------------------------------------------
# Parámetros de análisis dinámico
# Estoy configurando todas estas variables como variables globales
# para que puedan ser accedidas por un procedimiento
#                                 Silvia Mazzoni & Frank McKenna, 2006


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
