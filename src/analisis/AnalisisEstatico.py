import numpy as np
import pandas as pd
import openseespy.opensees as ops

class AnalisisEstatico:
    """
    Clase para realizar el análisis estático de una estructura, aplicando cargas sísmicas y calculando desplazamientos, fuerzas y drifts.
    """
    def __init__(self, modos, Z=0.45, U=1.5, S=1.0, Tp=0.4, Tl=2.5, R=8):
        """
        Inicializa los parámetros del análisis.
        :param modos: Lista de modos de vibración.
        :param Z: Factor de zona sísmica.
        :param U: Factor de uso de la estructura.
        :param S: Factor de suelo.
        :param Tp: Periodo corto.
        :param Tl: Periodo largo.
        :param R: Factor de reducción sísmica.
        """
        self.modos = modos  # Modos de vibración
        self.Z = Z          # Factor de zona sísmica
        self.U = U          # Factor de uso
        self.S = S          # Factor de suelo
        self.Tp = Tp        # Periodo corto
        self.Tl = Tl        # Periodo largo
        self.R = R          # Factor de reducción

    def calc_espectro(self):
        """
        Calcula el espectro de diseño según la norma E030.
        :return: Espectro E030 como un array numpy.
        """
        espectro = np.zeros(len(self.modos))  # Array para almacenar los valores del espectro
        for i, T in enumerate(self.modos):
            if T < 0:
                raise ValueError("El periodo no puede ser negativo!")
            if T < 0.2 * self.Tp:
                espectro[i] = 2.5
            elif T < self.Tp:
                espectro[i] = 2.5
            elif T < self.Tl:
                espectro[i] = 2.5 * (self.Tp / T)
            else:
                espectro[i] = 2.5 * (self.Tp * self.Tl / T**2)
        return espectro * self.Z * self.U * self.S / self.R

    def calc_cargas(self, coef_sismo, pesos, alturas):
        """
        Calcula las cargas estáticas para cada nivel basado en el coeficiente sísmico y el peso por nivel.
        :param coef_sismo: Coeficiente sísmico.
        :param pesos: Peso por nivel.
        :param alturas: Altura por nivel.
        :return: Array de fuerzas por nivel y el exponente k utilizado en el cálculo.
        """
        V_total = coef_sismo * np.sum(pesos)  # Carga total sísmica
        # Definir el exponente k según el primer modo de vibración
        k = 1.0 if self.modos[0] <= 0.5 else 0.75 + 0.5 * self.modos[0]

        divisor = np.sum(pesos * alturas**k)  # Suma ponderada por alturas
        fuerzas = pesos * alturas**k / divisor * V_total  # Cálculo de las fuerzas por nivel
        return fuerzas, k

    def aplicar_y_analizar(self, fuerzas, diaf, nx, ny, nz, dx, dy, dz, Ro, dir='x'):
        """
        Aplica las fuerzas estáticas a los nodos y ejecuta el análisis.
        :param fuerzas: Array de fuerzas por nivel.
        :param diaf: Nodos de los diafragmas.
        :param nx: Número de divisiones en X.
        :param ny: Número de divisiones en Y.
        :param nz: Número de niveles en Z.
        :param dx: Distancia entre columnas en X.
        :param dy: Distancia entre columnas en Y.
        :param dz: Altura entre niveles en Z.
        :param Ro: Factor de sobrerresistencia.
        :param dir: Dirección del análisis ('x' o 'y').
        :return: Corte acumulado y un DataFrame con los resultados del análisis estático.
        """
        # Configurar las cargas estáticas
        ops.loadConst('-time', 0.0)
        ops.remove('timeSeries', 1)
        ops.remove('loadPattern', 1)
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)

        # Aplicar las fuerzas en los nodos según la dirección (X o Y)
        for i, fuerza in enumerate(fuerzas):
            if dir == 'x':
                brazo_x = ny * dy * 0.05
                ops.load(int(diaf[i][0]), fuerza, 0., 0., 0., 0., fuerza * brazo_x)
            elif dir == 'y':
                brazo_y = nx * dx * 0.05
                ops.load(int(diaf[i][0]), 0., fuerza, 0., 0., 0., fuerza * brazo_y)

        # Configurar y ejecutar el análisis estático
        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('Plain')
        ops.system('FullGeneral')
        ops.algorithm('Linear')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')
        ops.analyze(1)

        # Calcular esfuerzos acumulados en los niveles
        corte_acum = np.cumsum(fuerzas[::-1])[::-1]

        # Recopilar los resultados de desplazamientos y drifts
        resultados = []
        tempX, tempY = 0., 0.  # Variables temporales para los drifts

        for i in range(nz):
            desX = ops.nodeDisp(int(diaf[i][0]), 1)  # Desplazamiento en X
            desY = ops.nodeDisp(int(diaf[i][0]), 2)  # Desplazamiento en Y
            rotZ = ops.nodeDisp(int(diaf[i][0]), 6)  # Rotación en Z
            desX += abs(rotZ * ny * dy / 2)  # Ajustar por rotación
            desY += abs(rotZ * nx * dx / 2)
            desX, desY = desX * 0.75 * Ro, desY * 0.75 * Ro  # Escalar con Ro
            driftX = 1000. * (desX - tempX) / dz  # Calcular drift en X
            driftY = 1000. * (desY - tempY) / dz  # Calcular drift en Y
            tempX, tempY = desX, desY  # Actualizar valores temporales

            # Agregar resultados de nivel
            fila = {
                'Nivel': i + 1,
                'FuerzaCortante(kN)': corte_acum[i] / 1000,
                'DesplazamientoX(cm)': desX * 100,
                'DesplazamientoY(cm)': desY * 100,
                'DriftX(‰)': driftX,
                'DriftY(‰)': driftY
            }
            resultados.append(fila)

        # Crear un DataFrame con los resultados del análisis
        result = pd.DataFrame(resultados)

        # Mostrar resultados del análisis
        print(f'\nANÁLISIS ESTÁTICO EN {dir.upper()}')
        print(result.round(4).to_string(index=False))

        return corte_acum, result