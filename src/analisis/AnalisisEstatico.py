import numpy as np
import pandas as pd
import openseespy.opensees as ops

class AnalisisEstatico:
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
        self.modos = modos
        self.Z = Z
        self.U = U
        self.S = S
        self.Tp = Tp
        self.Tl = Tl
        self.R = R

    def calc_espectro(self):
        """
        Calcula el espectro de diseño según la norma E030.
        :return: Espectro E030.
        """
        espectro = np.zeros(len(self.modos))
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
        Calcula las cargas estáticas para cada nivel.
        :param coef_sismo: Coeficiente sísmico.
        :param pesos: Peso por nivel.
        :param alturas: Altura por nivel.
        :return: Fuerza por nivel y exponente k.
        """
        V_total = coef_sismo * np.sum(pesos)
        if self.modos[0] <= 0.5:
            k = 1.0
        else:
            k = 0.75 + 0.5 * self.modos[0]

        divisor = np.sum(pesos * alturas**k)
        fuerzas = pesos * alturas**k / divisor * V_total
        return fuerzas, k

    def aplicar_y_analizar(self, fuerzas, diaf, nx, ny, nz, dx, dy, dz, Ro, dir='x'):
        """
        Aplica las fuerzas estáticas a los nodos y ejecuta el análisis.
        :param fuerzas: Fuerzas por nivel.
        :param diaf: Nodos de los diafragmas.
        :param nz: Número de niveles.
        :param dx: Distancia entre columnas en X.
        :param dy: Distancia entre columnas en Y.
        :param dz: Altura entre niveles.
        :param Ro: Factor de sobrerresistencia.
        :param dir: Dirección del análisis ('x' o 'y').
        :return: DataFrame con resultados del análisis estático.
        """
        # Configurar las cargas
        ops.loadConst('-time', 0.0)
        ops.remove('timeSeries', 1)
        ops.remove('loadPattern', 1)
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        
        # Aplicar las fuerzas en función de la dirección
        for i, fuerza in enumerate(fuerzas):
            if dir == 'x':
                brazo_x = ny * dy * 0.05                
                ops.load(int(diaf[i][0]), fuerza, 0., 0., 0., 0., fuerza * brazo_x)
            elif dir == 'y':
                brazo_y = nx * dx * 0.05
                ops.load(int(diaf[i][0]), 0., fuerza, 0., 0., 0., fuerza * brazo_y)

        # Configurar el análisis
        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('Plain')
        ops.system('FullGeneral')
        ops.algorithm('Linear')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')
        ops.analyze(1)

        # Calcular esfuerzos acumulados
        corte_acum = np.cumsum(fuerzas[::-1])[::-1]

        # Recopilar resultados
        resultados = []
        tempX, tempY = 0., 0.  # Desplazamiento temporal para el drift

        for i in range(nz):
            desX = ops.nodeDisp(int(diaf[i][0]), 1)
            desY = ops.nodeDisp(int(diaf[i][0]), 2)
            rotZ = ops.nodeDisp(int(diaf[i][0]), 6)
            desX += abs(rotZ * ny * dy / 2)
            desY += abs(rotZ * nx * dx / 2)
            desX, desY = desX * 0.75 * Ro, desY * 0.75 * Ro
            driftX = 1000. * (desX - tempX) / dz
            driftY = 1000. * (desY - tempY) / dz
            tempX, tempY = desX, desY

            # Agregar la fila como diccionario a la lista de filas
            fila = {
                'Nivel': i + 1,
                'FuerzaCortante(kN)': corte_acum[i] / 1000,
                'DesplazamientoX(cm)': desX * 100,
                'DesplazamientoY(cm)': desY * 100,
                'DriftX(‰)': driftX,
                'DriftY(‰)': driftY
            }
            resultados.append(fila)

        # Crear el DataFrame directamente desde la lista de resultados
        result = pd.DataFrame(resultados)

        print(f'\nANÁLISIS ESTÁTICO EN {dir.upper()}')
        print(result.round(4).to_string(index=False))
        return result
