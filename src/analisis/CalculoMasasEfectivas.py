import numpy as np
import pandas as pd
import openseespy.opensees as ops

class CalculoMasasEfectivas:
    """
    Clase para calcular las masas efectivas del sistema estructural a partir de los modos de vibración.
    """
    def __init__(self, n_modos, n_niveles, masa, periodos):
        """
        Inicializa los parámetros necesarios para el cálculo de las masas efectivas.
        :param n_modos: Número de modos a considerar.
        :param n_niveles: Número de niveles en la estructura.
        :param masa: Matriz de masas del sistema estructural.
        :param periodos: Lista de periodos de los modos de vibración.
        """
        self.n_modos = n_modos  # Número de modos
        self.n_niveles = n_niveles  # Número de niveles
        self.masa = masa  # Matriz de masas
        self.periodos = periodos  # Periodos de los modos de vibración

    def obtener_modos(self, nodos):
        """
        Obtiene los vectores de modos de vibración del sistema estructural.
        :param nodos: Lista de etiquetas de los nodos.
        :return: Matriz con los modos de vibración.
        """
        modos = np.zeros((self.n_modos, 3 * self.n_niveles))  # Matriz para almacenar los modos
        for j in range(1, self.n_modos + 1):
            idx = 0
            # Itera sobre los últimos niveles para obtener los vectores modales
            for nodo in nodos[-self.n_niveles:]:
                temp = ops.nodeEigenvector(nodo, j)  # Función de OpenSees para obtener el vector modal
                modos[j-1, [idx, idx+1, idx+2]] = temp[0], temp[1], temp[-1]
                idx += 3
        return modos

    def calcular_masas(self, nodos):
        """
        Calcula las masas efectivas del sistema estructural.
        :param nodos: Lista de etiquetas de los nodos.
        :return: DataFrame con los resultados de las masas efectivas y el número mínimo de modos.
        """
        modos = self.obtener_modos(nodos)  # Obtener los modos de vibración

        # Definir vectores unitarios en X, Y, y en rotación Z
        vec_x, vec_y, vec_rz = np.zeros(3 * self.n_niveles), np.zeros(3 * self.n_niveles), np.zeros(3 * self.n_niveles)
        vec_x[0::3], vec_y[1::3], vec_rz[2::3] = 1, 1, 1  # Define los vectores unitarios en X, Y y rotación Z

        sum_x, sum_y, sum_rz = 0., 0., 0.  # Inicializar las sumas de masas efectivas
        min_modos = 0  # Número mínimo de modos para alcanzar el 90% de la masa

        # Masas totales en X, Y y rotación Z
        masa_x = np.sum(self.masa[0::3, 0::3])
        masa_y = np.sum(self.masa[1::3, 1::3])
        masa_rz = np.sum(self.masa[2::3, 2::3])

        # Lista para almacenar los resultados
        resultados = []

        # Cálculo de las masas efectivas para cada modo
        for j in range(1, self.n_modos + 1):
            # Producto vectorial de cada modo con la matriz de masas
            fx = modos[j-1].T @ self.masa @ vec_x
            fy = modos[j-1].T @ self.masa @ vec_y
            frz = modos[j-1].T @ self.masa @ vec_rz

            # Sumas acumulativas de masas efectivas
            sum_x += (fx ** 2) / masa_x
            sum_y += (fy ** 2) / masa_y
            sum_rz += (frz ** 2) / masa_rz

            # Determinar el número mínimo de modos a considerar para el 90% de masa efectiva
            if min(sum_x, sum_y, sum_rz) >= 0.90 and min_modos == 0:
                min_modos = j

            # Agregar los resultados de este modo a la lista
            resultados.append({
                'Modo': j,
                'Periodo(s)': self.periodos[j-1],
                'SumX': sum_x,
                'SumY': sum_y,
                'SumRz': sum_rz
            })

        # Crear un DataFrame con los resultados de las masas efectivas
        Masas_Efec = pd.DataFrame(resultados)

        # Imprimir los resultados
        print(f'\nMASAS EFECTIVAS')
        print(Masas_Efec.round(5).to_string(index=False))
        print(f'\nN° mínimo de Modos a considerar:', min_modos)

        return Masas_Efec, min_modos  # Retornar el DataFrame y el número mínimo de modos