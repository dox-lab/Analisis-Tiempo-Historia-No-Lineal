import numpy as np
from matplotlib import pyplot as plt

class Espectro:
    """
    Clase padre para representar un espectro de respuesta genérico.
    """
    def __init__(self):
        """
        Inicializa los parámetros comunes que pueden ser utilizados por las clases hijas.
        """
        pass

    def calcular_sa(self, *args, **kwargs):
        """
        Método que debe ser implementado por las clases hijas para calcular el espectro de respuesta.
        """
        raise NotImplementedError("Este método debe ser implementado por las clases hijas")

    def graficar_espectro(self, T, Sa):
        """
        Grafica el espectro de respuesta.
        :param T: Lista de períodos.
        :param Sa: Lista de aceleraciones espectrales.
        """
        plt.figure(figsize=(8, 6))
        plt.plot(T, Sa, label="Espectro de Diseño", color="b")
        plt.title("Espectro de Diseño")
        plt.xlabel("Período (s)")
        plt.ylabel("Aceleración espectral Sa (g)")
        plt.grid(True)
        plt.legend()
        plt.show()


class EspectroE030(Espectro):
    """
    Clase hija que calcula el espectro de respuesta según la norma E.030.
    """
    def __init__(self, Z=0.45, U=1.5, S=1.0, R=1, Tp=0.4, Tl=2.5):
        """
        Inicializa los parámetros específicos de la norma E.030.
        :param Z: Factor de zona sísmica.
        :param U: Factor de uso de la estructura.
        :param S: Factor de tipo de suelo.
        :param R: Factor de reducción sísmica.
        :param Tp: Período límite inferior (en segundos).
        :param Tl: Período límite superior (en segundos).
        """
        super().__init__()
        self.Z = Z
        self.U = U
        self.S = S
        self.R = R
        self.Tp = Tp
        self.Tl = Tl

    def calcular_sa(self, T_max=6, paso=0.05):
        """
        Calcula el espectro de respuesta sísmica según la norma E.030.
        :param T_max: Periodo máximo en segundos (por defecto es 6 segundos).
        :param paso: Paso entre los periodos (por defecto 0.05 segundos).
        :return: Lista de períodos y aceleraciones espectrales (Sa).
        """
        periodos = list(np.arange(0, T_max + paso, paso))
        Sa = []

        for T in periodos:
            if T == 0:
                C = 2.5  # Valor máximo para períodos cercanos a 0, evitando división por cero
            elif T <= self.Tp:
                C = 2.5
            elif self.Tp < T <= self.Tl:
                C = 2.5 * self.Tp / T
            else:
                C = 2.5 * self.Tp * self.Tl / T ** 2

            # Calcular Sa para el período actual usando los factores Z, U, S, R
            Sa_T = (self.Z * self.U * self.S / self.R) * C
            Sa.append(Sa_T)

        return periodos, Sa