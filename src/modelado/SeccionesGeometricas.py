class SeccionesGeometricas:
    """
    Clase base para calcular las propiedades geométricas de diferentes tipos de secciones.
    """
    def __init__(self):
        pass
    
    def calcular_area(self):
        """
        Método para calcular el área de la sección.
        Se implementa en las clases hijas.
        """
        raise NotImplementedError("Este método debe ser implementado en la clase hija.")
    
    def calcular_momentos_inercia(self):
        """
        Método para calcular los momentos de inercia de la sección.
        Se implementa en las clases hijas.
        """
        raise NotImplementedError("Este método debe ser implementado en la clase hija.")
    
    def calcular_momento_polar(self):
        """
        Método para calcular el momento polar de la sección.
        Se implementa en las clases hijas.
        """
        raise NotImplementedError("Este método debe ser implementado en la clase hija.")

class SeccionRectangular(SeccionesGeometricas):
    """
    Clase para calcular las propiedades geométricas de una sección rectangular.
    """
    def __init__(self, b, h):
        """
        Inicializa la sección rectangular con base (b) y altura (h).
        :param b: base de la sección en metros
        :param h: altura de la sección en metros
        """
        self.b = b
        self.h = h
    
    def calcular_area(self):
        """
        Calcula el área de la sección rectangular.
        :return: Área de la sección (m^2)
        """
        return self.b * self.h
    
    def calcular_momentos_inercia(self):
        """
        Calcula los momentos de inercia de la sección rectangular respecto a los ejes Z e Y.
        :return: Momento de inercia Iz (m^4), Momento de inercia Iy (m^4)
        """
        Iz = self.b * self.h**3 / 12
        Iy = self.b**3 * self.h / 12
        return Iz, Iy
    
    def calcular_momento_polar(self):
        """
        Calcula el momento polar de inercia de la sección rectangular.
        :return: Momento polar de inercia Jxx (m^4)
        """
        aa, bb = max(self.b, self.h), min(self.b, self.h)
        beta = 1 / 3 - 0.21 * bb / aa * (1 - (bb / aa)**4 / 12)
        Jxx = beta * bb**3 * aa
        return Jxx
