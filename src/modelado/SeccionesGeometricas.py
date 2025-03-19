from unidades import *

class SeccionesGeometricas:
    """
    Clase base para calcular las propiedades geométricas de diferentes tipos de secciones.
    Esta clase actúa como una plantilla que debe ser heredada por otras clases que 
    representen tipos específicos de secciones geométricas.
    """
    def __init__(self):
        pass
    
    def calcular_area(self):
        """
        Método abstracto para calcular el área de la sección.
        Debe ser implementado en las clases hijas.
        """
        raise NotImplementedError("Este método debe ser implementado en la clase hija.")
    
    def calcular_momentos_inercia(self):
        """
        Método abstracto para calcular los momentos de inercia de la sección.
        Debe ser implementado en las clases hijas.
        """
        raise NotImplementedError("Este método debe ser implementado en la clase hija.")
    
    def calcular_momento_polar(self):
        """
        Método abstracto para calcular el momento polar de la sección.
        Debe ser implementado en las clases hijas.
        """
        raise NotImplementedError("Este método debe ser implementado en la clase hija.")

class SeccionRectangular(SeccionesGeometricas):
    """
    Clase para calcular las propiedades geométricas de una sección rectangular.
    Hereda de la clase SeccionesGeometricas.
    """
    def __init__(self, b, h):
        """
        Inicializa la sección rectangular con base (b) y altura (h).
        :param b: base de la sección en metros
        :param h: altura de la sección en metros
        """
        self.b = b  # Base de la sección rectangular
        self.h = h  # Altura de la sección rectangular
    
    def calcular_area(self):
        """
        Calcula el área de la sección rectangular.
        :return: Área de la sección en metros cuadrados (m^2).
        """
        return self.b * self.h  # El área de un rectángulo es base * altura
    
    def calcular_momentos_inercia(self):
        """
        Calcula los momentos de inercia de la sección rectangular respecto a los ejes Z e Y.
        :return: Momento de inercia Iz (respecto al eje Z) en metros a la cuarta potencia (m^4),
                 Momento de inercia Iy (respecto al eje Y) en metros a la cuarta potencia (m^4).
        """
        Iz = self.b * self.h**3 / 12  # Momento de inercia respecto al eje Z
        Iy = self.b**3 * self.h / 12  # Momento de inercia respecto al eje Y
        return Iz, Iy
    
    def calcular_momento_polar(self):
        """
        Calcula el momento polar de inercia de la sección rectangular.
        Este momento mide la resistencia de la sección a la torsión.
        :return: Momento polar de inercia Jxx en metros a la cuarta potencia (m^4).
        """
        aa, bb = max(self.b, self.h), min(self.b, self.h)  # Determina cuál es la mayor y la menor dimensión
        beta = 1 / 3 - 0.21 * bb / aa * (1 - (bb / aa)**4 / 12)  # Factor de corrección beta
        Jxx = beta * bb**3 * aa  # Momento polar calculado usando el factor beta
        return Jxx


class Muro(SeccionesGeometricas):
    """
    Clase para calcular las propiedades geométricas de un muro.
    Hereda de SeccionesGeometricas.
    """
    def __init__(self, l, t):
        self.l = l  # Longitud del muro
        self.t = t  # Espesor del muro

    def calc_area(self):
        return self.l * self.t  # Área del muro

    def calc_inercia(self):
        Iz = self.t * self.l**3 / 12  # Momento de inercia respecto al eje Z
        Iy = self.l * self.t**3 / 12  # Momento de inercia respecto al eje Y
        return Iz, Iy

    def calc_mom_polar(self):
        aa, bb = max(self.l, self.t), min(self.l, self.t)
        beta = 1 / 3 - 0.21 * bb / aa * (1 - (bb / aa)**4 / 12)
        Jxx = beta * bb**3 * aa
        return Jxx

    def prop_muros(self):
        """
        Calcula las propiedades específicas del muro en la dirección X.
        """
        lm = self.l * m  # Longitud del muro en metros
        tm = self.t * m  # Espesor del muro en metros
        w = 50 * cm  # Ancho de discretización
        nf = int(round(lm / w))
        
        # Crear listas para espesores, anchos, cuantías y materiales
        thk = [tm for _ in range(nf)]
        w_list = [lm / nf for _ in range(nf)]
        rho = [0.01 if i == 0 or i == nf - 1 else 0.0064 for i in range(nf)]
        conc = [4 if i == 0 or i == nf - 1 else 5 for i in range(nf)]
        steel = [6 for _ in range(nf)]
        
        return nf, thk, w_list, rho, conc, steel
