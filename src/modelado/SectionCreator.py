import openseespy.opensees as ops
import matplotlib.pyplot as plt
import opsvis as opsv

class Section:
    def __init__(self, tag, seccion, G, cover, material_tag, As):
        """
        :param seccion: Instancia de SeccionRectangular que contiene las propiedades de la geometría.
        """
        self.tag = tag
        self.seccion = seccion
        self.G = G
        self.cover = cover
        self.material_tag = material_tag
        self.As = As
        self.fibers = []

    def create_section(self):
        # Usar el objeto seccion para obtener las propiedades geométricas
        Ac = self.seccion.calcular_area()
        Iz, Iy = self.seccion.calcular_momentos_inercia()
        J = self.seccion.calcular_momento_polar()

        y1 = self.seccion.h / 2.0
        z1 = self.seccion.b / 2.0

        n_y = int((y1) / self.cover * 2)
        n_z = int((z1) / self.cover * 2)

        self.fibers = [['section', 'Fiber', self.tag, '-GJ', self.G * J],
                       ['patch', 'rect', 4, n_y, n_z, self.cover - y1, self.cover - z1, y1 - self.cover, z1 - self.cover],
                       ['patch', 'rect', 5, n_y + 2, 1, -y1, z1 - self.cover, y1, z1],
                       ['patch', 'rect', 5, n_y + 2, 1, -y1, -z1, y1, self.cover - z1],
                       ['patch', 'rect', 5, 1, n_z, -y1, self.cover - z1, self.cover - y1, z1 - self.cover],
                       ['patch', 'rect', 5, 1, n_z, y1 - self.cover, self.cover - z1, y1, z1 - self.cover],
                       ['layer', 'straight', 6, 3, self.As, y1 - self.cover, z1 - self.cover, y1 - self.cover, self.cover - z1],
                       ['layer', 'straight', 6, 2, self.As, 0.0, z1 - self.cover, 0.0, self.cover - z1],
                       ['layer', 'straight', 6, 3, self.As, self.cover - y1, z1 - self.cover, self.cover - y1, self.cover - z1]]

        for li in self.fibers:
            if li[0] == 'section':
                eval('ops.%s("%s",%s,"%s",%s)' % tuple(li))
            else:
                eval('ops.%s("%s",%s,%s,%s,%s,%s,%s,%s)' % tuple(li))

    def plot_section(self):
        matcolor = ['r', 'lightgrey', 'gold', 'r', 'lightgrey', 'gold']
        opsv.plot_fiber_section(self.fibers, matcolor=matcolor)
        plt.axis('equal')
        plt.show()

class StructureModel:
    def __init__(self, nsec_col, nsec_vig, Col, Vig, G, cover):
        """
        Col y Vig son listas de instancias de SeccionRectangular
        """
        self.nsec_col = nsec_col
        self.nsec_vig = nsec_vig
        self.Col = Col
        self.Vig = Vig
        self.G = G
        self.cover = cover

    def create_structure(self):
        for i in range(self.nsec_col):
            cuant = 0.01
            As = self.Col[i].b * self.Col[i].h * cuant / 8
            col_section = Section(i + 1, self.Col[i], self.G, self.cover, 4, As)
            col_section.create_section()
        col_section.plot_section()

        for i in range(self.nsec_vig):
            cuant = 0.01
            As = self.Vig[i].b * self.Vig[i].h * cuant / 8
            vig_section = Section(self.nsec_col + i + 1, self.Vig[i], self.G, self.cover, 5, As)
            vig_section.create_section()
        vig_section.plot_section()
