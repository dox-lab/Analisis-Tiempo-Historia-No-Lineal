class ModeloRegular:
    def __init__(self, dx, dy, dz, nx, ny, nz):
        """
        Inicializa el modelo estructural con los parámetros básicos.
        :param dx: Distancia entre columnas en dirección X.
        :param dy: Distancia entre columnas en dirección Y.
        :param dz: Altura entre pisos (dirección Z).
        :param nx: Número de divisiones en X.
        :param ny: Número de divisiones en Y.
        :param nz: Número de divisiones en Z.
        """
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.Lx = self.nx * self.dx  # Longitud total en X
        self.Ly = self.ny * self.dy  # Longitud total en Y

    def calcular_factor(self, i, j, k):
        """
        Calcula el factor del nodo basado en el número de nodos adyacentes en su plano.
        :param i: índice en la dirección X
        :param j: índice en la dirección Y
        :param k: índice en la dirección Z
        :return: factor (0.25 para nodos esquineros, 0.5 para fronteras, 1 para nodos centrales)
        """
        adyacentes = 0
        
        # Revisar adyacentes en el plano XY (para un k fijo)
        if i > 0:  # Adyacente a la izquierda
            adyacentes += 1
        if i < self.nx:  # Adyacente a la derecha
            adyacentes += 1
        if j > 0:  # Adyacente abajo
            adyacentes += 1
        if j < self.ny:  # Adyacente arriba
            adyacentes += 1

        # Determinar el factor según el número de adyacentes
        if adyacentes == 2:
            return 0.25  # Nodo esquinero
        elif adyacentes == 3:
            return 0.5  # Nodo en frontera
        else:
            return 1.0  # Nodo central

    def generar_nodos(self):
        """
        Genera los nodos del modelo y devuelve la lista de nodos.
        :return: Lista de nodos con [tag, x, y, z, factor].
        """
        Nodes = []
        for i in range(self.nx + 1):
            for j in range(self.ny + 1):
                for k in range(self.nz + 1):
                    nodeID = (i + 1) * 100 + (j + 1) * 10 + k
                    x, y, z = i * self.dx, j * self.dy, k * self.dz
                    factor = self.calcular_factor(i, j, k)
                    Nodes.append([nodeID, x, y, z, factor])
        return Nodes

    def generar_elementos_vigas_columnas(self):
        """
        Genera los elementos (vigas y columnas) y devuelve la lista de elementos.
        :return: Lista de elementos.
        """
        Elems = []
        for i in range(self.nx + 1):
            for j in range(self.ny + 1):
                for k in range(self.nz + 1):
                    nodeID = (i + 1) * 100 + (j + 1) * 10 + k

                    # Crear columnas en la dirección Z
                    if k > 0:  # Elementos verticales (columna)
                        Elems.append([f'{nodeID - 1}{nodeID}', nodeID - 1, nodeID, 1])
                    # Crear vigas en la dirección X
                    if i > 0 and k != 0:  # Elementos horizontales en X (viga)
                        Elems.append([f'{nodeID - 100}{nodeID}', nodeID - 100, nodeID, 2])
                    # Crear vigas en la dirección Y
                    if j > 0 and k != 0:  # Elementos horizontales en Y (viga)
                        Elems.append([f'{nodeID - 10}{nodeID}', nodeID - 10, nodeID, 2])
        return Elems

    def generar_diaphragmas(self):
        """
        Genera los centros de diafragmas rígidos y devuelve la lista de diafragmas.
        :return: Lista de diafragmas.
        """
        Diap = [[i + 1001, self.Lx / 2.0, self.Ly / 2.0, self.dz * (i + 1)] for i in range(self.nz)]
        return Diap
