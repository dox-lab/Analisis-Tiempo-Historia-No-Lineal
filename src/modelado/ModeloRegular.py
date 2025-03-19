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
        # Definición de dimensiones y divisiones en las tres direcciones
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # Cálculo de la longitud total en las direcciones X y Y
        self.Lx = self.nx * self.dx  # Longitud total en X
        self.Ly = self.ny * self.dy  # Longitud total en Y

    def calcular_factor(self, i, j, k):
        """
        Calcula el factor del nodo basado en el número de nodos adyacentes en su plano.
        :param i: Índice en la dirección X.
        :param j: Índice en la dirección Y.
        :param k: Índice en la dirección Z.
        :return: Factor del nodo (0.25 para nodos esquineros, 0.5 para fronteras, 1 para nodos centrales).
        """
        adyacentes = 0
        
        # Contar el número de nodos adyacentes en el plano XY para un nivel k dado
        if i > 0:  # Nodo adyacente a la izquierda
            adyacentes += 1
        if i < self.nx:  # Nodo adyacente a la derecha
            adyacentes += 1
        if j > 0:  # Nodo adyacente abajo
            adyacentes += 1
        if j < self.ny:  # Nodo adyacente arriba
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
        Cada nodo contiene su identificador, coordenadas (x, y, z) y su factor.
        :return: Lista de nodos con formato [nodeID, x, y, z, factor].
        """
        Nodes = []
        for i in range(self.nx + 1):
            for j in range(self.ny + 1):
                for k in range(self.nz + 1):
                    # Calcular el identificador del nodo (nodeID)
                    nodeID = (i + 1) * 100 + (j + 1) * 10 + k
                    
                    # Calcular las coordenadas (x, y, z) del nodo
                    x, y, z = i * self.dx, j * self.dy, k * self.dz
                    
                    # Calcular el factor del nodo
                    factor = self.calcular_factor(i, j, k)
                    
                    # Añadir nodo a la lista
                    Nodes.append([nodeID, x, y, z, factor])
        return Nodes

    def generar_elementos_vigas_columnas(self):
        """
        Genera los elementos (vigas y columnas) y devuelve la lista de elementos.
        :return: Lista de elementos en formato [elementID, nodo_inicio, nodo_final, tipo].
        Tipo 1: Columna (dirección Z), Tipo 2: Viga (direcciones X o Y).
        """
        Elems = []
        for i in range(self.nx + 1):
            for j in range(self.ny + 1):
                for k in range(self.nz + 1):
                    nodeID = (i + 1) * 100 + (j + 1) * 10 + k

                    # Crear columnas en la dirección Z
                    if k > 0:  # Solo si no es el primer nivel (Z > 0)
                        Elems.append([f'{nodeID - 1}{nodeID}', nodeID - 1, nodeID, 1])  # Columna
                    
                    # Crear vigas en la dirección X
                    if i > 0 and k != 0:  # Solo si no es el primer nodo en X y no es el primer nivel (Z > 0)
                        Elems.append([f'{nodeID - 100}{nodeID}', nodeID - 100, nodeID, 2])  # Viga en X
                    
                    # Crear vigas en la dirección Y
                    if j > 0 and k != 0:  # Solo si no es el primer nodo en Y y no es el primer nivel (Z > 0)
                        Elems.append([f'{nodeID - 10}{nodeID}', nodeID - 10, nodeID, 2])  # Viga en Y
        return Elems

    def generar_diaphragmas(self):
        """
        Genera los centros de diafragmas rígidos y devuelve la lista de diafragmas.
        Cada diafragma se ubica en el centro del plano X-Y para cada nivel Z.
        :return: Lista de diafragmas en formato [ID, centro_x, centro_y, z].
        """
        Diap = [[i + 1001, self.Lx / 2.0, self.Ly / 2.0, self.dz * (i + 1)] for i in range(self.nz)]
        return Diap
