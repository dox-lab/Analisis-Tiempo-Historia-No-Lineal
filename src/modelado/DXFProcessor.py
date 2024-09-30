import ezdxf

def DXFProcessor(title_file):
    # Leer archivo DXF y obtener el espacio del modelo
    doc = ezdxf.readfile(title_file)
    msp = doc.modelspace()

    # Inicializar listas para almacenar diferentes tipos de datos
    points, lines, polylines, elements, diaphragms = [], [], [], [], []
    nodes, rest = [], []
    tag_nodes = 0
    n_elems_init = 0

    # Procesar entidades del DXF
    for entity in msp:
        if entity.dxftype() == "POINT":
            location = list(entity.dxf.location)
            points.append([entity.dxf.layer, *location])
        elif entity.dxftype() == "LINE":
            lines.append([entity.dxf.layer, list(entity.dxf.start), list(entity.dxf.end)])
        elif entity.dxftype() == "LWPOLYLINE":
            vertices = list(entity.vertices_in_wcs())
            polylines.append([entity.dxf.layer] + [list(vertex) for vertex in vertices])

    # Crear nodos a partir de puntos
    for point in points:
        # Identificar tipo de nodo a partir del nombre de la capa
        node_type = point[0].split("_")[2]
        node_mass = {"centrales": 1, "esquinero": 0.5, "externo": 0.25}.get(node_type, 1)
        
        # Agregar nodo a la lista
        nodes.append([tag_nodes, point[1], point[2], point[3], node_mass])
        tag_nodes += 1

    # Función para redondear coordenadas de una lista
    def round_coordinates(coords):
        return [round(coord, 3) for coord in coords]

    # Asignar nodos a líneas
    for line in lines:
        line[1] = round_coordinates(line[1])
        line[2] = round_coordinates(line[2])
        
        for node in nodes:
            coord = round_coordinates(node[1:4])
            if line[1] == coord:
                line[1] = int(node[0])
            elif line[2] == coord:
                line[2] = int(node[0])

    # Asignar nodos a polilíneas
    for polyline in polylines:
        for i in range(1, 5):
            polyline[i] = round_coordinates(polyline[i])
        
        for node in nodes:
            coord = round_coordinates(node[1:4])
            for i in range(1, 5):
                if polyline[i] == coord:
                    polyline[i] = int(node[0])

    # Procesar elementos de columnas y vigas
    columns, beams = [], []
    for line in lines:
        element_type = line[0].split("_")[1]
        dimensions = [float(line[0].split("_")[2]), float(line[0].split("_")[3])]

        if element_type == "C":
            columns.append(dimensions)
            elements.append([n_elems_init, line[1], line[2], *dimensions, 1, 1])
        elif element_type == "V":
            beams.append(dimensions)
            elements.append([n_elems_init, line[1], line[2], *dimensions, 2, 2])
        
        n_elems_init += 1

    # Procesar elementos de muros
    for polyline in polylines:
        wall_type = polyline[0].split("_")[1]
        if wall_type in ["Mx", "My"]:
            elements.append([n_elems_init, polyline[1], polyline[2], polyline[3], polyline[4], 3 if wall_type == "Mx" else 4, float(polyline[0].split("_")[2])])
            n_elems_init += 1

    # Redondear y obtener valores absolutos de las coordenadas de los nodos
    for node in nodes:
        node[1:4] = [abs(round(coord, 2)) for coord in node[1:4]]

    # Identificar y guardar nodos en el nivel cero
    rest = [node for node in nodes if node[3] == 0]

    # Calcular el centro de masa para cada nivel
    levels = sorted(set(node[3] for node in nodes))  # Obtener niveles únicos (coordenada Z)
    for i, level in enumerate(levels):
        level_nodes = [node for node in nodes if node[3] == level]
        total_mass = sum(node[4] for node in level_nodes)
        x_cm = sum(node[1] * node[4] for node in level_nodes) / total_mass
        y_cm = sum(node[2] * node[4] for node in level_nodes) / total_mass
        diaphragms.append([1000 + i, round(x_cm, 2), round(y_cm, 2), level])

    diaphragms.pop(0)  # Eliminar el primer elemento de la lista de diafragmas

    # Calcular el área de planta
    def calcular_area_planta(nodes):
        # Usar el método del polígono (Shoelace formula) para calcular el área
        level_nodes = sorted(nodes, key=lambda x: (x[1], x[2]))  # Ordenar nodos por coordenadas
        x_coords = [node[1] for node in level_nodes]
        y_coords = [node[2] for node in level_nodes]
        
        area = 0.5 * abs(sum(x_coords[i] * y_coords[i+1] - y_coords[i] * x_coords[i+1] for i in range(-1, len(level_nodes)-1)))
        return round(area, 2)

    area_planta = calcular_area_planta([node for node in nodes if node[3] == levels[0]])  # Área del nivel base

    # Retornar los resultados como un diccionario
    return {
        "Nodes": nodes,
        "Elems": elements,
        "sec_vig": beams,
        "sec_col": columns,
        "nsec_vig": len(beams),
        "nsec_col": len(columns),
        "Diap": diaphragms,
        "AreaPlanta": area_planta,
        "Rest_node": rest,
        "Levels": levels
    }