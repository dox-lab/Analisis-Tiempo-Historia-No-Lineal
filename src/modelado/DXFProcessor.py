import ezdxf

def DXFProcessor(title_file):
    doc = ezdxf.readfile(title_file)
    msp = doc.modelspace()

    points = []
    lines = []
    polines = []
    Elems = []
    n_elems_init = 0
    for e in msp:
        if e.dxftype() == "POINT":
            v = list(e.dxf.location)
            points.append([e.dxf.layer, v[0], v[1], v[2]])
        if e.dxftype() == "LINE":
            lines.append([e.dxf.layer, list(e.dxf.start), list(e.dxf.end)])
        if e.dxftype() == "LWPOLYLINE":
            vpoli = list(e.vertices_in_wcs())
            polines.append([e.dxf.layer, list(vpoli[0]), list(vpoli[1]), list(vpoli[2]), list(vpoli[3])])

    Nodes = []
    tag_nodes = 0

    for p in points:
        if p[0].split("_")[2] == "centrales":
            Nodes.append([tag_nodes, p[1], p[2], p[3], 1])
            tag_nodes += 1
        elif p[0].split("_")[2] == "esquinero":
            Nodes.append([tag_nodes, p[1], p[2], p[3], 0.5])
            tag_nodes += 1
        elif p[0].split("_")[2] == "externo":
            Nodes.append([tag_nodes, p[1], p[2], p[3], 0.25])
            tag_nodes += 1

    for li in lines:
        li[1][0] = round(li[1][0], 3)
        li[1][1] = round(li[1][1], 3)
        li[1][2] = round(li[1][2], 3)
        li[2][0] = round(li[2][0], 3)
        li[2][1] = round(li[2][1], 3)
        li[2][2] = round(li[2][2], 3) 
        for nd in Nodes:
            coord = [round(float(nd[1]), 3), round(float(nd[2]), 3), round(float(nd[3]), 3)]
            if li[1] == coord:
                li[1] = int(nd[0])
            elif li[2] == coord:
                li[2] = int(nd[0])

    for pol in polines:
        pol[1][0] = round(pol[1][0], 3)
        pol[1][1] = round(pol[1][1], 3)
        pol[1][2] = round(pol[1][2], 3)
        pol[2][0] = round(pol[2][0], 3)
        pol[2][1] = round(pol[2][1], 3)
        pol[2][2] = round(pol[2][2], 3)
        pol[3][0] = round(pol[3][0], 3)
        pol[3][1] = round(pol[3][1], 3)
        pol[3][2] = round(pol[3][2], 3)
        pol[4][0] = round(pol[4][0], 3)
        pol[4][1] = round(pol[4][1], 3)
        pol[4][2] = round(pol[4][2], 3) 
        for nd in Nodes:
            coord = [round(float(nd[1]), 3), round(float(nd[2]), 3), round(float(nd[3]), 3)]
            if pol[1] == coord:
                pol[1] = int(nd[0])
            elif pol[2] == coord:
                pol[2] = int(nd[0])
            elif pol[3] == coord:
                pol[3] = int(nd[0])
            elif pol[4] == coord:
                pol[4] = int(nd[0])

    Col = []
    Vig = []
    for l in lines:
        if l[0].split("_")[1] == "C":
            Col.append([float(l[0].split("_")[2]), float(l[0].split("_")[3])])
            Elems.append([n_elems_init, l[1], l[2], float(l[0].split("_")[2]), float(l[0].split("_")[3]), 1, 1])
            n_elems_init = n_elems_init +1
        if l[0].split("_")[1] == "V":
            Vig.append([float(l[0].split("_")[2]), float(l[0].split("_")[3])])
            Elems.append([n_elems_init, l[1], l[2], float(l[0].split("_")[2]), float(l[0].split("_")[3]), 2, 2])
            n_elems_init = n_elems_init +1

    nsec_vig = len(Vig) 
    nsec_col = len(Col)

    for pl in polines:
        if pl[0].split("_")[1] == "Mx":
            Elems.append([n_elems_init, pl[1], pl[2], pl[3], pl[4], 3, float(pl[0].split("_")[2])])
            n_elems_init = n_elems_init +1
        if pl[0].split("_")[1] == "My":
            Elems.append([n_elems_init, pl[1], pl[2], pl[3], pl[4], 4, float(pl[0].split("_")[2])])
            n_elems_init = n_elems_init +1

    for node in Nodes:
        node[1:4] = [abs(round(coord, 2)) for coord in node[1:4]]

    return {
        "Nodes": Nodes,
        "Elems": Elems,
        "sec_vig": Vig,
        "sec_col": Col,
        "nsec_sec_vig": nsec_vig,
        "nsec_sec_col": nsec_col
    }
