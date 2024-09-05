
import ezdxf

def function_dxf(title_file):
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
    Diap = []
    tag_diag = 1000
    Rest = []

    for p in points:
        if p[0].split("_")[2] == "cent" or p[0].split("_")[2] == "restr":
            Nodes.append([tag_nodes, p[1], p[2], p[3]])
            tag_nodes = tag_nodes +1
        elif p[0].split("_")[2] == "CM":
            Diap.append([tag_diag, p[1], p[2], p[3]])
            tag_diag = tag_diag +1

    for p in Nodes:
        if p[3] == 0:
            Rest.append(p)

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

    Mx = []
    My = []
    for pl in polines:
        if pl[0].split("_")[1] == "Mx":
            Mx.append([float(pl[0].split("_")[2]), pl[1], pl[2], pl[3], pl[4]])
            Elems.append([n_elems_init, pl[1], pl[2], pl[3], pl[4], 3, float(pl[0].split("_")[2])])
            n_elems_init = n_elems_init +1
        if pl[0].split("_")[1] == "My":
            My.append([float(pl[0].split("_")[2]), pl[1], pl[2], pl[3], pl[4]])
            Elems.append([n_elems_init, pl[1], pl[2], pl[3], pl[4], 4, float(pl[0].split("_")[2])])
            n_elems_init = n_elems_init +1

    nz = len(Diap)
    dz = Diap[1][3] - Diap[0][3]

    return Nodes, Elems, Diap, Vig, Col, Mx, My, nsec_vig, nsec_col, nz, dz


Nodes, Elems, Diap, Vig, Col, Mx, My, nsec_vig, nsec_col, nz, dz= function_dxf("EDIFICIO_LP_v4.dxf")

for elem in Elems:
    print(elem)