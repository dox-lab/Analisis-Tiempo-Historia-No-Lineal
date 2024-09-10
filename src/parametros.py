from unidades import *
import numpy as np
import ezdxf

# PROPIEDADES DE MATERIALES Y SECCIONES
# fc = 210*kgf/cm**2
# E = 15100*(fc/(kgf/cm**2))**0.5*kgf/cm**2
# G = 0.5*E/(1+0.2)

fy = 4200*kgf/cm**2  # Fluencia del acero
Es = 2.1*10**6*kgf/cm**2  # Módulo de elasticidad del acero
fc = 210  # kg/cm2             #Resistencia a la compresión del concreto
E = 15100*fc**0.5*kgf/cm**2  # Módulo de elasticidad del concreto
G = 0.5*E/(1+0.2)  # Módulo de corte del concreto
fc = fc*kgf/cm**2
cover = 4*cm  # Recubrimiento de vigas y columnas

# Parametros no lineales de comportamiento del concreto
fc1 = -fc  # Resistencia a la compresión del concreto
Ec1 = E  # Módulo de elasticidad del concreto
nuc1 = 0.2  # Coeficiente de Poisson
Gc1 = Ec1/(2*(1+nuc1))  # Módulo de corte del concreto

# Concreto confinado
Kfc = 1.0 #1.3               # ratio of confined to unconfined concrete strength
Kres = 0.2                    # ratio of residual/ultimate to maximum stress
fpc1 = Kfc*fc1
epsc01 = 2*fpc1/Ec1
fpcu1 = Kres*fpc1
epsU1 = 5*epsc01  # 20
lambda1 = 0.1
# Concreto no confinado
fpc2 = fc1
epsc02 = -0.003
fpcu2 = Kres*fpc2
epsU2 = -0.006  # -0.01
# Propiedades de resistencia a la tracción
ft1 = -0.14*fpc1
ft2 = -0.14*fpc2
Ets = ft2/0.002
#print(E/10**8, Ets/10**8)

# Densidad del concreto
ρ = 2400*kg/m**3


def prop_vig(bv, hv):

    Av = bv*hv
    ρlv = 2400*Av*m**2  # Densidad lineal de la viga
    Izv = hv**3*bv/12
    Iyv = hv*bv**3/12
    k = 1/3-0.21*bv/hv*(1-(bv/hv)**4/12)
    Jv = k*bv**3*hv
    return Av, ρlv, Izv, Iyv, k, Jv


def prop_col(bc, hc):

    Ac = bc*hc
    ρlc = 2400*Ac*m**2 	# Densidad lineal de la columna
    Izc = hc**3*bc/12
    Iyc = hc*bc**3/12
    k = 1/3-0.21*bc/hc*(1-(bc/hc)**4/12)
    Jc = k*bc**3*hc
    return Ac, ρlc, Izc, Iyc, k, Jc

def funtion_mx(L, t):
    Lmx = L*m   # Longitud del muro
    t = t*m  # Espesor del muro
    ancho = 50*cm  # Discretización del muro
    mufx = int(round(Lmx/ancho))
    ttx = []
    for i in range(0,mufx):
        ttx.append(t)
    # ttx = np.zeros(mufx) 
    # ttx[:] = t # Arreglo de espesores
    wwx = []
    for i in range(0,mufx):
        wwx.append(Lmx/(mufx))
    # wwx = np.zeros(mufx)
    # wwx[:] = Lmx/(mufx) # Arreglo de anchos
    ρρx = []
    for i in range(0,mufx):
        if i == 0:
            ρρx.append(0.01)
        elif i == mufx:
            ρρx.append(0.01)
        else:
            ρρx.append(0.0064)
    # ρρx = np.zeros(mufx)
    # ρρx[:] = 0.0064  # Cuantía vertical en muros
    # ρρx[0], ρρx[-1] = 0.01, 0.01  # Cuantía vertical en núcleo
    concx = []
    for i in range(0,mufx):
        if i == 0:
            concx.append(int(4))
        elif i == mufx:
            concx.append(int(4))
        else:
            concx.append(int(5))
    # concx = np.zeros(mufx)
    # concx[:] = int(5) # Concreto sin confinar
    # concx[0], concx[-1] = int(4), int(4) # Concreto confinado
    acerox = []
    for i in range(0,mufx):
        acerox.append(int(6))
    # acerox = np.zeros(mufx) 
    # acerox[:] = int(6) # Modelo No Lineal del acero
    return mufx, ttx, wwx, ρρx, concx, acerox

def funtion_my(L, t):
    Lmy = L*m
    t = t*m
    ancho = 50*cm
    mufy = int(round(Lmy/ancho))
    tty = []
    for i in range(0,mufy):
        tty.append(t)
    # tty = np.zeros(mufy)
    # tty[:] = t
    wwy = []
    for i in range(0,mufy):
        wwy.append(Lmy/(mufy))
    # wwy = np.zeros(mufy)
    # wwy[:] = Lmy/(mufy)
    ρρy = []
    for i in range(0,mufy):
        if i == 0:
            ρρy.append(0.01)
        elif i == mufy:
            ρρy.append(0.01)
        else:
            ρρy.append(0.0064)
    # ρρy = np.zeros(mufy)
    # ρρy[:] = 0.0064  # Cuantía vertical en muros
    # ρρy[0], ρρy[-1] = 0.01, 0.01  # Cuantía vertical en núcleo
    concy = []
    for i in range(0,mufy):
        if i == 0:
            concy.append(int(4))
        elif i == mufy:
            concy.append(int(4))
        else:
            concy.append(int(5))
    # concy = np.zeros(mufy)
    # concy[:] = int(5)
    # concy[0], concy[-1] = int(4), int(4)
    aceroy = []
    for i in range(0,mufy):
        aceroy.append(int(6))
    # aceroy = np.zeros(mufy)
    # aceroy[:] = int(6)
    return mufy, tty, wwy, ρρy, concy, aceroy


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
    i, ii = 0, 0
    for l in lines:
        if l[0].split("_")[1] == "C":
            Col.append([float(l[0].split("_")[2]), float(l[0].split("_")[3])])
            Elems.append([n_elems_init, l[1], l[2], float(l[0].split("_")[2]), float(l[0].split("_")[3]), 1, i+1])
            i = i + 1
            n_elems_init = n_elems_init +1
        if l[0].split("_")[1] == "V":
            Vig.append([float(l[0].split("_")[2]), float(l[0].split("_")[3])])
            Elems.append([n_elems_init, l[1], l[2], float(l[0].split("_")[2]), float(l[0].split("_")[3]), 2, ii+1])
            n_elems_init = n_elems_init +1
            ii = ii + 1

    nsec_vig = len(Vig) 
    nsec_col = len(Col)

    for ele in Elems:
        if ele[5] == 2:
            ele[6] =ele[6] + nsec_col

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

    return Nodes, Elems, Diap, Vig, Col, Mx, My, nsec_vig, nsec_col, nz, dz, Rest

# Aplicando Cargas vivas y muertas
wLive = 250*kg/m**2
wLosa = 300*kg/m**2
wAcab = 100*kg/m**2
wTabi = 150*kg/m**2
wTotal = 1.0*(wLosa+wAcab+wTabi)+0.25*wLive

aPlanta = 307.5003