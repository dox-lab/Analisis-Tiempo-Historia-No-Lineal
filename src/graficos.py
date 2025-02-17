import vfo.vfo as vfo

# Definición de la ruta del modelo
model_path = "resultados/3D_Building"

# Graficar el modelo 3D con etiquetas de nodos y elementos
vfo.plot_model(model=model_path, show_nodes="yes", show_eletags="yes", line_width=5)

# Graficar los primeros tres modos de vibración
for mode in range(1, 4):
    vfo.plot_modeshape(model=model_path, modenumber=mode, scale=2500, contour="x", line_width=3)

# Graficar la forma deformada bajo la carga gravitacional
vfo.plot_deformedshape(model=model_path, loadcase="Gravity", scale=1000, line_width=3)

# Descomentar la siguiente línea para animar la forma deformada bajo una carga dinámica
# ani = vfo.animate_deformedshape(model=model_path, loadcase="Dynamic_GM1", line_width=1.5, scale=1000, moviename="Building_Dynamic_GM1")