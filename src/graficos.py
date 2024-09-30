
import vfo.vfo as vfo
vfo.plot_model(model="src/modelado/3D_Building", show_nodes="yes", show_eletags="yes", line_width=5)
vfo.plot_modeshape(model="src/modelado/3D_Building", modenumber=1, scale=2500,contour="x", line_width=3)
vfo.plot_modeshape(model="src/modelado/3D_Building", modenumber=2, scale=2500,contour="x", line_width=3)
vfo.plot_modeshape(model="src/modelado/3D_Building", modenumber=3, scale=2500,contour="x", line_width=3)
vfo.plot_deformedshape(model="src/modelado/3D_Building", loadcase="Gravity", scale=1000, line_width=3)
ani = vfo.animate_deformedshape(model="src/modelado/3D_Building", loadcase="Dynamic_GM1", line_width=1.5, scale=1000,  moviename="Building_Dynamic_GM1")

# vfo.plot_modeshape(modenumber=1, scale=2500, line_width=3, model="src/modelado/3D_Building")

# opsplt.animate_deformedshape(Model="src/modelado/3D_Building", dt=0.1)
