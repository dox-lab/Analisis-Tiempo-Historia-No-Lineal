from modelado.DXFProcessor import DXFProcessor
from modelado.SectionCreator import StructureModel
from modelado.SeccionesGeometricas import SeccionRectangular, Muro
from analisis.ReadSMDFile import ReadSMDFile
import openseespy.opensees as ops
import math
import matplotlib.pyplot as plt
import numpy as np
import vfo.vfo as vfo
from parametros import *

# Set up ground-motion parameters
GMdirection = 1  # ground-motion direction
GMfile = "H-e12140"  # ground-motion filenames
GMdir = "data/sismos"  # ground-motion directory
GMfact = 1.5  # ground-motion scaling factor

# Set up ground-motion-analysis parameters
DtAnalysis = 0.01  # time-step for lateral analysis
TmaxAnalysis = 10.0  # maximum duration of ground-motion analysis

#  --------------------------------- Perform Dynamic Ground-Motion Analysis
IDloadTag = 400  # for uniformSupport excitation
inFile = f"{GMdir}/{GMfile}.at2"
outFile = f"{GMdir}/{GMfile}.g3"  # Output file

# Definir el valor de dt, que será actualizado por la función
dt = None
# Read the ground-motion file (This would be a custom function to handle file reading and formatting)
ReadSMDFile(inFile, outFile, dt)  # Ensure this function is implemented in Python
