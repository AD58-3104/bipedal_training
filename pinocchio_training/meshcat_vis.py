import pinocchio as pin
from pinocchio.visualize import meshcat_visualizer
import numpy as np
import sys

vis = meshcat_visualizer.MeshcatVisualizer
args = ' '.join(sys.argv[1:])
model = meshcat_visualizer.loadMesh(args)
vis.loadViewerModel(model)
vis.initViewer()
vis.display()