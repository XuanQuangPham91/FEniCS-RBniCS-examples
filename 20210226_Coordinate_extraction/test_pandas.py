from sys import path
from dolfin import *
from mshr import *
from Python_module_Quang import *

import matplotlib.pyplot as plt
# %matplotlib inline

import sys
import numpy as np
import pandas as pd
from numpy.core.records import array
from numpy.lib.function_base import append

lab_computer = True
if lab_computer:
    # path of lab's computer
    sys.path.append(
        '/media/xuanquang/Gaumap Lab data/05_Git_project/FEniCS-RBniCS-examples/20210218_2D_tangential_load/'
    )
else:
    # path of MSI laptop
    sys.path.append(
        '/home/xuanquang/Project_Git/FEniCS-RBniCS-examples/20210218_2D_tangential_load/'
    )

try:
    import file
except:
    print("fail to import file")

try:
    import tangential_load
except:
    print("fail to import tangential_load.py")

format = "png"

from mpl_toolkits import mplot3d


def cal_magnitude(ux, uy):
    u_magnitude = []
    for i in range(len(ux)):
        norm = sqrt((ux[i])**2 + (uy[i])**2)
        u_magnitude.append(norm)
    return u_magnitude


# import mesh
mesh = Mesh("data/elastic_block.xml")
V = VectorFunctionSpace(mesh, "Lagrange", 1)

u_FE = load_HDF5(V, mesh, title='u_FE')
u_magnitude = cal_u_magnitude(u=u_FE, mesh=mesh)
u_mag = u_magnitude.vector().get_local()

x, y, ux, uy, nodal_values = file.coordinates_operator(V, u_FE)

# u_mag = cal_magnitude(ux, uy)

print(f"len(u_mag): {len(u_mag)}")
print(type(x))
print(type(y))
print(type(u_mag))
print(u_mag)

df = pd.DataFrame({
    'x': array(x),
    'y': array(y),
    'z': u_mag,
})

df.to_csv("solution/sample.csv", index=False)

from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

interpolation_point = 500
xline = np.linspace(-1, 1, interpolation_point)
yline = np.linspace(-1, 1, interpolation_point)
grid_x, grid_y = np.meshgrid(xline, yline)
inter_method = 'cubic'
R_list = np.linspace(0.8, 0.9, 11)  # Choose 11 or 21

for R in R_list:
    FE_data = pd.read_csv('solution/sample.csv')
    FE_data = FE_data.dropna()
    FE_data.columns = ['x_coord', 'y_coord', 'u_mag']
    FE_points = FE_data[['x_coord', 'y_coord']].values
    FE_mag = FE_data[['u_mag']].values

print(FE_mag.shape)

max_AE_list = []


def visualize_abs_error(FE_mag):
    # Visualize absolute error
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    abs_error = FE_mag
    max_AE = max(abs_error)
    min_AE = min(abs_error)
    max_AE_list.append(max_AE)
    vmin_AE, vmax_AE = min(min_AE), max(max_AE)

    # Visualize FE solution
    FE_grid = griddata(FE_points,
                       FE_mag, (grid_x, grid_y),
                       method=inter_method)
    FE_grid = FE_grid.reshape(interpolation_point, interpolation_point)
    img1 = ax1.imshow(FE_grid,
                      extent=(-1, 1, -1, 1),
                      cmap='jet',
                      origin='lower')
    # ax1.set_xlabel(r'$x_1$', fontsize=14)
    # ax1.set_ylabel(r'$x_2$', fontsize=14)
    # ax1.set_title(r'$u^{\\mathcal{N}}(\\boldsymbol{\\mu})$',
    #                 fontsize=14,
    #                 y=-0.2)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar1 = fig.colorbar(img1, cax=cax)
    cbar1.formatter.set_powerlimits((0, 0))