from sys import path
from dolfin import *
from mshr import *
from Python_module_Quang import *
from numpy.core.records import array
from numpy.lib.function_base import append

import matplotlib.pyplot as plt
import sys
import numpy as np

sys.path.append(
    '/media/xuanquang/Gaumap Lab data/05_Git_project/FEniCS-RBniCS-examples/20210218_2D_tangential_load/'
)

try:
    import tangential_load
except:
    print("fail to import tangential_load.py")

format = "png"

# mesh = RectangleMesh(Point(0., 0.), Point(1., 1.), 2, 2, "right")
# # V =
# V = VectorFunctionSpace(mesh, "CG", 1)

# --------------------------------------------------------------
# coordinates_mesh = mesh.coordinates()
# coordinates_V_space = V.tabulate_dof_coordinates()
# # coordinates = mesh.coordinates()
# print("nodes coordinates: \n %d", coordinates_mesh)
# print("dof coordinates: \n %d", coordinates_V_space)
# # print(coordinates_V_space)

# plt.figure()
# title = "mesh"
# plot(mesh)
# # plt.savefig('solution/%s.%s' % (title, format), format='%s' % format)
# plt.savefig(f'solution/{title}.{format}', format=format)
# plt.close()

# dim = F.dim()
# N = mesh.geometry().dim()

# coor = F.dofmap().tabulate_all_coordinates(mesh).reshape(dim, N)
# fx_dofs = F.sub(0).dofmap().dofs()
# fy_dofs = F.sub(1).dofmap().dofs()


def coordinates_operator(mesh, u):
    nodal_values = u.vector().get_local()
    coordinates = mesh.coordinates()
    # print(coordinates.shape)
    x = coordinates[:, 0]
    y = coordinates[:, 1]
    print(f"len(x): {len(x)}")
    print(f"len(y): {len(y)}")

    # print(f"coordinates[1]: {coordinates[1]}")
    # print(f"u(coordinates[1]): {u(coordinates[1])}")
    # vertex_values = u.compute_vertex_values()
    # print(f"vertex_values.shape: {vertex_values.shape}")

    ux, uy = u.split(True)
    uX = ux.vector().get_local()
    uY = uy.vector().get_local()
    return x, y, uX, uY, nodal_values
    # return x, y, nodal_values


def extract_ux_uy(u):
    u_x = []
    u_y = []
    u = array(u.vector().get_local())
    print(u[0])
    for i in range(len(u)):
        # coordinates_V_space.append(u_FE.vector().get_local[i])
        if i % 2 == 0:
            u_x.append(u[i])
        else:
            u_y.append(u[i])
    print(f"len of u_x: {len(u_x)}, {u_x[0]}")
    print(f"len of u_y: {len(u_y)}, {u_y[0]}")
    return u_x, u_y


if __name__ == "__main__":
    mesh = Mesh("data/elastic_block.xml")
    # subdomains = MeshFunction("size_t", mesh,
    #                           "data/elastic_block_physical_region.xml")
    # boundaries = MeshFunction("size_t", mesh,
    #                           "data/elastic_block_facet_region.xml")

    # problem = tangential_load.ElasticBlock(N=20,
    #                                        mesh=mesh,
    #                                        subdomains=subdomains,
    #                                        boundaries=boundaries)

    # FEniCS_plot_mode(u=problem.u,
    #                  number_of_figure=True,
    #                  title='u_FE',
    #                  mode='displacement',
    #                  savefig=True)

    V = VectorFunctionSpace(mesh, "Lagrange", 1)

    u_FE = load_HDF5(V, mesh, title='u_FE')
    plot(u_FE, "displacement")

    x, y, ux, uy, nodal_values = coordinates_operator(mesh, u_FE)
    print(x)
    print(f"ux (len = {len(ux)})")
    print(f"uy (len = {len(uy)})")
    # print(f"nodal_values coordinate: {nodal_values}")
    print(f"Number of nodal_values: {len(nodal_values)}")
    print(f"Shape of nodal_values: {len(nodal_values)}")

    coordinates_mesh = mesh.coordinates()
    print(f"coordinates_mesh: {coordinates_mesh}")
    print(f"x, y: \n {x, y}")
    title = "coordinates_mesh"
    np.savetxt(f"solution/{title}.txt", np.array(coordinates_mesh))

    # print("nodes coordinates: \n %d", coordinates_mesh)
    print(f"Number of nodes: {len(coordinates_mesh)}")
    print(f"node shape: {coordinates_mesh.shape}")
    coordinates_V_space = V.tabulate_dof_coordinates()
    # print("nodes coordinates: \n %d", coordinates_mesh)
    print(f"Number of dof: {len(coordinates_V_space)}")
    print(f"dof shape: {coordinates_V_space.shape}")
    print(f"dof: {coordinates_V_space}")

    ux, uy = extract_ux_uy(u=u_FE)

    # df = pd.DataFrame(data=numpy_data,
    #                   index=["row1", "row2"],
    #                   columns=["column1", "column2"])

    # plt.show()
