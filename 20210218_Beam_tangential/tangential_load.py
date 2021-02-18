from dolfin import *
from mshr import *
from ufl import Jacobian, diag
# from plotting import plot
from ufl.core.multiindex import indices
from Python_module_Quang import *

import matplotlib.pyplot as plt
import time
start = time.time()


def generate_mesh_package(N):
    # N = 10  # segment for generate mesh
    R = 1
    rectangle1 = Rectangle(Point(0., 0.), Point(1., 1.))
    rectangle2 = Rectangle(Point(0., 0.), Point(1., -1.))
    rectangle3 = Rectangle(Point(0., 0.), Point(-1., -1.))
    rectangle4 = Rectangle(Point(0., 0.), Point(-1., 1.))

    circle = Circle(Point(0., 0.), R, segments=5 * N)
    domain = circle - rectangle2 - rectangle3 - rectangle4
    domain.set_subdomain(1, circle - rectangle2 - rectangle3 - rectangle4)
    # domain.set_subdomain(2, rectangle - circle)

    mesh = generate_mesh(domain, N)

    # Create subdomains
    subdomains = MeshFunction("size_t", mesh,
                              mesh.topology().dim(), mesh.domains())

    # Create boundaries
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0.)

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], 0.)

    class Curve(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near((x[0]**2 + x[1]**2), 1.)

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)
    left = Left()
    left.mark(boundaries, 1)
    bottom = Bottom()
    bottom.mark(boundaries, 2)
    curve = Curve()
    curve.mark(boundaries, 3)

    #save
    File("data/elastic_block.xml") << mesh
    File("data/elastic_block_physical_region.xml") << subdomains
    File("data/elastic_block_facet_region.xml") << boundaries
    XDMFFile("data/elastic_block.xdmf").write(mesh)
    XDMFFile("data/elastic_block_physical_region.xdmf").write(subdomains)
    XDMFFile("data/elastic_block_facet_region.xdmf").write(boundaries)


class ElasticBlock():  #, young_modulus_E, Nu):
    def __init__(self, N, mesh, subdomains, boundaries, **kwargs):
        # Define mesh, subdomains and boundaries
        self.mesh = mesh
        self.subdomains = subdomains
        self.boundaries = boundaries

        # Define material properties
        self.E = 80e3  # MPa
        self.nu = 0.2  # MPa
        self.lambda_1 = self.E * self.nu / ((1.0 + self.nu) *
                                            (1.0 - 2.0 * self.nu))
        self.lambda_2 = self.E / (2.0 * (1.0 + self.nu))
        self.C1111 = self.lambda_1 + 2 * self.lambda_2
        self.C1122 = self.lambda_1
        self.C1212 = self.lambda_2

        # self.t = self.tangent2D(mesh)

        # Call functions -------------------------------------------------------
        self.assemble_operator()
        self.solve()
        # self.von_mises()

    # dx = Measure("dx")(subdomain_data=self.subdomains)
    # ds = Measure("ds")(subdomain_data=self.boundaries)

    # def tangent2D(self):
    #     norm = FacetNormal(self.mesh)
    #     return as_vector([norm[1], -norm[0]])

    # def tangent3D(self):
    #     t = Jacobian(self.mesh)
    #     return as_vector([t[0, 0], t[1, 0], t[2, 0]]) / sqrt(inner(t, t))

    def assemble_operator(self):

        # Add mesh -------------------------------------------------------------
        dx = Measure("dx")(subdomain_data=self.subdomains)
        ds = Measure("ds")(subdomain_data=self.boundaries)

        self.V = VectorFunctionSpace(self.mesh, "Lagrange", 1)
        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)

        self.bc = [
            DirichletBC(self.V.sub(0), 0.0, self.boundaries, 1),
            DirichletBC(self.V.sub(1), 0.0, self.boundaries, 2),
            # DirichletBC(self.V, Constant((0., 0.)), self.boundaries, 1),
            # DirichletBC(self.V, Constant((0., 0.)), self.boundaries, 2),
        ]
        self.f = Constant(10.0)
        # self.f = Expression("pow(x[0],2) + pow(x[1],2)", degree=1)

        # define tangential vector
        # self.Vv = VectorFunctionSpace(self.mesh, "Lagrange", 1)
        # self.Vf = FunctionSpace(self.mesh, "Lagrange", 1)
        # self.Vt = TensorFunctionSpace(self.mesh, "Lagrange", 1)
        # self.tau = TestFunction(self.Vv)
        # self.phi = TrialFunction(self.Vf)

        self.normal = FacetNormal(self.mesh)
        # self.tangential = as_vector([self.norm[1], -self.norm[0]])

        # self.tau_n = dot(self.tau, self.norm)
        # self.tau_t = dot(self.tau, self.tangential)
        # self.P_n = dot(self.f, self.norm)
        # self.P_t = dot(self.f, self.tangential)

        # Weak form
        # i, j = indices(2)
        self.a = (
            self.C1111 * (self.u[0].dx(0) * self.v[0].dx(0) + self.u[1].dx(1) * self.v[1].dx(1)) \
            + self.C1122 * (self.u[1].dx(1) * self.v[0].dx(0) + self.u[0].dx(0) * self.v[1].dx(1)) \
            + self.C1212 * (
                self.u[0].dx(1) * self.v[0].dx(1) +
                self.u[0].dx(1) * self.v[1].dx(0) +
                self.u[1].dx(0) * self.v[0].dx(1) +
                self.u[1].dx(0) * self.v[1].dx(0))
            ) * dx
        # self.l = self.f[i] * self.v[i] * ds
        self.l = self.f * dot(self.normal, self.v) * ds
        # self.l = dot(self.P_n, self.tau_n) * ds \
        #         + dot(self.P_n, self.tau_n) * ds

    def solve(self):
        ## Call solver
        self.u = Function(self.V, name='displacement')
        solve(self.a == self.l, self.u, self.bc)
        # return self.u


def FE_tangent(N):
    mesh = Mesh("data/elastic_block.xml")
    subdomains = MeshFunction("size_t", mesh,
                              "data/elastic_block_physical_region.xml")
    boundaries = MeshFunction("size_t", mesh,
                              "data/elastic_block_facet_region.xml")
    problem = ElasticBlock(N=N,
                           mesh=mesh,
                           subdomains=subdomains,
                           boundaries=boundaries)

    FEniCS_plot_mode(u=problem.u,
                     number_of_figure=True,
                     title='u_FE',
                     mode='displacement',
                     savefig=True)
    save_HDF5(u=problem.u, mesh=problem.mesh, title='u_FE')


if __name__ == "__main__":
    N = 20
    new_geometry = True
    if new_geometry:
        generate_mesh_package(N)
        FE_tangent(N)
    else:
        FE_tangent(N)
