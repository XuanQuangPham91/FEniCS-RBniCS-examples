from dolfin import *

mesh = UnitSquareMesh(2, 2)

f = Expression('sin(pi*x[0])*x[1]', degree=1)
V = FunctionSpace(mesh, 'DG', 1)

f_proj = project(f, V)

F = f_proj.vector().get_local()
X = V.tabulate_dof_coordinates()
X.resize((V.dim(), 2))

print('dof index | dof coordinate |  dof value')
for i, (x, v) in enumerate(zip(X, F)):
    print(i, x, v)
