import itertools
import matplotlib.pyplot as plt
import numpy as np
import pprint
import sympy
import sympy.plotting
import sympy.printing.cxxcode

plt.close("all")

(x, y, z) = sympy.symbols("x, y, z", real=True)
ex = sympy.Matrix([1, 0, 0])
ey = sympy.Matrix([0, 1, 0])
ez = sympy.Matrix([0, 0, 1])

def dot(u, v):
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

def cross(u, v):
	return sympy.Matrix([
		u[1] * v[2] - u[2] * v[1],
		u[2] * v[0] - u[0] * v[2],
		u[0] * v[1] - u[1] * v[0],
	])

def grad(v):
	return sympy.Matrix([sympy.diff(v, x), sympy.diff(v, y), sympy.diff(v, z)])

def div(v):
	return sympy.diff(v[0], x) + sympy.diff(v[1], y) + sympy.diff(v[2], z)

def curl(v):
	return sympy.Matrix([
		sympy.diff(v[2], y) - sympy.diff(v[1], z),
		sympy.diff(v[0], z) - sympy.diff(v[2], x),
		sympy.diff(v[1], x) - sympy.diff(v[0], y),
	])

def tensor(xx, yy, zz):
	return sympy.Matrix([
		[xx,  0,  0],
		[ 0, yy,  0],
		[ 0,  0, zz],
	])

def inner_product(u, v):
	return (v.transpose() * u)[0] if v.is_Matrix else v * u
	#return (v.adjoint() * u)[0] if v.is_Matrix else sympy.conjugate(v) * u

def inner_product_slice_div(u, v):
	return -inner_product(u, grad(v)) + sympy.diff(inner_product(dot(u, ez), v), z)

def inner_product_slice_curl(u, v):
	return inner_product(u, curl(v)) - sympy.diff(inner_product(cross(u, ez), v), z)

(delta_x, delta_y) = sympy.symbols("delta_x, delta_y", real=True, positive=True)

# material parameters
(vacuum_permittivity, vacuum_permeability) = sympy.symbols("vacuum_permittivity, vacuum_permeability", real=True)
(permittivity_x, permittivity_y, permittivity_z) = sympy.symbols("permittivity_x, permittivity_y, permittivity_z", complex=True)
(permeability_x, permeability_y, permeability_z) = sympy.symbols("permeability_x, permeability_y, permeability_z", complex=True)
(conductivity_x, conductivity_y, conductivity_z) = sympy.symbols("conductivity_x, conductivity_y, conductivity_z", complex=True)
impedance = sympy.symbols("impedance", complex=True)

# material tensors
permittivity_tensor = tensor(permittivity_x, permittivity_y, permittivity_z)
permeability_tensor = tensor(1 / permeability_x, 1 / permeability_y, 1 / permeability_z)

omega = sympy.symbols("omega", real=True)
gamma = sympy.symbols("gamma", complex=True)

typical_freq = 1e9
typical_values = {
	delta_x: 1e-4,
	delta_y: 1e-4,
	vacuum_permittivity: 8.854e-12,
	vacuum_permeability: 1.256e-6,
	permittivity_x: 8.854e-12,
	permittivity_y: 8.854e-12,
	permittivity_z: 8.854e-12,
	permeability_x: 1.256e-6,
	permeability_y: 1.256e-6,
	permeability_z: 1.256e-6,
	conductivity_x: 5.96e7,
	conductivity_y: 5.96e7,
	conductivity_z: 5.96e7,
	impedance: np.sqrt(1j * 2 * np.pi * typical_freq * 1.256e-6 / 5.96e7),
	omega: 2 * np.pi * typical_freq,
	gamma: 1j * 2 * np.pi * typical_freq / 3e8,
}

# basis functions
basis_xline_nodes = [
	(1 - x / delta_x),
	(    x / delta_x),
]
basis_yline_nodes = [
	(1 - y / delta_y),
	(    y / delta_y),
]
basis_rect_nodes = [
	(1 - x / delta_x) * (1 - y / delta_y),
	(    x / delta_x) * (1 - y / delta_y),
	(1 - x / delta_x) * (    y / delta_y),
	(    x / delta_x) * (    y / delta_y),
#	(1 - x / delta_x) * (    x / delta_x) * (1 - y / delta_y),
#	(1 - x / delta_x) * (    x / delta_x) * (    y / delta_y),
#	(1 - x / delta_x) * (1 - y / delta_y) * (    y / delta_y),
#	(    x / delta_x) * (1 - y / delta_y) * (    y / delta_y),
]
basis_rect_xedges = [
	(1 - y / delta_y),
	(    y / delta_y),
#	(x / delta_x - 1 / 2) * (1 - y / delta_y),
#	(x / delta_x - 1 / 2) * (    y / delta_y),
]
basis_rect_yedges = [
	(1 - x / delta_x),
	(    x / delta_x),
#	(1 - x / delta_x) * (y / delta_y - 1 / 2),
#	(    x / delta_x) * (y / delta_y - 1 / 2),
]

def generate_fem_matrix(basis_functions, bilinear_form):
	
	n = len(basis_functions)
	mat = np.empty((n, n), dtype=object)
	
	for i in range(n):
		for j in range(n):
			print("Generating FEM matrix ... %d / %d" % (i * n + j, n**2), end="\r")
			mat[i, j] = bilinear_form(basis_functions[j], basis_functions[i])
	print("Generating FEM matrix ... %d / %d" % (n**2, n**2))
	
	#pprint.pprint(mat)
	
	return mat

def check_symmetry(mat, collect_var=None, collect_num=1, scale_factors=None):
	
	n = len(mat)
	if scale_factors is None:
		scale_factors = np.ones(n, dtype=object)
	
	symm = np.zeros((n, n), dtype=int)
	for i in range(n):
		for j in range(n):
			aa = mat[i, j]
			bb = mat[j, i] #.conjugate()
			#print(aa)
			#print(bb)
			symm[i, j] = 1 if aa == bb else -1 if aa == -bb else 0
	
	print("-- symmetry %d / %d" % (abs(symm).sum(), n**2))
	
	if True:
		
		pprint.pprint(symm)
		
		def color(x):
			if x == 0:
				return 0
			return 1
		
		matc = mat if collect_var is None else np.array([[mat[i, j].collect(collect_var) for j in range(n)] for i in range(n)], dtype=object)
		
		plt.figure(figsize=(3 * (collect_num + 1), 8))
		plt.subplot(2, collect_num + 1, 1)
		plt.imshow(symm, vmin=-1, vmax=1, extent=(0, n, n, 0))
		plt.grid()
		plt.xticks(np.arange(n + 1))
		plt.yticks(np.arange(n + 1))
		plt.gca().tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
		plt.colorbar(orientation="horizontal", fraction=0.05, ticks=np.arange(-1, 2))
		for k in range(collect_num):
			plt.subplot(2, collect_num + 1, k + 2)
			matcc = matc if collect_var is None else np.array([[matc[i, j].coeff(collect_var, k) for j in range(n)] for i in range(n)], dtype=object)
			col = np.array([[color(matcc[i][j]) for j in range(n)] for i in range(n)])
			plt.imshow(col, vmin=0, vmax=1, cmap=plt.cm.rainbow, extent=(0, n, n, 0))
			plt.grid()
			plt.xticks(np.arange(n + 1))
			plt.yticks(np.arange(n + 1))
			plt.gca().tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
			plt.colorbar(orientation="horizontal", fraction=0.05, ticks=np.arange(2))
		for k in range(collect_num + 1):
			plt.subplot(2, collect_num + 1, collect_num + k + 2)
			matcc = matc if collect_var is None or k == 0 else np.array([[matc[i, j].coeff(collect_var, k - 1) * collect_var**(k-1) for j in range(n)] for i in range(n)], dtype=object)
			val = np.array([[np.nan if matcc[i][j] == 0 else complex((matcc[i][j] * scale_factors[i] * scale_factors[j]).evalf(subs=typical_values)) for j in range(n)] for i in range(n)], dtype=complex)
			col = np.log10(np.abs(val))
			#col = np.degrees(np.angle(val))
			if k == 0:
				vmin = col[np.isfinite(col)].min()
				vmax = col[np.isfinite(col)].max()
				#vmin = -180
				#vmax = 180
			plt.imshow(col, vmin=vmin, vmax=vmax, cmap=plt.cm.rainbow, extent=(0, n, n, 0))
			plt.grid()
			plt.xticks(np.arange(n + 1))
			plt.yticks(np.arange(n + 1))
			plt.gca().tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
			plt.colorbar(orientation="horizontal", fraction=0.05)
		plt.tight_layout()
		plt.show()

def generate_code(funcname, mat, collect_var=None, collect_num=1):
	
	n = len(mat)
	
	def matrix_print(target, mat):
		for i in range(n):
			for j in range(n):
				if mat[i, j] != 0:
					print("\tFemMatrix_Insert(%s, vars[%d], vars[%d], %s);" % (target, i, j, sympy.printing.cxxcode(mat[i, j])))
	
	print("void %s() {" % (funcname))
	print("\tconstexpr complex_t I(0.0, 1.0);")
	if collect_var is None:
		matrix_print("matrix", mat)
	else:
		matc = np.array([[mat[i, j].collect(collect_var) for j in range(n)] for i in range(n)], dtype=object)
		for k in range(collect_num):
			matcc = np.array([[matc[i, j].coeff(collect_var, k) for j in range(n)] for i in range(n)], dtype=object)
			matrix_print("matrix[%d]" % (k), matcc)
	print("}")

def calc_static_epot_rect():
	
	basis_functions = basis_rect_nodes
	
	def bilinear_form(u, v):
		u_e = u
		v_e = v
		res = inner_product(permittivity_tensor * grad(u_e), grad(v_e))
		return sympy.integrate(res, (x, 0, delta_x), (y, 0, delta_y)).expand()
	
	print("---- static_epot_rect ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat)
	#generate_code("FemMatrix_StaticEPot_Rect", mat)
	
	return mat

def calc_static_mpot_rect():
	
	basis_functions = basis_rect_nodes
	
	def bilinear_form(u, v):
		u_m = sympy.Matrix([0, 0, u])
		v_m = sympy.Matrix([0, 0, v])
		res = inner_product(permeability_tensor * curl(u_m), curl(v_m))
		return sympy.integrate(res, (x, 0, delta_x), (y, 0, delta_y)).expand()
	
	print("---- static_mpot_rect ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat)
	#generate_code("FemMatrix_StaticMPot_Rect", mat)
	
	return mat

def calc_empot_rect():
	
	basis_functions = (
		[(b, 0, 0, 0) for b in basis_rect_nodes] +
		[(0, b, 0, 0) for b in basis_rect_xedges] +
		[(0, 0, b, 0) for b in basis_rect_yedges] +
		[(0, 0, 0, b) for b in basis_rect_nodes]
	)
	scale_factors = (
		4 * [1 / sympy.sqrt(vacuum_permittivity)] +
		2 * [sympy.sqrt(vacuum_permeability)] +
		2 * [sympy.sqrt(vacuum_permeability)] +
		4 * [1 / sympy.sqrt(vacuum_permittivity)]
	)
	
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_e = sympy.I * u[0] * sympy.exp(-gamma * z)
		v_e = sympy.I * v[0] * sympy.exp( gamma * z)
		u_m = sympy.Matrix([u[1], u[2],  gamma / omega * u[3]]) * sympy.exp(-gamma * z)
		v_m = sympy.Matrix([v[1], v[2], -gamma / omega * v[3]]) * sympy.exp( gamma * z)
		res_e = inner_product_slice_div(-permittivity_tensor * (grad(u_e) + s * u_m), v_e)
		res_m = inner_product_slice_curl(permeability_tensor * curl(u_m), v_m) + inner_product(permittivity_tensor * (s * grad(u_e) + s**2 * u_m), v_m)
		res = (res_e + res_m).subs(z, 0)
		return sympy.integrate(res, (x, 0, delta_x), (y, 0, delta_y)).expand()
	
	print("---- empot_rect ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat, collect_var=gamma, collect_num=3, scale_factors=scale_factors)
	
	generate_code("FemMatrix_EMPot_Rect", mat, collect_var=gamma, collect_num=3)
	
	return mat

def calc_empot_xline():
	
	basis_functions = (
		[(b, 0, 0, 0) for b in basis_xline_nodes] +
		[(0, 1, 0, 0)] +
		[(0, 0, 0, b) for b in basis_xline_nodes]
	)
	scale_factors = (
		2 * [1 / sympy.sqrt(vacuum_permittivity)] +
		1 * [sympy.sqrt(vacuum_permeability)] +
		2 * [1 / sympy.sqrt(vacuum_permittivity)]
	)
	
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_e = sympy.I * u[0] * sympy.exp(-gamma * z)
		v_e = sympy.I * v[0] * sympy.exp( gamma * z)
		u_m = sympy.Matrix([u[1], u[2],  gamma / omega * u[3]]) * sympy.exp(-gamma * z)
		v_m = sympy.Matrix([v[1], v[2], -gamma / omega * v[3]]) * sympy.exp( gamma * z)
		impedance_tensor = tensor(1 / impedance, 0, 1 / impedance)
		res_e = inner_product_slice_div(-impedance_tensor * (grad(u_e) / s + u_m), v_e)
		res_m = inner_product(impedance_tensor * (grad(u_e) + s * u_m), v_m)
		res = (res_e + res_m).subs(z, 0)
		return sympy.integrate(res, (x, 0, delta_x)).expand()
	
	print("---- empot_xline ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat, collect_var=gamma, collect_num=3, scale_factors=scale_factors)
	
	generate_code("FemMatrix_EMPot_XLine", mat, collect_var=gamma, collect_num=3)
	
	return mat

def calc_empot_yline():
	
	basis_functions = (
		[(b, 0, 0, 0) for b in basis_yline_nodes] +
		[(0, 0, 1, 0)] +
		[(0, 0, 0, b) for b in basis_yline_nodes]
	)
	scale_factors = (
		2 * [1 / sympy.sqrt(vacuum_permittivity)] +
		1 * [sympy.sqrt(vacuum_permeability)] +
		2 * [1 / sympy.sqrt(vacuum_permittivity)]
	)
	
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_e = sympy.I * u[0] * sympy.exp(-gamma * z)
		v_e = sympy.I * v[0] * sympy.exp( gamma * z)
		u_m = sympy.Matrix([u[1], u[2],  gamma / omega * u[3]]) * sympy.exp(-gamma * z)
		v_m = sympy.Matrix([v[1], v[2], -gamma / omega * v[3]]) * sympy.exp( gamma * z)
		impedance_tensor = tensor(0, 1 / impedance, 1 / impedance)
		res_e = inner_product_slice_div(-impedance_tensor * (grad(u_e) / s + u_m), v_e)
		res_m = inner_product(impedance_tensor * (grad(u_e) + s * u_m), v_m)
		res = (res_e + res_m).subs(z, 0)
		return sympy.integrate(res, (y, 0, delta_y)).expand()
	
	print("---- empot_yline ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat, collect_var=gamma, collect_num=3, scale_factors=scale_factors)
	
	generate_code("FemMatrix_EMPot_YLine", mat, collect_var=gamma, collect_num=3)
	
	return mat

def calc_empot_rect_gauge():
	
	basis_functions = (
		[(b, 0, 0) for b in basis_rect_nodes] +
		[(0, b, 0) for b in basis_rect_xedges] +
		[(0, 0, b) for b in basis_rect_yedges]
	)
	scale_factors = (
		4 * [1] +
		2 * [1] +
		2 * [1]
	)
	
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_g = u[0] * sympy.exp(-gamma * z)
		v_g = v[0] * sympy.exp( gamma * z)
		u_m = sympy.Matrix([u[1], u[2], 0]) * sympy.exp(-gamma * z)
		v_m = sympy.Matrix([v[1], v[2], 0]) * sympy.exp( gamma * z)
		res_g = inner_product_slice_div(tensor(1, 1, 0) * (grad(u_g) + u_m), v_g)
		res = (res_g).subs(z, 0)
		return sympy.integrate(res, (x, 0, delta_x), (y, 0, delta_y)).expand()
	
	print("---- empot_rect_gauge ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat, collect_var=gamma, collect_num=3, scale_factors=scale_factors)
	
	generate_code("FemMatrix_EMPot_Rect_Gauge", mat, collect_var=gamma, collect_num=3)
	
	return mat

#mat_static_epot_rect = calc_static_epot_rect()
#mat_static_mpot_rect = calc_static_mpot_rect()
#mat_empot_rect = calc_empot_rect()
#mat_empot_xline = calc_empot_xline()
#mat_empot_yline = calc_empot_yline()
mat_empot_rect_gauge = calc_empot_rect_gauge()
