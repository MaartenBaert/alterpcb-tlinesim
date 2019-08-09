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
speed_of_light = sympy.symbols("SPEED_OF_LIGHT")
(permittivity_x, permittivity_y, permittivity_z) = sympy.symbols("permittivity_x, permittivity_y, permittivity_z", complex=True)
(permeability_x, permeability_y, permeability_z) = sympy.symbols("permeability_x, permeability_y, permeability_z", complex=True)
(conductivity_x, conductivity_y, conductivity_z) = sympy.symbols("conductivity_x, conductivity_y, conductivity_z", complex=True)
impedance = sympy.symbols("impedance", complex=True)

# material tensors
permittivity_tensor = tensor(permittivity_x, permittivity_y, permittivity_z)
permeability_tensor = tensor(1 / permeability_x, 1 / permeability_y, 1 / permeability_z)

omega = sympy.symbols("omega", real=True)
effective_index = sympy.symbols("effective_index", complex=True)

typical_freq = 1e6
typical_values = {
	delta_x: 1e-4,
	delta_y: 1e-4,
	speed_of_light: 299792458.0,
	#vacuum_permittivity: 8.854e-12,
	#vacuum_permeability: 1.256e-6,
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
	effective_index: 1.4,
}

# shorthands for basis functions
bx0 = 1 - x / delta_x
bx1 = x / delta_x
bx2 = x / delta_x - sympy.S(1) / 2
by0 = 1 - y / delta_y
by1 = y / delta_y
by2 = y / delta_y - sympy.S(1) / 2

# H(grad) phi/Az basis functions
basis_hz_xline = [
	bx0,
	bx1,
	bx0 * bx1,
]
basis_hz_yline = [
	by0,
	by1,
	by0 * by1,
]
basis_hz_rect = [
	bx0 * by0,
	bx1 * by0,
	bx0 * by1,
	bx1 * by1,
	bx0 * bx1 * by0,
	bx0 * bx1 * by1,
	bx0 * by0 * by1,
	bx1 * by0 * by1,
	#bx0 * bx1 * by0 * by1,
]

# H(curl) Ax basis functions
basis_hx_xline = [
	1,
	bx2,
]
basis_hx_rect = [
	by0,
	by1,
	#by0 * by1,
	bx2 * by0,
	bx2 * by1,
	#bx2 * by0 * by1,
]

# H(curl) Ay basis functions
basis_hy_yline = [
	1,
	by2,
]
basis_hy_rect = [
	bx0,
	bx1,
	#bx0 * bx1,
	bx0 * by2,
	bx1 * by2,
	#bx0 * bx1 * by2,
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

def check_symmetry(mat, collect_var=None, collect_num=1):
	
	n = len(mat)
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
		
		#pprint.pprint(symm)
		
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
			val = np.array([[np.nan if matcc[i][j] == 0 else complex(matcc[i][j].evalf(subs=typical_values)) for j in range(n)] for i in range(n)], dtype=complex)
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

def tocpp(f):
	code = sympy.printing.cxxcode(f)
	code = code.replace("std::pow(omega, 2)", "square(omega)")
	code = code.replace("std::pow(SPEED_OF_LIGHT, 2)", "square(SPEED_OF_LIGHT)")
	code = code.replace("1/", "1.0/")
	return code

def generate_code(funcname, mat, collect_vars=None):
	
	n = len(mat)
	
	def matrix_print(target, mat):
		for i in range(n):
			for j in range(n):
				if mat[i, j] != 0:
					code = tocpp(mat[i, j])
					print("\tFemMatrix_Insert(%s, vars[%d], vars[%d], %s);" % (target, i, j, code))
	
	print("void %s() {" % (funcname))
	print("\tconstexpr complex_t I(0.0, 1.0);")
	if collect_vars is None:
		matrix_print("matrix", mat)
	else:
		for k in range(len(collect_vars)):
			matc = np.array([[mat[i, j].collect(collect_vars[k][0]).coeff(collect_vars[k][0], collect_vars[k][1]) for j in range(n)] for i in range(n)], dtype=object)
			matrix_print("matrix[%d]" % (k), matc)
	print("}")

def calc_static_epot_rect():
	
	basis_functions = basis_hz_rect
	
	def bilinear_form(u, v):
		u_e = u
		v_e = v
		res = inner_product_slice_div(-permittivity_tensor * grad(u_e), v_e)
		return sympy.integrate(res, (x, 0, delta_x), (y, 0, delta_y)).expand()
	
	print("---- static_epot_rect ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat)
	
	generate_code("FemMatrix_StaticEPot_Rect", mat)
	
	return mat

def calc_static_mpot_rect():
	
	basis_functions = basis_hz_rect
	
	def bilinear_form(u, v):
		u_m = sympy.Matrix([0, 0, u])
		v_m = sympy.Matrix([0, 0, v])
		res = inner_product_slice_curl(permeability_tensor * curl(u_m), v_m)
		return sympy.integrate(res, (x, 0, delta_x), (y, 0, delta_y)).expand()
	
	print("---- static_mpot_rect ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat)
	
	generate_code("FemMatrix_StaticMPot_Rect", mat)
	
	return mat

def calc_static_mpot_xline():
	
	basis_functions = (
		[(b, 0) for b in basis_hz_xline] +
		[(0, 1)]
	)
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_m = sympy.Matrix([0, 0, u[0] - u[1]])
		v_m = sympy.Matrix([0, 0, v[0] - v[1]])
		impedance_tensor = tensor(1 / impedance, 0, 1 / impedance)
		res = inner_product(impedance_tensor * (s * u_m), v_m)
		return sympy.integrate(res, (x, 0, delta_x)).expand()
	
	print("---- static_mpot_xline ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat)
	
	generate_code("FemMatrix_StaticMPot_XLine", mat)
	
	return mat

def calc_static_mpot_yline():
	
	basis_functions = (
		[(b, 0) for b in basis_hz_yline] +
		[(0, 1)]
	)
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_m = sympy.Matrix([0, 0, u[0] - u[1]])
		v_m = sympy.Matrix([0, 0, v[0] - v[1]])
		impedance_tensor = tensor(0, 1 / impedance, 1 / impedance)
		res = inner_product(impedance_tensor * (s * u_m), v_m)
		return sympy.integrate(res, (y, 0, delta_y)).expand()
	
	print("---- static_mpot_yline ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat)
	
	generate_code("FemMatrix_StaticMPot_YLine", mat)
	
	return mat

def calc_empot_rect():
	
	basis_functions = (
		[(b, 0, 0, 0) for b in basis_hz_rect] +
		[(0, b, 0, 0) for b in basis_hx_rect] +
		[(0, 0, b, 0) for b in basis_hy_rect] +
		[(0, 0, 0, b) for b in basis_hz_rect]
	)
	collect_vars = [
		(effective_index, 0),
		(effective_index, 2),
	]
	
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_e = u[0] * sympy.exp(-effective_index / speed_of_light * s * z)
		v_e = v[0] * sympy.exp( effective_index / speed_of_light * s * z)
		u_m = sympy.Matrix([u[1], u[2],  effective_index * u[3]]) / speed_of_light * sympy.exp(-effective_index / speed_of_light * s * z)
		v_m = sympy.Matrix([v[1], v[2], -effective_index * v[3]]) / speed_of_light * sympy.exp( effective_index / speed_of_light * s * z)
		res_e = inner_product_slice_div(-permittivity_tensor * (grad(u_e) + s * u_m), v_e)
		res_m = inner_product_slice_curl(permeability_tensor * curl(u_m), v_m) + inner_product(permittivity_tensor * (s * grad(u_e) + s**2 * u_m), v_m)
		res = (res_e + res_m).subs(z, 0)
		return sympy.integrate(res, (x, 0, delta_x), (y, 0, delta_y)).expand()
	
	print("---- empot_rect ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat, collect_var=effective_index, collect_num=3)
	
	generate_code("FemMatrix_EMPot_Rect", mat, collect_vars=collect_vars)
	
	return mat

def calc_empot_xline():
	
	basis_functions = (
		[(b, 0, 0, 0) for b in basis_hz_xline] +
		[(0, b, 0, 0) for b in basis_hx_xline] +
		[(0, 0, 0, b) for b in basis_hz_xline]
	)
	collect_vars = [
		(effective_index, 0),
		(effective_index, 2),
	]
	
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_e = u[0] * sympy.exp(-effective_index / speed_of_light * s * z)
		v_e = v[0] * sympy.exp( effective_index / speed_of_light * s * z)
		u_m = sympy.Matrix([u[1], u[2],  effective_index * u[3]]) / speed_of_light * sympy.exp(-effective_index / speed_of_light * s * z)
		v_m = sympy.Matrix([v[1], v[2], -effective_index * v[3]]) / speed_of_light * sympy.exp( effective_index / speed_of_light * s * z)
		impedance_tensor = tensor(1 / impedance, 0, 1 / impedance)
		res_e = inner_product_slice_div(-impedance_tensor * (grad(u_e) / s + u_m), v_e)
		res_m = inner_product(impedance_tensor * (grad(u_e) + s * u_m), v_m)
		res = (res_e + res_m).subs(z, 0)
		return sympy.integrate(res, (x, 0, delta_x)).expand()
	
	print("---- empot_xline ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat, collect_var=effective_index, collect_num=3)
	
	generate_code("FemMatrix_EMPot_XLine", mat, collect_vars=collect_vars)
	
	return mat

def calc_empot_yline():
	
	basis_functions = (
		[(b, 0, 0, 0) for b in basis_hz_yline] +
		[(0, 0, b, 0) for b in basis_hy_yline] +
		[(0, 0, 0, b) for b in basis_hz_yline]
	)
	collect_vars = [
		(effective_index, 0),
		(effective_index, 2),
	]
	
	def bilinear_form(u, v):
		s = sympy.I * omega
		u_e = u[0] * sympy.exp(-effective_index / speed_of_light * s * z)
		v_e = v[0] * sympy.exp( effective_index / speed_of_light * s * z)
		u_m = sympy.Matrix([u[1], u[2],  effective_index * u[3]]) / speed_of_light * sympy.exp(-effective_index / speed_of_light * s * z)
		v_m = sympy.Matrix([v[1], v[2], -effective_index * v[3]]) / speed_of_light * sympy.exp( effective_index / speed_of_light * s * z)
		impedance_tensor = tensor(0, 1 / impedance, 1 / impedance)
		res_e = inner_product_slice_div(-impedance_tensor * (grad(u_e) / s + u_m), v_e)
		res_m = inner_product(impedance_tensor * (grad(u_e) + s * u_m), v_m)
		res = (res_e + res_m).subs(z, 0)
		return sympy.integrate(res, (y, 0, delta_y)).expand()
	
	print("---- empot_yline ----")
	mat = generate_fem_matrix(basis_functions, bilinear_form)
	check_symmetry(mat, collect_var=effective_index, collect_num=3)
	
	generate_code("FemMatrix_EMPot_YLine", mat, collect_vars=collect_vars)
	
	return mat

#mat_static_epot_rect = calc_static_epot_rect()
#mat_static_mpot_rect = calc_static_mpot_rect()
#mat_static_mpot_xline = calc_static_mpot_xline()
#mat_static_mpot_yline = calc_static_mpot_yline()
mat_empot_rect = calc_empot_rect()
mat_empot_xline = calc_empot_xline()
mat_empot_yline = calc_empot_yline()
