from sympy import Matrix, pprint, Integer
import itertools
class BasisElement(Matrix):
	"""docstring for BasisElement"""
	linear_factors = []

	def __init__(self, *args, **kwargs):
		super(BasisElement, self).__init__(*args, **kwargs)
		self.linear_factors = [0] * (self.cols - 1)
		self.linear_factors.append(1)

	def __add__(self, other):
		z = super(BasisElement, self).__add__(other)
		z.linear_factors = [x + y for x,y in itertools.izip(self.linear_factors, other.linear_factors)]
		return z


	def __sub__(self, other):
		z = super(BasisElement, self).__sub__(other)
		z.linear_factors = [x - y for x,y in itertools.izip(self.linear_factors, other.linear_factors)]
		return z


	def __mul__(self, other):
		# if multiplying by a matrix, do whatever the superclass does
		# should not be needed
		if getattr(other, 'is_Matrix', False):
			return super(BasisElement, self).__mul__(other)


		# otherwise multiply the linear factors as well
		z = super(BasisElement, self).__mul__(other)
		z.linear_factors = [other * x for x in self.linear_factors]
		return z


	def lift_single_choice(self, row_to_add):
		print "self"
		pprint(self)
		row_to_add = row_to_add[:, -1]
		print "linear factors"
		print self.linear_factors
		linear_factors_vector = Matrix([self.linear_factors])
		print "row to add"
		pprint(row_to_add)
		print "linear factos vec"
		pprint(linear_factors_vector)

		last_el = (linear_factors_vector * row_to_add)[-1]
		linear_factors = self.linear_factors
		self = self.row_join(Matrix([last_el]))
		self.linear_factors = linear_factors
		return self

	def lift_multiple_choice(self, row_to_add):
		self.linear_factors.append(0)
		# linear_factors = self.linear_factors
		# last_el = 0
		# pprint(row_to_add)
		# while True:
		# 	h = self.row_join(Matrix([last_el]))
		# 	if (row_to_add.inv()*h.T).has(Integer):
		# 		self = h
		# 		break
		# 	last_el = last_el + 1

		linear_factors_vector = Matrix([self.linear_factors])
		print "row to add"
		pprint(row_to_add)
		last_el = (linear_factors_vector * row_to_add)[-1]
		while last_el < 0:
			last_el = last_el + row_to_add[-1]
			self.linear_factors[-1] = self.linear_factors[-1] +1
		pprint(last_el)
		if last_el >= row_to_add[-1]:
			self.linear_factors[-1] = -1 * (last_el/row_to_add[-1])
			last_el = last_el % row_to_add[-1]
		
		print "last_el"
		pprint(last_el)
		#crazy mutation stuff ahppening here
		linear_factors = self.linear_factors
		self = self.row_join(Matrix([last_el]))
		self.linear_factors = linear_factors
		return self
		


		


