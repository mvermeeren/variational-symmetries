#!/usr/bin/env sage

sys.path.append('tools')

def range(t):
	"""Returns the time variables with index at most t"""
	if even:
		return [1..t]
	else:
		if double:
			if t <= numtimes:
				return [1,3..t] 
			else:
				return [1,3..numtimes] + [numtimes+1,numtimes+3..t]
		else:
			return [1,3..t]

if double:
	dim = 2*numtimes
	dt = 2
else:
	dim = numtimes
	dt = 1


def weight(multiindex):
	"""Returns a list of weights corresponding to each time.
Usually the weights are the orders of the PDEs for the respective times."""
	return sum([multiindex[i-1]*weights[i-1] for i in [1..len(multiindex)]])

def indices_sub(wgt, length):
	"""Returns a list of all multi-indices of a given weight and length"""
	if length == 1:
		yield [int(wgt)]
	else:
		if length in range(length):
			for i in [0..int(wgt/weights[length-1])]:
				for index in indices_sub(wgt-i*weights[length-1],length-1):
					yield index+[int(i)]
		else:
			for index in indices_sub(wgt,length-1):
				yield index+[0]
				

def indices(wgt, length):
	"""Returns a list of all multi-indices of a given length and at most a given weight"""
	for w in [0..wgt]:
		for index in indices_sub(w,length):
			yield index

# Initialize t and v_ij.. variables
tstring = ""
t = []
for i in [1..dim]:
	var('t' + str(i))
	t += [eval('t' + str(i))]
	tstring += 't' + str(i) + ','
	
function('u',nargs=dim+1)
def ueval(i):
	return eval("u(" + tstring + str(i) + ")")
	
for index in indices(order,dim):
	for j in [1..components]:
		if components == 1:
			string = 'v_'
		else:
			string = 'v' + str(j) + '_'
		try:
			varnames
		except NameError:
			latexstring = string + '{'
		else:
			latexstring = varnames[j-1] + '_{'
			
		for i in [1..numtimes]:
			string += str(i)*index[i-1]
			latexstring += str(i)*index[i-1]
		if double:
			string += '_'
			for i in [1..numtimes]:
				string += str(i)*index[numtimes+i-1]
				latexstring += ('\\overline{' + str(i) + '}')*index[numtimes+i-1]
		latexstring += '}'
		var(string, latex_name=latexstring)
	
	
def dfield(index1,index2, i=1):
	"""Returns the v-variable corresponding to a given pair of multiindeces 
(and component) for a doubly infinite hierarchy"""
	if components == 1:
		string = 'v_'
	else:
		string = 'v' + str(i) + '_'
	for i in [1..len(index1)]:
		string += str(i)*index1[i-1]
	string += '_'
	for i in [1..len(index2)]:
		string += str(i)*index2[i-1]
	return eval(string)

def field(index, i=1):
	"""Returns the v-variable corresponding to a given multiindex (and component)"""
	l = len(index)
	if double:
		return dfield(index[0:l/2],index[l/2:l],i)
	else:
		if components == 1:
			string = 'v_'
		else:
			string = 'v' + str(i) + '_'
		for i in [1..len(index)]:
			string += str(i)*index[i-1]
		return eval(string)



def fieldtoindex(func):
	"""Returns multiindex corresponding to given v-variable"""
	for index in indices(order,dim):
		for i in [1..components]:
			if field(index,i) == func:
				return index
				

def deri(index):
	"""Returns the array [t_1,...,t_1,t_2...t_2,...] corresponding to a multiindex.
Can be used as input for differentiation functions."""
	deri = []
	for i in [1..len(index)]:
		deri += index[i-1]*[eval('t' + str(i))]
	return deri
    
    
def dertovar(func):
	"""Converts a function of partial derivatives of u into the corresponding function of v-variables."""
	for index in indices(order,dim):
		for i in [1..components]:
			func = func.subs(ueval(i).diff(deri(index))==field(index,i))
		if not('u' in str(func)): # stop if no more u variables present
			return func
	return func
	

def vartoder(func):
	"""Converts a function of v-variables into the corresponding function of partial derivatives of u"""
	for index in indices(order,dim):
		for i in [1..components]:
			func = func.subs(field(index,i)==ueval(i).diff(deri(index)))
		if not('_' in str(func)): # stop if no more v_ variables present
			return func
	return func
    

def vdiff(func,*var):
	"""Returns time-derivatives of a function of v-variables with respect to times *var."""
	func = vartoder(func)
	func = diff(func,*var)
	func = dertovar(func)
	return func
    
def cleantrig(func):
	"""Simplify using trigonometric identities """
	new = expand(func)
	old = 0
	for component in [1..components]:
		v = field([],component)
		iterations = 0
		while iterations < numtimes:
			if hash(new) == hash(old):
				return old
			else:
				old = new
				for k in reversed([2..2*numtimes]):
					new = new.subs( cos(v)^k == (1 - sin(v)^2)*cos(v)^(k-2) )
		return func
	else:
		return new


def summands(func):
	"""Returns terms of a sum."""
	if func == 0:
		return []
	elif type(func) == type(1+a) and func.operator() == (1+a).operator():
		return func.operands()
	else:
		return [func]
		
### integrate away terms with time derivatives as far as possible, and
### 
def vintegrate_part(func):
	"""Returns a potential boundary term for integration by pars (e.g. v_111*v_1 yields v_11*v_1)"""
	out = 0
	for j in (range(numtimes)[1:] + ['']):
		if j == '':
			addorder = order
		else:
			addorder = order-weights[j-1]
		for i in reversed([1..addorder]): 
			for component in [1..components]:
				for t in summands(func):
					v = field([],component)
					vstring = str(v)
					new = 0
					dhighest = diff(t, eval(vstring+i*'1'+str(j)))
					intcandidate = dhighest * eval(vstring+(i-1)*'1'+str(j))
					#check that it is highest derivative
					if (i+1)*'1'+str(j) in str(t) or bool(intcandidate == 0):
						new = 0
					#catch logaritmic derivatives
					elif not( vstring in str(intcandidate) ):
						new = intcandidate*log(eval(vstring+(i-1)*'1'+str(j)))
					elif components == 1 and not( vstring in str( expand(intcandidate*sin(v)/cos(v)/v) ) ):
						new = intcandidate*sin(v)/cos(v)/v * log(sin(eval(vstring+(i-1)*'1'+str(j))))
					elif components == 1 and 'v_*wzeta(2*v' in str(intcandidate):
						new = intcandidate/(v_*wzeta(2*v)) * (1/2)*log(wsigma(2*v)) #NOT ROBUST
					#highest derivative has to occur linearly
					elif not(diff(dhighest, eval(vstring+i*'1'+str(j))) == 0):
						new = 0
					else:
						newterms = summands(intcandidate)
						for term in newterms:
							if not(term == 1):
								coeff = diff(term, eval(vstring+(i-1)*'1'+str(j))) / term*eval(vstring+(i-1)*'1'+str(j))
								if not(coeff == 0):
									new += term/coeff
					#update function, add term to integral
					func = expand(func - vdiff(new,t1))
					out += new
	return out
	
def vintegrate(func):
	"""Returns the t1-integral of func. If unsuccesfull, returns 0 and sets a global warning flag."""
	out = vintegrate_part(func)
	if expand(vdiff(out,t1) - func) == 0:
		return out
	else:
		global warning
		warning = True
		textadd("Integration failed!")
		return 0
	

var('target')
def varder(dir,func,var,component=1):
	"""Returns the variational derivative of func, in the directions listed in 
"dir" (list of integers, usually of length 2)
with respect to "var", which can be given either as a field (e.g. v_12) or as a multi-index (e.g. [1,1,0,0,0])."""
	result = 0
	if type(var) == type(t1):
		var = fieldtoindex(var)
	freeorder = order - weight(var)
	for index in indices(freeorder,dim):
		if sum(index) == sum(index[d-1] for d in list(set(dir))):
			realindex = [index[i-1]+var[i-1] for i in [1..dim]]
			term = diff(func.subs(field(realindex,component)==target),target)
			term = vdiff(term.subs(target==field(realindex,component)),deri(index))
			result += (-1)^sum(index) * term
	return result


@parallel
def pvarder(i,j,func,index,component=1):
	"""Parallelizable version of varder, where dir=[i,j] and var=index."""
	if func == 0:
		return 0
	else:
	    return varder([i,j],func,index,component).combine()


def varders(matrix,index,edge=False,component=1, lagnumtimes=numtimes):
	"""Take variational derivatives of all entries of a matrix."""
	eindex = [copy(index) for k in [1..lagnumtimes]]
	if double:
		for i in range(lagnumtimes):
			eindex[i-1][numtimes+i-1] = eindex[i-1][numtimes+i-1]+1
		if edge:
			vgen = pvarder(flatten([[(numtimes+i,j,matrix[i-1,j-1],eindex[i-1],component) for i in range(lagnumtimes)] for j in range(lagnumtimes)],max_level=1))
		else:
			vgen = pvarder(flatten([[(numtimes+i,j,matrix[i-1,j-1],index,component) for i in range(lagnumtimes)] for j in range(lagnumtimes)],max_level=1))
	else:
		for i in range(lagnumtimes):
			eindex[i-1][i-1] = eindex[i-1][i-1]+1
		if edge:
			vgen = pvarder(flatten([[(i,j,matrix[i-1,j-1],eindex[i-1],component) for i in range(lagnumtimes)] for j in range(lagnumtimes)],max_level=1))
		else:
			vgen = pvarder(flatten([[(i,j,matrix[i-1,j-1],index,component) for i in range(lagnumtimes)] for j in range(lagnumtimes)],max_level=1))

	out = 0*copy(matrix)
	for i in vgen:
		if double:
			out[i[0][0][0]-numtimes-1,i[0][0][1]-1] = i[1]
		else:
			out[i[0][0][0]-1,i[0][0][1]-1] = i[1]
	return out
