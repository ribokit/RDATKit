import pickle
import pdb
from numpy import *
import numpy as np
from cvxmod import *
from cvxmod.atoms import *
from cvxmod.sets import *
from rdatkit import secondary_structure
from rdatkit import settings
import scipy.stats as stats

"""
Next best structure addition (i.e. keep adding structures from the ensemble to optimize)
"""
def by_nbs_addition(sequence, data, database='default.db.dists', use_data=False, nstructs=20):
    db = pickle.load(open('%s/models/%s' % (settings.MAPPING_DATABASE_PATH, database)))
    if use_data:
	structures = secondary_structure.fold(sequence, mapping_data=data, nstructs=nstructs)
    else:
	structures = secondary_structure.fold(sequence, nstructs=nstructs)
    if nstructs > len(structures):
        print 'WARNING: Found %d non trivial structures in the ensemble, not intended %d' % (len(structures), nstructs) 
    struct_probs = []
    struct_lh = []
    for struct in structures:
	probs, lh = struct.likelihood(data, db_obj=db)  
	struct_probs.append(probs)
	struct_lh.append(-log(lh))
    currsol = 0.
    prevsol = float('-Inf')
    numstructs = 1
    while  numstructs <= len(structures):
	"""
	 Solve
	    minimize -sum( ci*log(p(Mi|D)))  over models Mi, given data D, for c
	    0 <= c <= 1 for all i

	    p(Mi|D) are given in struct_lh
	"""
	dim = numstructs
	p = cvxopt.matrix(struct_lh[:dim])
	c = variable(dim, 'c')
	C1 = (c <= 1)
	C2 = (c >= 0)
	C3 = (sum(c) == 1)
	lp = op(min(p.trans()*c), [C1,C2,C3])
	lp.solve()
	print 'Status is %s' % lp.status
	"""
	G = cvxopt.spmatrix(1., range(dim), range(dim))
	h = cvxopt.matrix([1.]*dim)
	A = cvxopt.spmatrix(-1., range(dim), range(dim))
	b = cvxopt.matrix([0.]*dim)
	sol = solvers.lp(c, G, h, A, b)
	opt = array(sol['x'])
	"""
	opt = array(c.value)
	prevsol = currsol
        print 'Solution %s' % opt
        print 'Struct lh %s' % struct_lh[:numstructs]
	currsol = sum(array([opt[i]*struct_lh[i] for i in range(numstructs)]))
	numstructs += 1
    numstructs -= 1
    print 'Finished...'
    print 'Used %d of %d structures' % (numstructs, len(structures))
    return structures[:numstructs], opt, currsol 
    

def by_objective_minimization(sequence, data, t, objfun='norm_l1', structures=[], use_data=False, nstructs=1000):
    print 'Generating structures, this can take a long time...'
    if use_data:
	structures = secondary_structure.sample(sequence, mapping_data=data, nstructs=nstructs, unique=True)
    else:
	if not structures:
	    structures = secondary_structure.sample(sequence, nstructs=nstructs, unique=True)
	else:
	    structures = structures
    if nstructs > len(structures):
        print 'WARNING: Found %d non trivial structures in the ensemble, not intended %d' % (len(structures), nstructs) 
    print 'Preparing inputs for the "Dantzig Selector"'
    A = zeros([len(sequence), len(structures)])
    Aa = np.zeros([len(sequence), len(structures)])
    # CVXMOD having trouble converting to numpy arrays?
    for j, struct in enumerate(structures):
	for i in range(len(sequence)):
	    if struct.dbn[i] == '.':
		A[i,j] = 1.
		Aa[i,j] = 1.
    """
     Solve
	minimize obj_fun( beta )  subject to
	||A'(data-A*beta)||_inf <= (1+(1/t))*sqrt(2*log(len(sequence)))*sigma

	Where sigma is the standard deviation of the data (viewed as a gamma distribution)

    """
    p = len(sequence)
    galpha, gloc, gbeta = stats.gamma.fit(data)
    sigma = sqrt(galpha*(1/gbeta))
    if t != 0:
	l = (1+(1/t))*sqrt(2*log(300))
    else:
	l = 0.
    data[isnan(data)] = 1.
    beta = optvar('b', len(structures))
    Am = param('A', value=A)
    Amt = param('At', value=transpose(A))
    y = param('y', value=matrix(data.tolist()))
    lam = param('lam', value=matrix([l*sigma]))
    con = optvar('con', 1)
    lp = problem()
    print 'Executing the LP'
    if objfun == 'norm_l1':
	lp.constr = [norminf(Amt*(y - Am*beta)) <= lam, beta >= 0 ]
	lp.objective = minimize(norm1(beta))
    if objfun == 'entropy':
	lp.constr = [norminf(Amt*(y - Am*beta)) <= lam, beta >= 0, sum(beta) == 1]
	lp.objective = maximize(entropy(beta))
    if objfun == 'lasso':
	lp.constr = [beta >= 0]
	lp.objective = minimize(norm2(y - Am*beta) + t*norm1(beta))
    print lp
    lp.solve()
    opt = array(beta.value)
    print 'Solution (coefficients) %s' % opt
    error = float(np.linalg.norm(array(data) - np.dot(Aa, opt)))
    predicted = np.dot(Aa, opt)
    print '%d structures have non-zero coefficients' % (np.sum([x > 0 for x in opt[:,0]]))
    return opt, predicted, structures, error

