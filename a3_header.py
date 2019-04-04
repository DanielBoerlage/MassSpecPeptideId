
import math
from collections import namedtuple
from bisect import bisect_left, bisect_right

# protein database = protein -> tag
# peptide tag database = peptide -> {tag}
# peptide mass database = sorted([(mass, peptide)])

Spectrum = namedtuple('Spectrum', ['index', 'mz', 'z', 'n_scans', 'retention_time', 'spec_data'])

proton_mass = 1.0073
h2o_mass = 18.0105
aa_mass = { 'A':  71.03711
          , 'R': 156.10111
	    , 'N': 114.04293
	    , 'D': 115.02694
	    , 'C': 160.03065
	    , 'E': 129.04259
	    , 'Q': 128.05858
	    , 'G':  57.02146
	    , 'H': 137.05891
	    , 'I': 113.08406
	    , 'L': 113.08406
	    , 'K': 128.09496
	    , 'M': 131.04049
	    , 'F': 147.06841
	    , 'P':  97.05276
	    , 'S':  87.03203
	    , 'T': 101.04768
	    , 'W': 186.07931
	    , 'Y': 163.06333
	    , 'V':  99.06841
	    }

aa = list(aa_mass.keys())

def points_near(x, points, tolerance):
	lower_bound = x - tolerance
	upper_bound = x + tolerance
	i = bisect_left (points, (lower_bound, ))
	j = bisect_right(points, (upper_bound, ))
	return points[i:j]

def build_normdist_table(max_d, n):
	# these dependencies are not required unless building the table
	import scipy.stats as st
	import numpy as np
	return [round(st.norm.pdf(d) / .398942280, 5) for d in np.linspace(0, max_d, n)]

# build_normdist_table(3.5, 100) generated this table:
normdist_table = [1.0, 0.99938, 0.9975, 0.99439, 0.99005, 0.9845, 0.97775, 0.96984, 0.96079,
0.95064, 0.93942, 0.92717, 0.91394, 0.89977, 0.88472, 0.86883, 0.85216, 0.83476, 0.8167, 0.79804,
0.77882, 0.75912, 0.73899, 0.7185, 0.6977, 0.67666, 0.65543, 0.63408, 0.61266, 0.59122, 0.56982,
0.5485, 0.52733, 0.50634, 0.48557, 0.46508, 0.44489, 0.42505, 0.40559, 0.38654, 0.36792, 0.34976,
0.33208, 0.3149, 0.29823, 0.2821, 0.2665, 0.25146, 0.23696, 0.22302, 0.20964, 0.19682, 0.18455,
0.17283, 0.16165, 0.15101, 0.14089, 0.13128, 0.12218, 0.11356, 0.10542, 0.09775, 0.09051, 0.08371,
0.07732, 0.07134, 0.06573, 0.06049, 0.05559, 0.05103, 0.04679, 0.04284, 0.03918, 0.03578, 0.03264,
0.02974, 0.02706, 0.02459, 0.02232, 0.02024, 0.01832, 0.01657, 0.01496, 0.0135, 0.01216, 0.01094,
0.00983, 0.00883, 0.00791, 0.00708, 0.00633, 0.00566, 0.00504, 0.00449, 0.004, 0.00355, 0.00315,
0.00279, 0.00247, 0.00219]

# same n as above
def normdist_func(tol, n = 100):
	if tol == 0:
		return lambda x, y: 0
	ifactor = n / tol
	def nd(x, y):
		index = int(abs(x - y) * ifactor)
		if index >= len(normdist_table):
			return 0
		return normdist_table[index]
	return nd

def lognormalize_points(points):
	# Assumes all important values of y roughly >= 1/e
	# This seems to be the case with the spectra provided
	return [(x, max(0, math.log(y) + 1)) for (x, y) in points]

def estimate_density_func(points, tol):
	nd = normdist_func(tol)
	def density(x):
		return sum([point[1] * nd(point[0], x) for point in points_near(x, points, tol)])
	return density

