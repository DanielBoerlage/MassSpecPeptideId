
import re, math, random
from collections import defaultdict
from itertools import chain

import matplotlib.pyplot as plt

from a3_header import Spectrum, aa_mass, proton_mass, h2o_mass, points_near, estimate_density_func, normdist_func, lognormalize_points
from a3_io import read_protein_database, read_mgf_file


# 1.	For each spectrum, you calculate the mass of the peptide by mz×z − 1.0073×z.
#	Here mz is the mass to charge ratio of the precursor and z is the charge state.
#	The –1.0073×z in the formula is to subtract the mass of the extra protons due to the charge.

def spectrum_peptide_mass(spec):
	return float(spec.mz) * spec.z - proton_mass * spec.z


# 2.	Each protein can be digested with trypsin rule: after R or K, and not before P.

trypsin = re.compile(r'(?<=[RK])(?# snip )(?=[^P])')
def digest_protein(protein, enzyme = trypsin):
	return enzyme.split(protein)

def build_peptide_tag_database(protein_db):
	peptide_db = defaultdict(set)
	for (protein, tag) in protein_db.items():
		for peptide in digest_protein(protein):
			peptide_db[peptide].add(tag)
	return peptide_db


# 3.	For each tryptic peptide, the peptide mass is calculated as the total residue mass +18.0105.
#	The +18.0105 is because of the extra water on the peptide.

def peptide_mass(peptide):
	return sum([aa_mass[aa] for aa in peptide]) + h2o_mass

def build_peptide_mass_database(peptide_tag_db):
	return sorted([(peptide_mass(peptide), peptide) for peptide in peptide_tag_db])


# 4.	If a peptide mass matches the spectrum’s peptide mass with in error ±0.1Da,
#	then evaluate the peptide-spectrum match with your scoring function.

def spectrum_matching_peptides(spec, peptide_mass_db, tolerance = .1):
	spec_mass = spectrum_peptide_mass(spec)
	return points_near(spec_mass, peptide_mass_db, tolerance)


# 5.	After all proteins are evaluated, output the peptide with the highest matching score to the spectrum.

def best_matching_peptide(spec, peptide_mass_db):
	matches = spectrum_matching_peptides(spec, peptide_mass_db)
	return max([(score_1(peptide, spec), peptide) for (mass, peptide) in matches], default = (0, None))



# The first scoring function
# ---


# 1.	Calculate the y-ion mass of a peptide as the total residue mass of the suffix, plus 19.0178.

def y_ion_mz(residue, charge = 1):
	return (peptide_mass(residue) + charge * proton_mass) / charge

def y_ions(peptide, charge = 1):
	residues = [peptide[i:] for i in range(1, len(peptide))]
	return [(residue, y_ion_mz(residue, charge)) for residue in residues]


# 2.	Find the tallest peak within ±0.5Da of the calculated y-ion mass.
#	Suppose it’s relative intensity (defined as the ratio between the peak’s intensity and the
#	intensity of the tallest peak in the spectrum) is x, add max{0, log_10(100⋅x)} to the total score.

def tallest_peak_in_spectrum(spec_data):
	return max(spec_data, key = lambda peak: peak[1], default = (None, 0))

def tallest_peak_near_mz(mz, spec_data, tolerance = .5):
	near_mz = points_near(mz, spec_data, tolerance)
	return tallest_peak_in_spectrum(near_mz)

def score_1(peptide, spec, significant_relative_intensity = .01):
	(tallest_peak, tallest_peak_intensity) = tallest_peak_in_spectrum(spec.spec_data)
	assert (tallest_peak_intensity > 0)

	score = 0
	for (y_ion, y_ion_mz) in y_ions(peptide):
		(matching_peak, matching_peak_intensity) = tallest_peak_near_mz(y_ion_mz, spec.spec_data)
		relative_intensity = matching_peak_intensity / tallest_peak_intensity
		if relative_intensity > significant_relative_intensity:
			score = score + math.log10(100 * relative_intensity)

	return score



# The second scoring function
# ---


# The default tolerance .05 was found experimentally
def precursor_score(peptide, spec, tolerance = .05):
	nd = normdist_func(tolerance)
	return nd(peptide_mass(peptide), spectrum_peptide_mass(spec))

def b_ion_mz(residue, charge = 1):
	return (sum([aa_mass[aa] for aa in residue]) + charge * proton_mass) / charge

def b_ions(peptide, charge = 1):
	residues = [peptide[:i] for i in range(1, len(peptide))]
	return [(residue, b_ion_mz(residue, charge)) for residue in residues]

# The default tolerance 40 was chosen experimentally so the density function is smooth
def anomaly_func(logspec, tolerance = 40):
	density = estimate_density_func(logspec, tolerance)
	def specialness(mz, intensity):
		return math.sqrt(intensity / density(mz))
	return specialness

def significance_of_peak(peak_intensity, tallest_peak_intensity):
	return math.sqrt(peak_intensity / tallest_peak_intensity)

# The default tolerance .5 was found experimentally
def fragment_score(peptide, spec, tolerance = .5):
	logspec = lognormalize_points(spec.spec_data)
	(_, tallest_peak_intensity) = tallest_peak_in_spectrum(logspec)
	assert (tallest_peak_intensity > 0)

	anomaly_spec = anomaly_func(logspec)
	nd = normdist_func(tolerance)

	score = 0
	for (yb_intensity, yb_mz) in y_ions(peptide) + b_ions(peptide):
		evidence_likelihoods = []
		for (peak_mz, peak_intensity) in points_near(yb_mz, logspec, tolerance):
			significance = significance_of_peak(peak_intensity, tallest_peak_intensity)  # This value will be between 0 and 1.  Typically around .7
			specialness = anomaly_spec(peak_mz, peak_intensity)      # This value will be between 0 and 1.  Typically around .4
			locality = nd(peak_mz, yb_mz)                            # This value will be between 0 and 1.  Typically around .9
			evidence_likelihoods.append(significance * locality * specialness)
		score = score + max(evidence_likelihoods, default = 0)

	return score

def score_2(peptide, spec, alpha = 6):
	precursor = precursor_score(peptide, spec)
	fragment = fragment_score(peptide, spec)
	return fragment + alpha * precursor



# Main
# ---


def a3_output_line(spec, peptide_mass_db, peptide_tag_db, pretty_print = False):
	(score1, matching_peptide) = best_matching_peptide(spec, peptide_mass_db)
	if matching_peptide:
		matching_protein = random.choice(tuple(peptide_tag_db[matching_peptide]))
		score2 = score_2(matching_peptide, spec)
		if (pretty_print):
			return '{} {} {} {} {} {:8.5f} {:8.5f}'.format(str(spec.index).rjust(4), spec.mz.rjust(11), spec.z, matching_peptide.rjust(25), matching_protein, score1, score2)
		return '{} {} {} {} {} {:.5f} {:.5f}'.format(spec.index, spec.mz, spec.z, matching_peptide, matching_protein, score1, score2)


def main(argv):
	if (len(argv) < 3):
		import sys
		print('Usage: {} spectrum-file.mgf protein-file.fasta'.format(argv[0]), file = sys.stderr)
		return

	fasta_file = argv[2]
	protein_db = read_protein_database(fasta_file)
	peptide_tag_db = build_peptide_tag_database(protein_db)
	peptide_mass_db = build_peptide_mass_database(peptide_tag_db)

	mgf_file = argv[1]

	for spec in read_mgf_file(mgf_file):
		out_line = a3_output_line(spec, peptide_mass_db, peptide_tag_db)
		if out_line:
			print(out_line)


def plot_spectrum(spec, peptide = ''):
	fig = plt.figure(num=None, figsize=(100, 100), dpi=80, facecolor='w', edgecolor='k')
	spec_data = spec.spec_data
	(_, tallest_peak) = tallest_peak_in_spectrum(spec_data)
	for (mz, intensity) in spec_data:
		plt.axvline(mz, 0, math.log10(100 * (intensity / tallest_peak)) / 2, color = 'k', lw = 1)
	for (_, mz) in y_ions(peptide):
		plt.axvline(mz, 0, 1, color = 'r', dashes = [3, 2], lw = 1)
	for (_, mz) in b_ions(peptide):
		plt.axvline(mz, 0, 1, color = 'b', dashes = [3, 2], lw = 1)
	fig.tight_layout()
	return fig



if __name__ == "__main__":
    import sys
    main(sys.argv)
