
import random, statistics
from collections import Counter

import matplotlib.pyplot as plt

from a3_header import aa
from a3_io import read_protein_database, read_mgf_file
from a3 import build_peptide_tag_database, build_peptide_mass_database, best_matching_peptide, score_1, score_2, precursor_score


def fake_protein(length):
	return ''.join(random.choices(aa, k=length))

def realistic_protein_len(protein_db):
	lengths = [len(protein) for protein in protein_db]
	return lambda: random.choice(lengths)

def build_td_protein_database(protein_db):
	gen_len = realistic_protein_len(protein_db)
	td_db = protein_db.copy()
	for i in range(len(protein_db)):
		td_db[fake_protein(gen_len())] = 'decoy{0:05d}'.format(i)
	return td_db

def td_scores(score_func, specs, peptide_mass_db, peptide_tag_db):
	decoy_scores = []
	target_scores = []

	for spec in specs:
		(_, matching_peptide) = best_matching_peptide(spec, peptide_mass_db)
		if matching_peptide:
			score = score_func(matching_peptide, spec)
			matching_protein = random.choice(tuple(peptide_tag_db[matching_peptide]))
			if (matching_protein.startswith('decoy')):
				decoy_scores.append(score)
			else:
				target_scores.append(score)

	decoy_scores.sort(reverse = True)
	target_scores.sort(reverse = True)
	return (target_scores, decoy_scores)

def tda_performance(score_func, specs, peptide_mass_db, peptide_tag_db, fdr_tolerance = .05):
	(target, decoy) = td_scores(score_func, specs, peptide_mass_db, peptide_tag_db)
	def fdr():
		if not target:
			return 0
		return len(decoy) / len(target)

	while fdr() > fdr_tolerance:
		if (decoy[-1] < target[-1]):
			decoy.pop()
		else:
			target.pop()

	return len(target)

def tda_performace_sample(sample_size, score_func, specs, protein_db):
	perfs = []
	for i in range(sample_size):
		print('\r[{:4d}/{}]'.format(i, sample_size), end = '')
		td_db = build_td_protein_database(protein_db)
		peptide_tag_db = build_peptide_tag_database(td_db)
		peptide_mass_db = build_peptide_mass_database(peptide_tag_db)
		perfs.append(tda_performance(score_func, specs, peptide_mass_db, peptide_tag_db))
	print('\r', end = '')
	return perfs

def tda_performance_statistics(perfs):
	mean = statistics.mean(perfs)
	stddev = statistics.stdev(perfs)
	return (mean, stddev, 'mean: {:.2f}, stdev: {:.2f}'.format(statistics.mean(perfs), statistics.stdev(perfs)))

def tda_performance_hist(perfs, buckets = 50):
	plt.hist(perfs, buckets)
	plt.show()

def optimize_score_parameter(score_func_of_param, sample_size, param_min, param_max, resolution, specs, protein_db):
	step_size = (param_max - param_min) / resolution
	param_space = [param_min + step_size * i for i in range(resolution)]
	means = []
	stddevs = []
	for param in param_space:
		score_func = score_func_of_param(param)
		perfs = tda_performace_sample(sample_size, score_func, specs, protein_db)
		(mean, stddev, _) = tda_performance_statistics(perfs)
		means.append(mean)
		stddevs.append(stddev)
		print('param = {:.4}, mean: {:.2f}, stddev: {:.2f}'.format(param, mean, stddev))

def main(argv):
	fasta_file = 'data/ups.fasta'
	mgf_file = 'data/larger-test.mgf'
	protein_db = read_protein_database(fasta_file)
	specs = read_mgf_file(mgf_file, keep_in_memory = True)
	for score_func in [score_2]:
		perfs = tda_performace_sample(1500, score_func, specs, protein_db)
		print(tda_performance_statistics(perfs)[2])
		tda_performance_hist(perfs)


if __name__ == "__main__":
    import sys
    main(sys.argv)

