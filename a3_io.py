
import re, itertools

from a3_header import Spectrum

def typed_regex_parser(tr, flags = 0, custom_types = {}):
	typed_capture_search = r'\(\?P<([^:]+):([^>]*)>([^\)]*)\)'
	typed_capture_replace = r'(?P<\1>\3)'
	type_functions = {'int': int, 'float': float, 'str': str, **custom_types}
	type_function_for_sym = {match[0]: type_functions[match[1]] for match in re.findall(typed_capture_search, tr)}
	r = re.compile(re.sub(typed_capture_search, typed_capture_replace, tr), flags)
	def matches(s):
		for match in r.finditer(s):
			yield {sym: type_function_for_sym.get(sym, str)(val) for (sym, val) in match.groupdict().items() if val is not None}
	return matches



# Reading fasta
# ---

def read_protein_database(fname):
	n_tag_chars = 10
	with open(fname) as f:
		return {protein.rstrip(): meta[:n_tag_chars] for (meta, protein) in zip(f, f)}



# Reading mgf
# ---

def spec_data_parser(spec_data_str):
	return sorted([tuple(map(float, line.split())) for line in spec_data_str.splitlines()])

def mgf_parse_error(capture):
	raise ValueError('mgf parse error: {}'.format(capture))

parse_mgf_file = typed_regex_parser(r"""
BEGIN\ IONS\n
(?:(?:
	  (?: TITLE=index = (?P<index:int>.*))
	| (?: PEPMASS     = (?P<mz:str>.*))
	| (?: CHARGE      = (?P<z:int>.*)\+)
	| (?: SCANS       = (?P<n_scans:int>.*))
	| (?: RTINSECONDS = (?P<retention_time:float>.*))
	| (?P<err:mgf_parse_error>\D.*)
)\n)*
(?P<spec_data:spec_data>[^E]*)\n
END\ IONS
""", flags = re.MULTILINE | re.VERBOSE, custom_types = {'spec_data': spec_data_parser, 'mgf_parse_error': mgf_parse_error})

def read_mgf_file(fname, limit = None, keep_in_memory = False):
	with open(fname) as f:
		specs = (Spectrum(**record) for record in parse_mgf_file(f.read()))
		return (list if keep_in_memory else iter)(itertools.islice(specs, limit))

