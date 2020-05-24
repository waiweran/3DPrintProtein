sizes = {
	'H': 1.10,
	'C': 1.70,
	'N': 1.55,
	'O': 1.52,
	'F': 1.47,
	'LI': 0.90,
	'BE': 0.59,
	'P': 1.80,
	'S': 1.80,
	'CL': 1.67,
	'NA': 1.16,
	'MG': 0.86,
	'AL': 0.68,
	'K': 1.52,
	'CA': 1.14,
	'FE': 0.75,
	'NI': 0.75,
	'CU': 0.90,
	'ZN': 0.88,
	'AS': 1.85,
	'BR': 1.82,
}

def get_space_size(name):
	'''Gets atom radius for space-filling model, by element name'''
	if name not in sizes:
		print('WARNING: Did not know atomic radius of {}'.format(name))
		return 2
	return sizes[name]