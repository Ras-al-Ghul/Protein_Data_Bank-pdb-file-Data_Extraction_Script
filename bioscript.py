import math
import sys

# Returns modulus of the vector
def get_Vector_Modulus(vector):
 	return math.sqrt(math.fabs(vector[0] * vector[0]) + math.fabs(vector[1] * vector[1]) + math.fabs(vector[2] * vector[2]) )

# Returns dihedral angle when four points are passed as params
def get_Dihedral_Angle(a, b, c, d):
	# a,b,c,d are the points
	# Given the coordinates of the four points, obtain the vectors b1, b2, and b3 by vector subtraction.
	b1 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
	b2 = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
	b3 = (d[0] - c[0], d[1] - c[1], d[2] - c[2])

	# Compute normal vectors to plane - normal1 and normal2
	sub1 = 0
	sub2 = 0
	sub3 = 0

	# Cross product
	sub1 = b1[1]*b2[2] - b1[2]*b2[1]
	sub2 = b1[2]*b2[0] - b1[0]*b2[2]
	sub3 = b1[0]*b2[1] - b1[1]*b2[0]
	normal1 = (sub1, sub2, sub3)

	# Cross product
	sub1 = b2[1]*b3[2] - b2[2]*b3[1]
	sub2 = b2[2]*b3[0] - b2[0]*b3[2]
	sub3 = b2[0]*b3[1] - b2[1]*b3[0]
	normal2 = (sub1, sub2, sub3)

	# Compute unit vector along normal1 and normal2
	modnormal1 = get_Vector_Modulus(normal1)
	normal1 = (normal1[0]/modnormal1, normal1[1]/modnormal1, normal1[2]/modnormal1)
	
	modnormal2 = get_Vector_Modulus(normal2)
	normal2 = (normal2[0]/modnormal2, normal2[1]/modnormal2, normal2[2]/modnormal2)

	# Compute unit vector along b2
	modb2 = get_Vector_Modulus(b2)
	b2 = (b2[0]/modb2, b2[1]/modb2, b2[2]/modb2)

	# normal1 X b2
	sub1 = normal1[1]*b2[2] - normal1[2]*b2[1]
	sub2 = normal1[2]*b2[0] - normal1[0]*b2[2]
	sub3 = normal1[0]*b2[1] - normal1[1]*b2[0]
	m1 = (sub1,sub2,sub3)

	x = normal1[0]*normal2[0] + normal1[1]*normal2[1] + normal1[2]*normal2[2]
	y = m1[0]*normal2[0] + m1[1]*normal2[1] + m1[2]*normal2[2]

	# Returns dihedral angle with correct sign
	return math.atan2(y, x) * 180/math.pi

def main():
	infile = open(sys.argv[1], 'rU')
	outfile = open(sys.argv[1].split('.')[0] + '_output.txt', 'w')
	
	chains = {}
	aminoacids = {}
	ligands = {}
	angles = []
	temp1 = {}
	temp2 = {}

	index_of_angle = 0
	acid_count = 0
	unknown_count = 0
	
	angle_tuple = (0,0,0)
	cur_chain = None

	temp1['N'] = temp2['N'] = (float('nan'), float('nan'), float('nan'))
	temp1['CA'] = temp2['CA'] = (float('nan'), float('nan'), float('nan'))
	temp1['C'] = temp2['C'] = (float('nan'), float('nan'), float('nan'))

	arr = infile.readlines()
	i = -1

	# Parse the pdb file and update variables
	while i < len(arr) - 1:
		i += 1
		line = arr[i]
		if line.startswith('TITLE'):
			name = line[5:].strip()
		if line.startswith('SEQRES'):
			list_from_line = line.split()
			chain_name = list_from_line[2].strip()
			if chain_name not in chains:
				chains[chain_name] = int(list_from_line[3])
				acid_count += chains[chain_name]
			for temps in list_from_line[4:]:
				temps = temps.strip()
				if temps not in aminoacids:
					aminoacids[temps] = 0
				aminoacids[temps] += 1
		if line.startswith('HET') and not line.startswith('HETATM'):
			list_from_line = line.split()
			ligand = list_from_line[1].strip()
			if ligand not in ligands:
				ligands[ligand] = 1
				
		if line.startswith('ATOM'):
			list_from_line = line.split()
			if list_from_line[4].strip() == cur_chain:
				if index_of_angle == 1:
					temp1['N'] = temp2['N']
					temp1['CA'] = temp2['CA']
					temp1['C'] = temp2['C']
					while True:
						if line.startswith('ATOM') == False or list_from_line[4].strip() != cur_chain:
							angles.append(('psi', float('nan')))
							angles.append(('omega', float('nan')))
							index_of_angle = (index_of_angle + 2) % 3
							break
						if list_from_line[2].strip() == 'N':
							temp2['N'] = (float(list_from_line[6].strip()), float(list_from_line[7].strip()), float(list_from_line[8].strip()))
							psi = get_Dihedral_Angle(temp1['N'], temp1['CA'], temp1['C'], temp2['N'])
							angle_tuple = ('psi', psi)
							angles.append(angle_tuple)
							index_of_angle = (index_of_angle + 1)%3
							break
						i += 1
						line = arr[i]
						list_from_line = line.split()
					continue

				if index_of_angle == 2:
					temp2['CA'] = (float(list_from_line[6].strip()), float(list_from_line[7].strip()), float(list_from_line[8].strip()))
					omega = get_Dihedral_Angle(temp1['CA'], temp1['C'], temp2['N'], temp2['CA'])
					angle_tuple = ('omega', omega)
					angles.append(angle_tuple)
					index_of_angle = (index_of_angle + 1)%3
					continue

				if index_of_angle == 0:
					temp2['C'] = (float(list_from_line[6].strip()), float(list_from_line[7].strip()), float(list_from_line[8].strip()))
					phi = get_Dihedral_Angle(temp1['C'], temp2['N'], temp2['CA'], temp2['C'])
					angle_tuple = ('phi', phi)
					angles.append(angle_tuple)
					index_of_angle = (index_of_angle + 1)%3
					continue

			else:
				cur_chain = list_from_line[4].strip()
				temp2['N'] = (float(list_from_line[6].strip()), float(list_from_line[7].strip()), float(list_from_line[8].strip()))
				i += 1
				line = arr[i]
				list_from_line = line.split()
				temp2['CA'] = (float(list_from_line[6].strip()), float(list_from_line[7].strip()), float(list_from_line[8].strip()))
				i += 1
				line = arr[i]
				list_from_line = line.split()
				temp2['C'] = (float(list_from_line[6].strip()), float(list_from_line[7].strip()), float(list_from_line[8].strip()))
				
				angle_tuple = ('phi', float('nan'))
				angles.append(angle_tuple)
				index_of_angle = (index_of_angle + 1)%3

	if 'UNK' in aminoacids:
		unknown_count = aminoacids['UNK']
		del aminoacids['UNK']

	# Printing
	outfile.write(name + '\n')
	outfile.write('LENGTH  '+  str(acid_count)+'\n')
	outfile.write('CHAINS  ' + str(len(chains.keys())) +'\t' + ','.join(sorted(chains.keys())) +'\n')
	acid_count = float(acid_count)
	for k in sorted(aminoacids.keys()):
		outfile.write(k + '\t' + str(aminoacids[k]/acid_count) + '\n')
	outfile.write('UNKNOWN ' + str(unknown_count) + '\n')
	outfile.write('LIGANDS ' + ','.join(sorted(ligands.keys())) + '\n')
	i = 0
	for k in sorted(chains.keys()):
		outfile.write('CHAIN-' + str(k) + '\n')
		while i < len(angles)-1:
			phi = angles[i][1]
			psi = angles[i+1][1]
			omega = angles[i+2][1]
			if math.isnan(phi):
				phi = 'NA'
			if math.isnan(psi):
				psi = 'NA'
			if math.isnan(omega):
				omega = 'NA'
			outfile.write(str(phi) + '\t\t' + str(psi) + '\t\t' + str(omega) + '\n')
			i += 3
			if psi == 'NA' and omega == 'NA':
				break

if __name__ == '__main__':
	main()