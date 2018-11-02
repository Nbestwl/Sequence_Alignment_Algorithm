# Author: Lei Wang
# Date: 10/12/2018
# EECS 730

import sys
from math import ceil

SCORING_MATRIX_FILE = "BLOSUM62.txt"
GAP_OPENNING = -11
GAP_EXTENSION = -1
DIAGONAL = 1
LEFT = 2
TOP = 3

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def read_score_matrix():
	label, scores = [], []
	with open(SCORING_MATRIX_FILE) as fp:
		l = [line for line in fp]
		label = l[0][:-2].split()
		scores = [map(int, line.split()[1:-1]) for line in l[1:-1]]
	return label, scores


def print_matrix(x, y, A):
    """Print the matrix with the (0,0) entry in the top left
    corner. Will label the rows by the sequence and add in the
    0-row if appropriate."""

    # decide whether there is a 0th row/column
    if len(x) == len(A):
        print "%5s" % (" "),
    else:
        print "%5s %5s" % (" ","*"),
        y = "*" + y

    # print the top row
    for c in x:
        print "%5s" % (c),
    print

    for j in xrange(len(A[0])):
        print "%5s" % (y[j]),
        for i in xrange(len(A)):
            print "%5.0f" % (A[i][j]),
        print


def matrix_init(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in range(sizex)]


def find_score(x, y, sequences, label, scores):
	letter1 = sequences[0][x-1]
	letter2 = sequences[1][y-1]
	index1 = label.index(letter1.upper())
	index2 = label.index(letter2.upper())
	score = scores[index1][index2]
	return score

def local_alignment(sequences, label, scores):
	local_sequence, sequence1, sequence2, symbols = [], [], [], []
	optiaml_score = float('-inf')
	optloc = (0,0)

	col = len(sequences[0])
	row = len(sequences[1])
	A = matrix_init(col + 1, row + 1)
	trace_back = matrix_init(col + 1, row + 1)

	# fill the table
	for j in range(1, row + 1):
		for i in range(1, col + 1):
			score = find_score(i, j, sequences, label, scores)
			A[i][j] = max(
				A[i][j-1] + GAP_OPENNING if trace_back[i][j-1] == 1 else A[i][j-1] + GAP_EXTENSION,
				A[i-1][j] + GAP_OPENNING if trace_back[i-1][j] == 1 else A[i-1][j] + GAP_EXTENSION,
				A[i-1][j-1] + score,
				0
			)
			if (A[i][j-1] + GAP_OPENNING == A[i][j]) or (A[i][j-1] + GAP_EXTENSION == A[i][j]):
				trace_back[i][j] = 2
			elif (A[i-1][j] + GAP_OPENNING == A[i][j]) or (A[i-1][j] + GAP_EXTENSION == A[i][j]):
				trace_back[i][j] = 3
			elif A[i-1][j-1] + score != 0 and A[i][j] == 0:
				trace_back[i][j] = 0
			else:
				trace_back[i][j] = 1

			if A[i][j] >= optiaml_score:
				optiaml_score = A[i][j]
				optloc = (i,j)

	# print "\n\nfilled matrix"
	# print_matrix(sequences[0], sequences[1], A)
	# print "\n\ntrack back matrix"
	# print_matrix(sequences[0], sequences[1], trace_back)
	# print "Max location in matrix =", optloc

	temp = (sequences[0][optloc[0]-1], sequences[1][optloc[1]-1])
	local_sequence.append(temp)
	while not (optloc[0]-1 == 0 and optloc[1]-1 == 0):
		direction = trace_back[optloc[0]][optloc[1]]
		if direction == 1:
			optloc = (optloc[0] - 1, optloc[1] - 1)
		elif direction == 2:
			optloc = (optloc[0], optloc[1] - 1)
		else:
			optloc = (optloc[0] - 1, optloc[1])
		temp = (sequences[0][optloc[0]-1], sequences[1][optloc[1]-1])
		local_sequence.append(temp)

	unzipped = zip(*local_sequence)
	sequence1 = list(unzipped[0])[::-1]
	sequence2 = list(unzipped[1])[::-1]

	for i in range(len(sequence1)):
		if sequence1[i] == sequence2[i]:
			symbols.append('|')
		elif sequence1[i] != '-' and sequence2[i] != '-' and sequence1[i] != sequence2[i]:
			symbols.append('*')
		else:
			symbols.append(' ')

	print "\n\nlocal alignment:"
	print "Score =", optiaml_score
	times = int(ceil(len(sequence1) / 50.0))
	for i in range(times):
		print '\n'
		if i + 1 != times:
			print ' '.join(sequence1[i*50:(i*50+50)])
			print ' '.join(symbols[i*50:(i*50+50)])
			print ' '.join(sequence2[i*50:(i*50+50)])
		else:
			print ' '.join(sequence1[i*50:])
			print ' '.join(symbols[i*50:])
			print ' '.join(sequence2[i*50:])


def global_alignment(sequences, label, scores):
	local_sequence, sequence1, sequence2, symbols = [], [], [], []
	optloc = (0,0)

	col = len(sequences[0])
	row = len(sequences[1])
	A = matrix_init(col + 1, row + 1)
	trace_back = matrix_init(col + 1, row + 1)

	# fill the table
	for i in range(row + 1):
		A[0][i] = 0 + GAP_EXTENSION * i
		A[i][0] = 0 + GAP_EXTENSION * i

	for j in range(1, row + 1):
		for i in range(1, col + 1):
			score = find_score(i, j, sequences, label, scores)
			A[i][j] = max(
				A[i][j-1] + GAP_OPENNING if trace_back[i][j-1] == 1 else A[i][j-1] + GAP_EXTENSION,
				A[i-1][j] + GAP_OPENNING if trace_back[i-1][j] == 1 else A[i-1][j] + GAP_EXTENSION,
				A[i-1][j-1] + score
			)
			if (A[i][j-1] + GAP_OPENNING == A[i][j]) or (A[i][j-1] + GAP_EXTENSION == A[i][j]):
				trace_back[i][j] = 2
			elif (A[i-1][j] + GAP_OPENNING == A[i][j]) or (A[i-1][j] + GAP_EXTENSION == A[i][j]):
				trace_back[i][j] = 3
			else:
				trace_back[i][j] = 1

		optiaml_score = A[i][j]
		optloc = (i,j)

	# print "\n\nfilled matrix"
	# print_matrix(sequences[0], sequences[1], A)
	# print "\n\ntrack back matrix"
	# print_matrix(sequences[0], sequences[1], trace_back)
	# print "Max location in matrix =", optloc

	temp = (sequences[0][optloc[0]-1], sequences[1][optloc[1]-1])
	local_sequence.append(temp)
	while not (optloc[0]-1 == 0 and optloc[1]-1 == 0):
		direction = trace_back[optloc[0]][optloc[1]]
		if direction == 1:
			optloc = (optloc[0] - 1, optloc[1] - 1)
		elif direction == 2:
			optloc = (optloc[0], optloc[1] - 1)
		else:
			optloc = (optloc[0] - 1, optloc[1])
		temp = (sequences[0][optloc[0]-1], sequences[1][optloc[1]-1])
		local_sequence.append(temp)

	unzipped = zip(*local_sequence)
	sequence1 = list(unzipped[0])[::-1]
	sequence2 = list(unzipped[1])[::-1]

	for i in range(len(sequence1)):
		if sequence1[i] == sequence2[i]:
			symbols.append('|')
		elif sequence1[i] != '-' and sequence2[i] != '-' and sequence1[i] != sequence2[i]:
			symbols.append('*')
		else:
			symbols.append(' ')

	print "\n\nglobal alignment"
	print "Score =", optiaml_score
	times = int(ceil(len(sequence1) / 50.0))
	for i in range(times):
		print '\n'
		if i + 1 != times:
			print ' '.join(sequence1[i*50:(i*50+50)])
			print ' '.join(symbols[i*50:(i*50+50)])
			print ' '.join(sequence2[i*50:(i*50+50)])
		else:
			print ' '.join(sequence1[i*50:])
			print ' '.join(symbols[i*50:])
			print ' '.join(sequence2[i*50:])


if __name__ == '__main__':
	filenames, sequences, names = [], [], []
	if len(sys.argv) == 3:
	    filenames.append(sys.argv[1])
	    filenames.append(sys.argv[2])
	else:
	    print "Usage: python " + sys.argv[0] + " [sequence1] [sequence2]"

	for filename in filenames:
		with open(filename) as fp:
	   		for (name, seq) in read_fasta(fp):
	   			names.append(name)
	   			sequences.append(seq)

	label, scores = read_score_matrix()
	local_alignment(sequences, label, scores)
	global_alignment(sequences, label, scores)



