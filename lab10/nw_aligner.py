#!/usr/bin/python

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""

import os
import sys

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.

        Example:

        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None
        reference_codes = []
        #count = 0;
        with open(fname) as fp:
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                    ### TODO ###
                    # Parse matrix file line-by-line and load into nested dictionaries.
                line = line.rstrip()
                if line.startswith("#"):
                    continue
                if 'A' in line:
                    for code in line.split(' '):
                        # if code == ' ' or code =='\n':
                            # continue
                        reference_codes.append(code)
                        score_matrix[code] = {}
                    for code in line.split():
                        for key in score_matrix.keys():
                            score_matrix[key][code] = 0;
                    #print(reference_codes)
                elif len(line.split()) < len(reference_codes):
                    gap_penalty = int(line)
                else:
                    curr = list(line.split())
                    for i in range(len(curr)):
                        score_matrix[reference_codes[i]][reference_codes[line_num-1]] = curr[i]
                    #count +=1
                # Last line of matrix contains the gap penalty which must be pulled
                # out and returned.
                #count = 0
                # print(score_matrix)
                # print(gap_penalty)
        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """

        seqs = []

        ### TODO ###
        # Load FASTA file and return list of sequences.
        with open(fname) as fp:
            for num, line in enumerate(fp):
                line = line.rstrip()
                if '>' in line:
                    continue;
                else:
                    seqs.append(line)
            if len(seqs) != 2:
                raise Exception('There were not 2 seqs in the fasta file')
        # Throw an error if there are more than two sequences in the file.
        # print(seqs)
        return seqs

    def align(self, seq_x, seq_y, print_matrix = True):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]

        ### TODO ###
        # Fill the top row of the matrix with scores for gaps.
        for i in range(len(matrix)):
            matrix[i][0] = self.gap_penalty
        # Fill the first column of the matrix with scores for gaps.
        for i in range((len(seq_y) + 1)):
            matrix[0][i] = self.gap_penalty
        #print(matrix)
        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        for x in range(1, len(seq_x)+1):
            for y in range(1, len(seq_y)+1):
                match_score = self.score_matrix[seq_x[x-1]][seq_y[y-1]]
                match_score = int(match_score)
                #print(match_score)
                ### TODO ###
                # Take the maximum score of three possibilities:
                #   1) The element in the matrix diagonal from this one
                #      plus the score of an exact match
                #   2) The element to the left plus a gap penalty
                #   3) The element above plus a gap penalty
                # ... and set the current element (matrix[x][y]) equal to that
                diagonal = int(matrix[x-1][y-1]) + match_score
                left = int(matrix[x-1][y]) + self.gap_penalty
                up = int(matrix[x][y-1]) + self.gap_penalty
                max_val = max(diagonal, left, up)
                matrix[x][y] = max_val
                # print(diagonal)
                # print(left)
                # print (up)
                # print('max is ' + str(max_val))
                if max_val == diagonal:
                    pointers[x][y]= 'D'
                elif max_val == left:
                    pointers[x][y] = 'U'
                else:
                    pointers[x][y]= 'L'
                # Keep track of which of these choices you made by setting
                #   the same element (i.e., pointers[x][y]) to some value that
                #   has meaning to you.

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print(" ".join(map(lambda i: str(int(i)), matrix[x])))

        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            # print(x)
            # print(y)
            move = pointers[x][y]
            if move == 'D':
                align_x.append(seq_x[x-1])
                align_y.append(seq_y[y-1])
                x -= 1
                y -= 1
            elif move == 'L':
                align_y.append('-')
                align_x.append(seq_x[x-1])
                x -= 1
            elif move == 'U':
                align_x.append('-')
                align_y.append(seq_y[y-1])
                y -= 1
            ### TODO ###
            # Follow pointers back through the matrix to the origin.
            # Depending on which "move" you made at each element in the
            #   matrix, you'll either align seq_x to seq_y, seq_x to a gap, or
            #   seq_y to a gap.

        # flip the alignments, as they're reversed
        return "".join(align_x[::-1]), "".join(align_y[::-1])

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

if __name__ == '__main__':
    def usage():
        print('usage: %s matrixfilename stringfilename')
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print('Can not open %s' % (fname,))
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
