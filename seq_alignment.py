import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

class InvalidPairException(Exception):
  pass

class Matrix:
  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)

  def _load_matrix(self, matrix_filename):
    with open(matrix_filename) as matrix_file:
      matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
      entries = row.split()
      row_name = entries.pop(0)
      matrix[row_name] = {}

      if len(entries) != len(columns):
        raise Exception('Improper entry number in row')
      for column_name in columns:
        matrix[row_name][column_name] = entries.pop(0)

    self._matrix = matrix

  def lookup_score(self, a, b):
    a = a.upper()
    b = b.upper()

    if a not in self._matrix or b not in self._matrix[a]:
      raise InvalidPairException('[%s, %s]' % (a, b))
    return float(self._matrix[a][b])

seq1 = sys.argv[1]
seq2 = sys.argv[2]
matrix_filename = "/work/langar2/scripts/weights/blosum62.txt"
matrix = Matrix(matrix_filename)

alignments = pairwise2.align.globalcs(seq1, seq2, matrix.lookup_score,-100,-100,penalize_end_gaps=False)#, score_only=True)
print format_alignment(*alignments[0]),