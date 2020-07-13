import glob

from Bio import motifs
from Bio import SeqIO
from functools import reduce


def get_best_score(seq, ppm):
    """
    This function calculates the highest score
    that can be obtained for the sequence based on the ppm
    :param seq: sequence
    :param ppm: position probability matrix
    :return:
    """
    n_seq = len(seq)
    n_ppm = ppm.length
    if n_ppm < n_seq:
        return 0

    score = 0
    for shift in range(n_ppm - n_seq + 1):
        score_pos = [ppm[s][i + shift] for i, s in enumerate(seq)]
        score = max(score, reduce((lambda x, y: x * y), score_pos))
    return score


path_jaspar = 'jaspar/'
path_motifs = 'motifs/'

# Read ppm in jaspar format
files_jaspar = glob.glob(path_jaspar + '*.jaspar')

ppms = []
pfms = []
for f in files_jaspar:
    with open(f) as fh:
        pfms += [motifs.read(fh, "jaspar")]
        ppms += [motifs.matrix.PositionWeightMatrix('ACGT', pfms[-1].counts)]

# Read motifs
files_motifs = glob.glob(path_motifs + '*.fasta')
records = [item for f in files_motifs for item in list(SeqIO.parse(f, "fasta"))]


# Get best ppm for each motif
best_motifs = []
for record in records:
    ppm_scores = [get_best_score(record, ppm) for ppm in ppms]
    max_id = ppm_scores.index(max(ppm_scores))
    best_motifs += [[pfms[max_id].matrix_id, pfms[max_id].name, max(ppm_scores)]]



