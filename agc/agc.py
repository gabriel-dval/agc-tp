#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import gzip
import textwrap
from pathlib import Path
from collections import Counter
from itertools import combinations
from typing import Iterator, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Gabriel Duval"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Gabriel Duval"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Gabriel Duval"
__email__ = "gabriel.duval@etu.u-paris.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, 'rt') as filout:
        current_seq = ''
        for line in filout:
            if line.strip().startswith('>') and len(current_seq) >= minseqlen:
                yield current_seq
                current_seq = ''
            elif line.strip().startswith('>') and len(current_seq) < minseqlen:
                current_seq = ''
                continue
            else:
                current_seq += line.strip()
        
        if len(current_seq) >= minseqlen:
            yield current_seq


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    all_seqs = [sequence for sequence in read_fasta(amplicon_file, minseqlen)]
    count_dict = Counter(all_seqs)
    count_dict = count_dict.most_common()

    for (key, count) in count_dict:
        if count >= mincount:
            yield[key, count]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    identities = 0
    for ci, cj in zip(alignment_list[0], alignment_list[1]):
        if ci == cj:
            identities += 1
    
    id = (identities/len(alignment_list[0])) * 100

    return id


def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    seqs_and_counts = [seq for seq in dereplication_fulllength(amplicon_file, minseqlen, mincount)]

    otu = [seqs_and_counts[0]]  # First sequence necessarily has most abundance
    for seq_and_count in seqs_and_counts:
        status = True
        for ref in otu:
            alignment = nw.global_align(seq_and_count[0], ref[0], gap_open=-1, gap_extend=-1, 
                                        matrix=str(Path(__file__).parent / "MATCH"))
            identity = get_identity(alignment)

            if identity > 97:
                status = False
                break
        
        if status:
            otu.append(seq_and_count)

    return otu


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, 'w') as filin:
        for i, (sequence, count) in enumerate(OTU_list):
            filin.write(f'>OTU_{i+1} occurrence:{count}\n')
            filin.write(f'{textwrap.fill(sequence, width=80)}\n')


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Args
    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = 0
    kmer_size = 0
    output_file = args.output_file
    
    # Programme
    otu_list = abundance_greedy_clustering(amplicon_file,
                                           minseqlen,
                                           mincount,
                                           chunk_size,
                                           kmer_size)
    
    # Save results
    write_OTU(otu_list, output_file)


if __name__ == '__main__':
    main()
