#!/usr/bin/env python2.7
# Chris Eisenhart chrisEisenhart1992@gmail.com

"""
"""


import sys
import operator
import argparse
import random

def parseArgs(args): 
    """ 
    Sets up the argparse command-line parser
    """
    parser = argparse.ArgumentParser(description = __doc__) # Some python magic to reference the docstring
    parser.add_argument ("--fa",
    help = " The input fa file",
    action = "store")
    parser.add_argument ("--fastq",
    help = " The output fastq file",
    action = "store")
    args = parser.parse_args() # Initializes the options specified above
    return args # Returns the parser as an object


def get_random_base():
    number = random.randint(0,3)
    if number == 0:
        return "A"
    elif number == 1:
        return "T"
    elif number == 2:
        return "C"
    elif number == 3:
        return "G"


def rev_comp_base(input_base):
    if input_base == "A": return "T"
    if input_base == "T": return "A"
    if input_base == "C": return "G"
    if input_base == "G": return "C"


def rev_comp(in_dna):
    result = ""
    for base in in_dna:
        result += rev_comp_base(base)
    return result[::-1] # Python list comprehension for reversing


def main(args):
    """
    """
    options = parseArgs(args)
    ref_seq = ""
    with open(options.fa, "r") as in_file:
        for line in in_file:
            if line.startswith(">"): continue # Skip headers/comments
            ref_seq += line.strip()

    r1_results = []
    r2_results = []
    for i in xrange(0, 1000):
        seq_len = 250
        start_base = random.randint(0, len(ref_seq) - seq_len)
        seq = ref_seq[start_base - 1 :start_base + seq_len - 1]
        header = ">{}".format(i)
        qual = "A" * seq_len
        r1_results.append(header + ".artificialseq R1")
        r1_results.append(seq)
        r1_results.append("+")
        r1_results.append(qual)
        r2_results.append(header + ".artificialseq R2")
        r2_results.append(rev_comp(seq))
        r2_results.append("+")
        r2_results.append(qual)

    with open(options.fastq + ".R1.fastq", "w") as out_file:
        out_file.write("\n".join(r1_results))

    with open(options.fastq + ".R2.fastq", "w") as out_file:
        out_file.write("\n".join(r2_results))


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
