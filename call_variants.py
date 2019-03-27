#!/usr/bin/env python2.7
# Chris Eisenhart chrisEisenhart1992@gmail.com

"""
"""


import sys
import operator
import argparse


def parseArgs(args): 
    """ 
    Sets up the argparse command-line parser
    """
    parser = argparse.ArgumentParser(description = __doc__) # Some python magic to reference the docstring
    parser.add_argument ("--fa",
    help = " The reference file",
    action = "store")
    parser.add_argument ("--bam",
    help = " The input (s)bam file",
    action = "store")
    parser.add_argument ("--vcf",
    help = " The output vcf file",
    action = "store")
    args = parser.parse_args() # Initializes the options specified above
    return args # Returns the parser as an object


class ref_base_pos:
    def __init__(self, base):
        self.ref_base = base
        self.base_qual_list = [] # list of tuples, (base, quality) corresponding to the bases seen in aligments
        self.insertion_list = [] # list of base_qual_lists, one per base in the insertion


    def concensus_pileup(self):
        total_A, total_T, total_C, total_G = 0.1,0.1,0.1,0.1
        qual_A, qual_T, qual_C, qual_G = 0,0,0,0
        for tup in self.base_qual_list:
            base = tup[0]
            qual = ord(tup[1]) - 33 # convert phred char into a numeric score
            if base == "A":
                total_A += 1
                qual_A += qual
            elif base == "T":
                total_T += 1
                qual_T += qual
            elif base == "G":
                total_G += 1
                qual_G += qual
            elif base == "C":
                total_C += 1
                qual_C += qual
        avg_A = qual_A/total_A
        avg_T = qual_T/total_T
        avg_C = qual_C/total_C
        avg_G = qual_G/total_G

        sorted_pileups = sorted([("A", avg_A), ("T", avg_T), ("C", avg_C), ("G", avg_G)],
                                key = operator.itemgetter(1),
                                reverse = True)
        return sorted_pileups[0][0]


def parse_md(md_string):
    """
    md_strings represent how the sequence aligns base by base. They can be simple, for instance
    250M would indicate that 250 bases match. They can also be complicated, ex 10S120M10D20M2I40M
    would suggest 10 skipped bases, 120 match, 10 deleted, 20 match, 2 insertions and 40 matches.

    The sum of all bases in the md_string is the length of the input sequence.

    Return a list of tuples representing the entire MD string;
        [(250, M)]
            or
        [(10, S), (120, M), (10, D), (20, M), (2, I), (40, M)]

    """
    result = []
    count = ""
    for char in md_string:
        if char.isdigit(): count += char
        else:
            result.append((count, char))
            count = ""
    return result


def main(args):
    """
    """
    options = parseArgs(args)

    # Load the ref file. TODO: Split on chroms, put in multithreading
    ref = ""
    with open(options.fa, "r") as ref_file:
        for line in ref_file:
            if line.startswith(">"): continue
            ref += line.strip()

    # Build the data structure (for one chrom at a time)
    data_struct = [] # This needs a way better name
    for base in ref:
        data_struct.append(ref_base_pos(base))

    for line in open(options.bam, "r"):
        if line.startswith("@"): continue
        split_line = line.strip().split()
        start_base = int(split_line[3])
        cigar = split_line[5]
        base_string = split_line[9]
        qual_string = split_line[10]
        parse_md_list = parse_md(split_line[12]) # List of tuples [(10, S), (150, M), (10, S)]
        count = 0
        for base, qual in zip(base_string, qual_string):
            count += 1
            data_struct[start_base + count - 2].base_qual_list.append((base, qual))
    
    count = 0
    for base in data_struct:
        count += 1
        called_base = base.concensus_pileup()
        if called_base is not base.ref_base:
            print ("Found a variant {} at pos {}, expected {}".format(called_base, count, base.ref_base))


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
