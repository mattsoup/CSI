#!/usr/bin/env python
"""
This script calculates the HII (published as CSI) using the .tbl files
generated after annotating Hodgkinia.

"""

import sys
import re

if len(sys.argv) != 3:
        print "Usage: HII.py <list of .tbl files> <list of fasta files>\n"
        quit()

tbl_dict = {}
tbl_list_file = open(sys.argv[1], "r")
fasta_list_dict = {}
fasta_list_file = open(sys.argv[2], "r")

# This will identify the species and corresponding .tbl files
for line in tbl_list_file:
    if "Reference" in line:
        regex = re.match("Reference: (.*?)\t(.*?)\n", line)
        ref = regex.group(1)
        tbl_dict[ref] = regex.group(2)
    else:
        regex = re.match("(.*?)\t(.*?)\n", line)
        sp = regex.group(1)
        tbl_dict[sp] = regex.group(2)
tbl_list_file.close()

# This will identify the species and corresponding fasta files
for line in fasta_list_file:
    regex = re.match("(.*?)\t(.*?)\n", line)
    sp = regex.group(1)
    file_loc = regex.group(2)
    fasta_list_dict[sp] = file_loc


def get_genes(species, file):
    """ This will go through a .tbl file and extract all protein/rRNA genes"""

    scaffold_dict = {}
    HGC_set = set([])
    HGC_list = []
    cicada_set = set([])
    file = open(file, "r")
    tbl = file.readlines()
    for x in range(0, len(tbl)):
        if tbl[x].startswith(">"):
            if len(HGC_list) != 0:
                HGC_list.sort()
                HGC_genes = ""
                for item in HGC_list:
                    HGC_genes += "%s_" % item
                HGC_set.add(HGC_genes)
                HGC_list = []
            regex = re.match(">Feature (.*?)\n", tbl[x])
            header = regex.group(1)
            scaffold_dict[header] = set([])
        elif (tbl[x].startswith("\t\t\tgene\t") and "pseudo" not in tbl[x + 2]
             and "tRNA" not in tbl[x] and
             ("CDS" in tbl[x-1] or "CDS" in tbl[x-2] or "rRNA" in tbl[x-1])):
            regex = re.match("\t\t\tgene\t(.*?)\n", tbl[x])
            gene = regex.group(1)
            scaffold_dict[header].add(gene)
            cicada_set.add(gene)
            HGC_list.append(gene)
    file.close()
    return scaffold_dict, HGC_set, cicada_set


def get_fasta(species, file):
    """This will go through all fasta file and put each entry into a dictionary"""

    temp_dict = {}
    file = open(file, "r")
    for line in file:
        if line.startswith(">"):
            header = line[1:].strip()
            temp_dict[header] = ""
        else:
            temp_dict[header] += line.strip()
    return temp_dict


def compare_genes(ref_scaffold, species, scaffold_dict):
    """
    This will compare gene presence/absence for one reference scaffold
    against all scaffolds from a target species

    """

    max = 0
    scaffold_dict = scaffold_dict
    closest_scaffold = ""
    num_ref_genes = len(ref_dict[ref][ref_scaffold])
    HII = 0
    print "\n\n\n\n#########################################################\n\
                   Working on %s" % ref_scaffold
    for scaffold in scaffold_dict[species]:
        num_scaffold_genes = len(scaffold_dict[species][scaffold])
        homologues = ref_dict[ref][ref_scaffold].intersection(scaffold_dict[species][scaffold])
        if (num_ref_genes + num_scaffold_genes) != 0:
            homologue_score = check_sister(homologues, ref_scaffold, scaffold, species)
            print "#######################################################\n\
                   Homologue score of %s is %s" % (scaffold, homologue_score)
            Jaccard = float(len(ref_dict[ref][ref_scaffold].intersection(scaffold_dict[species][scaffold]))) / len(ref_dict[ref][ref_scaffold].union(scaffold_dict[species][scaffold]))
            HII = Jaccard
            length_scalar = compare_length(ref, species, ref_scaffold, scaffold)
            HII = HII * length_scalar
            print "HII of %s is %s" % (scaffold, HII)
        if abs(HII) > abs(max):
            max = HII
            print HII, abs(HII), max
            closest_scaffold = scaffold
    print "The closest related circle of %s %s in %s is %s, HII = \
           %s" % (ref, ref_scaffold, species, closest_scaffold, max)
    return [closest_scaffold, max]


def compare_HGC(species):
    """
    This calculates the 'HGC HII' for two species...in the end I did not use
    this for the paper, and it currently is not called in the script.

    """
    print "doing HGC comparison of %s and %s" % (ref, species)
    print ref_HGC_dict[ref], HGC_dict[species]
    Jaccard = float(len(ref_HGC_dict[ref].intersection(HGC_dict[species]))) / len(ref_HGC_dict[ref].union(HGC_dict[species]))
    return Jaccard


def compare_cicada(species):
    """
    This computes the 'cicada HII' by comparing the total gene sets of the HGC
    between two cicada species.

    """
    print "doing cicada comparison of %s and %s" % (ref, species)
    Jaccard = float(len(ref_cicada_dict[ref].intersection(cicada_dict[species]))) / len(ref_cicada_dict[ref].union(cicada_dict[species]))
    return Jaccard


def compare_length(ref, species, ref_scaffold, scaffold):
    """
    This compares the length of two scaffolds (a reference and a target) to
    scale the HII score of those same scaffolds


    """
    ref_scaffold_length = len(fasta_dict[ref][ref_scaffold])
    scaffold_length = len(fasta_dict[species][scaffold])
    scalar = float(ref_scaffold_length) / scaffold_length
    if scalar <= 1:
        print "Length score is %s" % scalar
        return scalar
    else:
        print "Length score is %s (%s)" % (scalar, (1 / scalar))
        return (1 / scalar)


# This calls 'get_genes' for each species
ref_dict = {}
ref_cicada_dict = {}
ref_HGC_dict = {}
scaffold_dict = {}
cicada_dict = {}
HGC_dict = {}
for item in tbl_dict:
    if item == ref:
        ref_dict[ref], ref_HGC_dict[ref], ref_cicada_dict[ref] = get_genes(ref, tbl_dict[ref])
    else:
        scaffold_dict[item], HGC_dict[item], cicada_dict[item] = get_genes(item, tbl_dict[item])

# This calls 'get_fasta' for each species
fasta_dict = {}
for species in fasta_list_dict:
    fasta_dict[species] = get_fasta(species, fasta_list_dict[species])

# This calls 'compare_genes' for each scaffold in the reference species
comp_dict = {}
for ref_scaffold in ref_dict[ref]:
    comp_dict[ref_scaffold] = {}
    comp_dict[ref_scaffold][ref] = []
    comp_dict[ref_scaffold][ref] = compare_genes(ref_scaffold, ref, ref_dict)
    for species in scaffold_dict:
        comp_dict[ref_scaffold][species] = []
        comp_dict[ref_scaffold][species] = compare_genes(ref_scaffold, species, scaffold_dict)

# This calls "compare_HGC" and "compare_cicada" for all target species
comp_cicada_dict = {}
comp_HGC_dict = {}
for species in scaffold_dict:
    comp_cicada_dict[species] = compare_cicada(species)

# This prints the highest HGC and cicada HII for all target species
out_cicada = open("%s_cicada_HII.txt" % ref, "w")
out_HGC = open("%s_HGC_HII.txt" % ref, "w")
out_cicada.write("Reference species\tTarget species\tCicada HII\n")
out_HGC.write("Reference species\tTarget species\tHGC HII\n")
for species in scaffold_dict:
    cicada_HII = comp_cicada_dict[species]
    print "The cicada HII score between %s and %s is %s" % (ref, species, cicada_HII)
    out_cicada.write("%s\t%s\t%s\n" % (ref, species, cicada_HII))

# This prints out the highest HII for each reference scaffold
out = open("%s_HII.txt" % ref, "w")
out.write("Reference species\tReference scaffold\tTarget species\tClosest\
           scaffold\tHII\n")
for scaffold in comp_dict:
    for species in comp_dict[scaffold]:
        closest_scaffold = comp_dict[scaffold][species][0]
        HII = comp_dict[scaffold][species][1]
        print "The closest related circle of %s in %s is %s, HII =\
               %s" % (scaffold, species, closest_scaffold, HII)
        out.write("%s\t%s\t%s\t%s\t%s\n" % (ref, scaffold, species, closest_scaffold, HII))
