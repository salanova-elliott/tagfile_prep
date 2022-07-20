#! /home/lukedane/miniconda3/bin/python

import sys
import argparse
import pprint
import string
import contextlib

'''
UPDATED 2022-02-03

TO DO
Strip leading characters from input files for library names (so you can use full filepaths)
Clarify --notused usage
Provide useful summary of data loaded in tagfile
Take functionality from SplitTagTable
'''

# {primer_name : [forward_primer, rev_primer]}
primer_dict = {"gh" : ["GGGCAATCCTGAGCCAA", "CCATTGAGTCTCTGCACCTATC"], 
            "12s" : ["ACACCGCCCGTCACCCT", "GTAYRCTTACCWTGTTACGACTT"], 
            "16s" : ["CGAGAAGACCCTATGGAGCT", "CCGAGGTCRCCCCAACC"]}

# Argument parser
parser = argparse.ArgumentParser(description = "Creates a single tagfile from multiple libraries")
parser.add_argument("-t", "--tags", metavar="tags", dest="tags", type=str, 
                    help="Tag file", required=True)
parser.add_argument("-l", "--libraries", metavar="libraries", dest="libraries",
                    nargs="+", help="List of libraries to add.", required=True)
parser.add_argument("-p", "--primer", metavar="primer", dest="primer",
                    type=str, help="Primer used in sequencing (gh, 12S, 16S, or edit dictionary in script)")
parser.add_argument("-o", "--outdir", metavar="outdir", dest="outdir", const=".", nargs="?",
                    type=str, help="Name / filepath of output dir", required=True)
parser.add_argument("-n", "--notused", metavar="notlibs", dest="notlibs", nargs="+", 
                    help="List of libs with no samples for the primer")
args = parser.parse_args()

# Accepted characters for sample name
accepted_chars = set(string.ascii_lowercase + string.ascii_uppercase + "_-0123456789")

# Nested dictionary of tags {plate_num : {well_coord : tag}}
tag_dict = {}

# Nested dictionary of libraries {lib_name : {plate_num: {well_coord: sample_name}}}
lib_dict = {}

# Nested dictionary of sample name counts {lib_name : {sample_name : count}}
lib_count_dict = {}
lib_count_dict["not_used"] = 1

# Populates tag dictionary
with open(args.tags) as f:

    assert f.readline().split("\t")[0] != "1", "Tag file should have headers"

    for l in f:
        split_line = l.rstrip().split("\t")
        if split_line[0] not in tag_dict:
            tag_dict[split_line[0]] = {}
        tag_dict[split_line[0]][split_line[1]] = split_line[2] 

    print(f"Loaded {len(tag_dict)} tag plates", file=sys.stderr)

# Checks that sample name is kosher
def name_check(sample_name):
    assert set(sample_name) <= accepted_chars, f"{sample_name} contains invalid characters"

# Function for loading libraries into dictionary
def library_loader(lib_file):

    # Names library after filename (without extension)
    lib_name = lib_file.split("/")[-1].split(".")[0]
    lib_dict[lib_name] = {}

    #Parses through library file
    with open(lib_file) as f:

        # Starts plate count at zero
        plate_num = 0

        for l in f:
            # "Primer plate" marks where plate layout starts in spreadsheet
            if "Primer plate" in l:
                plate_num = int(l.split("\t")[0].split(" ")[-1])
                lib_dict[lib_name][plate_num] = {}

            # "Extract plate" is at the bottom of spreadsheet. This kills the loading
            if "Extract plate" in l:
                break

            # If reached plate layout section and there are samples present
            if plate_num > 0 and l.split("\t")[1]:
                # Important not to rstrip() here since it will remove empty wells
                split_line = l.split("\t")
                # Grabs well row letter
                well_row = split_line[1]

                # Loops through each sample in row
                for i, sample_name in enumerate(split_line[2:]):
                    # If it's not the end of a row
                    if i < 12:
		                # Removes any possible spaces after sample names
                        sample_name = sample_name.rstrip()
                        # Checks sample name for valid characters
                        name_check(sample_name)
                        # Only adds to dictionary if sample name is present
                        if sample_name:
                            lib_dict[lib_name][plate_num][f"{well_row}{i+1}"] = sample_name
        
# Adds repeat number to each sample name
def add_repeat(lib_name):
    
    lib_count_dict[lib_name] = {}

    for plate in lib_dict[lib_name]:

        for well_coord in lib_dict[lib_name][plate]:
            sample_name = lib_dict[lib_name][plate][well_coord]

            if sample_name not in lib_count_dict[lib_name]:
                lib_count_dict[lib_name][sample_name] = 1
            else:
                lib_count_dict[lib_name][sample_name] += 1
            
            lib_dict[lib_name][plate][well_coord] += f"_rpt{str(lib_count_dict[lib_name][sample_name])}"

# EXECUTION
# Loops through libraries and loads them
for libr in args.libraries:
    print(f"Loading library: {libr}", file=sys.stderr)
    library_loader(libr)
# Adds repeat info to libraries in dictionary
for lib in lib_dict:
    add_repeat(lib)

# Displays info on how many repeats per sample there are
pprint.pprint(lib_count_dict)

# Output real library files
for outl in args.libraries:
    out_filename = f"{args.outdir}/{outl.split('.')[0]}_{args.primer}_tag.tsv"
    with open(out_filename, "w") as outfile:
        outfile.write("#exp\tsample\ttags\tforward_primer\treverse_primer\n")     
        lib_dict_key = outl.split("/")[-1].split(".")[0]
        # Loops through wells
        for plate in range(1,5):
            for num in range(1,13):
                for let in "ABCDEFGH":
                    current_well = f"{let}{num}"                
                    # Writes line info
                    if plate in lib_dict[lib_dict_key] and current_well in lib_dict[lib_dict_key][plate]:
                        outfile.write(f"{lib_dict_key}\t{lib_dict[lib_dict_key][plate][current_well]}\t{tag_dict[str(plate)][current_well][-8:]}\t{primer_dict[args.primer][0]}\t{primer_dict[args.primer][1]}\n")
                    else:
                        outfile.write(f"{lib_dict_key}\tnot_used{str(lib_count_dict['not_used'])}\t{tag_dict[str(plate)][current_well][-8:]}\t{primer_dict[args.primer][0]}\t{primer_dict[args.primer][1]}\n")
                        lib_count_dict['not_used'] += 1
# Output empty 
if args.notlibs:
    for nlib in args.notlibs:
        out_filename = f"{args.outdir}/{nlib.split('.')[0]}_{args.primer}_tag.tsv"
        with open(out_filename, "w") as outfile:
            outfile.write("#exp\tsample\ttags\tforward_primer\treverse_primer\n")
            lib_dict_key = outl.split("/")[-1].split(".")[0]
            # Loops through wells
            for plate in range(1,5):
               for num in range(1,13):
                   for let in "ABCDEFGH":
                       current_well = f"{let}{num}"
                       outfile.write(f"{lib_dict_key}\tnot_used{str(lib_count_dict['not_used'])}\t{tag_dict[str(plate)][current_well][-8:]}\t{primer_dict[args.primer][0]}\t{primer_dict[args.primer][1]}\n")
                       lib_count_dict['not_used'] += 1
'''      
# OUTPUT
tag_num = 1
with open(args.outfile, "w") as f:
    
    # Headers
    f.write("Tag_number\tTag")
    for libn in lib_dict:
        f.write(f"\t{libn}")
    if args.notlibs:
        for libn in args.notlibs:
            f.write(f"\t{libn.split('.')[0]}")
    f.write("\n")

    # Loops through wells
    for plate in range(1,5):
        for num in range(1,13):
            for let in "ABCDEFGH":

                current_well = f"{let}{num}"
                
                # Writes tag number and tag sequence
                f.write(f"{tag_num}\t{tag_dict[str(plate)][current_well]}")
                tag_num += 1

                # Loops through libraries and outputs sample name (adding a not_used if no sample name or plate was loaded)
                for libn in lib_dict:
                    if plate in lib_dict[libn] and current_well in lib_dict[libn][plate]:
                        f.write(f"\t{lib_dict[libn][plate][current_well]}")
                    else:
                        f.write(f"\tnot_used{str(lib_count_dict['not_used'])}")
                        lib_count_dict['not_used'] += 1
                
                # Writes not_used for unused libraries
                if args.notlibs:
                    for _ in range(len(args.notlibs)):
                        f.write(f"\tnot_used{str(lib_count_dict['not_used'])}")
                        lib_count_dict['not_used'] += 1

                f.write("\n")
'''

