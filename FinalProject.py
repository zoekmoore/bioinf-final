# This program reads, stores, and determines
# the secondary structure of tRNA sequences.
# @authors Andrew Flynn and Zoe Moore
# @lastModified November 20, 2019

#################################################################################

# This is the main function
def main():
    gene_list = read_genes()
    results = []
    print(gene_list)
    for i in gene_list:
        results.append(find_secondary_structure("GGGUCGUUAGCUCAGUUGGUAGAGCAAUUGACUUUUAAUCAAUUGGUCGCAGGUUCGAAUCCUGCACGACCCACCA"))
    #write_alignments(results)
   
#################################################################################
 
# This reads genes from "FinalProject.txt" and stores them
def read_genes():
    # File that contains FASTA formatted sequences
    database = "FinalProject.txt"

    # Declare variables
    gene_count = 0  # NOTE FOR ANDREW: is this necessary?
    gene_list = []
    current_gene = []

    # Iterate through each line of the file and store the sequences
    with open(database) as fp:
        line = fp.readline()
        while line:
            if line[0] == '>':
                if gene_count == 0:
                    gene_count += 1
                    line = fp.readline()
                else:
                    gene_list.append(current_gene)
                    current_gene = ""
                    gene_count += 1
                    line = fp.readline()
            else:
                current_gene += line
                line = fp.readline()
        # Change all Ts in the sequences to Us
        for i in range(len(current_gene)):
            if current_gene[i] == 'T':
                current_gene[i] = 'U'
        gene_list.append("".join(current_gene))
    return gene_list

# This function finds the secondary structure of an RNA sequence
def find_secondary_structure(gene):
    # The following 8 variables are our region parameters
    a_stem_length = 7
    d_stem_length = [3, 4]
    anticodon_stem_length = 5
    t_stem = 5
    t_loop_length = [5, 7]
    v_loop_length = [3, 24]
    anticodon_loop_length = 7
    d_loop_length = [4, 11]

#####################################A STEM#####################################

    a_stem = []
    i = 0
    j = len(gene) - 1
    reset_count = 1

    while(i < 7):
        if match(gene[i], gene[j]):
            pair = [i, j]
            a_stem.append(pair)
            i += 1
            j -= 1
        else:
            i = 0
            j = len(gene) - 1 - reset_count
            a_stem = []
            reset_count += 1
    print("A STEM: " + str(a_stem))

################################T STEM & T LOOP#################################

    t_stem = []
    temp_i = j
    i = temp_i
    temp_j = i - 17
    j = temp_j
    reset_count = 1
    while(len(t_stem) < 5):
        if match(gene[i], gene[j]):
            pair = [j, i]               # j, i because swapped for indices
            t_stem.append(pair)
            i -= 1
            j += 1
        else:
            i = temp_i
            j = temp_j + reset_count
            t_stem = []
            reset_count += 1
    t_loop = []
    v_max = j - 6
    while (j <= i):
        t_loop.append(j)
        j += 1
    print("T STEM: " + str(t_stem))
    print("T LOOP: " + str(t_loop))

################################D STEM & D LOOP#################################

    d_stem = []
    temp_i = 9
    i = temp_i
    temp_j = i + 19
    j = temp_j
    reset_count = 1
    while(len(d_stem) < 4):
        if match(gene[i], gene[j]):
            pair = [i, j]
            d_stem.append(pair)
            i += 1
            j -= 1
        else:
            if (len(d_stem) == 3):
                break
            i = temp_i
            j = temp_j - reset_count
            d_stem = []
            reset_count += 1
    d_loop = []
    while (i <= j):
        d_loop.append(i)
        i += 1
    print("D STEM: " + str(d_stem))
    print("D LOOP: " + str(d_loop))

########################ANTICODON STEM & ANTICODON LOOP#########################

    anticodon_stem = []
    temp_i = len(a_stem) + 2 + len(d_stem) * 2 + len(d_loop) + 1
    i = temp_i
    temp_j = i + 17
    j = temp_j
    reset_count = 1
    while(len(anticodon_stem) < 5):
        if match(gene[i], gene[j]):
            pair = [i, j]
            anticodon_stem.append(pair)
            i += 1
            j -= 1
        else:
            i = temp_i
            j = temp_j - reset_count
            anticodon_stem = []
            reset_count += 1
    anticodon_loop = []
    while (i <= j):
        anticodon_loop.append(i)
        i += 1
    print("ANTICODON STEM: " + str(anticodon_stem))
    print("ANTICODON LOOP: " + str(anticodon_loop))

#####################################V LOOP#####################################

    v_loop = []
    j += 6
    while (j <= v_max):
        v_loop.append(j)
        j += 1
    print("V LOOP: " + str(v_loop))

    # Call a function to reformat the structure
    result = reformat(a_stem, t_stem, d_stem, anticodon_stem, t_loop, d_loop, anticodon_loop, v_loop)

    return t_stem # we'll just return this for now

#################################################################################

# this is a "lazy" match function which we'll improve upon
def match(i, j): # we're temporarily working with sequences with T's instead of U's --> we can easily change this later
    if i == 'A':
        if j == 'U':
            return True
    if i == 'U':
        if j == 'A':
            return True
    if i == 'C':
        if j == 'G':
            return True
    if i == 'G':
        if j == 'C':
            return True
    return False # if not a match, return false

#################################################################################

def reformat(a_stem, t_stem, d_stem, anticodon_stem, t_loop, d_loop, anticodon_loop, v_loop):
    

#################################################################################

if __name__ == '__main__':
    main()
