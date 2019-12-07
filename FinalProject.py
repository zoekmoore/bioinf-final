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
    # TODO: the output for the stems is of the form [low index, high index, 0/1 for no mismatch/yes mismatch]
    # but we need to update where we assign these values for the input in each of these segments in order to
    # properly have the third 0/1 business

    a_stem = []
    i = 0
    j = len(gene) - 1
    reset_count = 1

    while(i < 7):
        if match(gene[i], gene[j]):
            pair = [i, j, 0]                # No mismatch = 0
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
            pair = [j, i, 0]               # j, i because swapped for indices; no mismatch = 0
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
            pair = [i, j, 0]                # No mismatch = 0
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
            pair = [i, j, 0]                # No mismatch = 0
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
    result = reformat(gene, a_stem, t_stem, d_stem, anticodon_stem, t_loop, d_loop, anticodon_loop, v_loop)

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

def reformat(sequence, a_stem, t_stem, d_stem, anticodon_stem, t_loop, d_loop, anticodon_loop, v_loop):
    
    indexCount = 0

    # A Stem (left)
    nucleo1 = ''
    nucleo2 = ''
    style1 = ''
    style2 = ''

    for i in range(0, 7):
        nucleo1 += sequence[i]
        if(a_stem[i][2] == 0):
            style1 += '('
        else:
            style1 += '.'
        indexCount += 1
    nucleo2 += sequence[7]
    style2 += '.'
    nucleo2 += sequence[8]
    style2 += '.'
    indexCount += 2

    # D Stem (left)
    nucleo3 = ''
    style3 = ''
    for i in range(indexCount, indexCount + len(d_stem)):
        nucleo3 += sequence[i]
        if(d_stem[i - indexCount][2] == 0):
            style3 += '('
        else:
            style3 += '.'
        indexCount += 1
    
    # D Loop
    nucleo4 = ''
    style4 = ''
    for i in range(indexCount, indexCount + len(d_loop)):
        nucleo4 += sequence[i]
        style4 += '.'
        indexCount += 1

    # D Stem (right)
    nucleo5 = ''
    style5 = ''
    for i in range(indexCount, indexCount + len(d_stem)):
        nucleo5 += sequence[i]
        if(d_stem[i - indexCount][2] == 0):
            style5 += ')'
        else:
            style5 += '.'
        indexCount += 1

    # Gap between D Stem (right) and Anticodon Stem (left)
    nucleo6 = ''
    style6 = ''
    for i in range(indexCount, anticodon_stem[0][0]):
        nucleo6 += sequence[i]
        style6 += '.'
        indexCount += 1

    # Anticodon Stem (left)
    nucleo7 = ''
    style7 = ''
    for i in range(indexCount, indexCount + len(anticodon_stem)):
        nucleo7 += sequence[i]
        if(anticodon_stem[i - indexCount][2] == 0):
            style7 += '('
        else:
            styel7 += '.'
        indexCount += 1

    # Anticodon Loop
    nucleo8 = ''
    style8 = ''
    for i in range(indexCount, indexCount + len(anticodon_loop)):
        nucleo8 += sequence[i]
        style8 += '.'
        indexCount += 1

    # Anticodon Stem (right)
    nucleo9 = ''
    style9 = ''
    for i in range(indexCount, indexCount + len(anticodon_stem)):
        nucleo9 += sequence[i]
        if(anticodon_stem[i - indexCount][2] == 0):
            style9 += ')'
        else:
            style9 += '.'
        indexCount += 1
    
    # V Loop
    nucleo10 = ''
    style10 = ''
    for i in range(indexCount, indexCount + len(v_loop)):
        nucleo10 += sequence[i]
        style10 += '.'
        indexCount += 1

    # T Stem (left)
    nucleo11 = ''
    style11 = ''
    for i in range(indexCount, indexCount + len(t_stem)):
        nucleo11 += sequence[i]
        if(t_stem[i - indexCount][2] == 0):
            style11 += '('
        else:
            style11 += '.'
        indexCount += 1

    # T Loop
    nucleo12 = ''
    style12 = ''
    for i in range(indexCount, indexCount + len(t_loop)):
        nucleo12 += sequence[i]
        style12 += '.'
        indexCount += 1

    # T Stem (right)
    nucleo13 = ''
    style13 = ''
    for i in range(indexCount, indexCount + len(t_stem)):
        nucleo13 += sequence[i]
        if(t_stem[i - indexCount][2] == 0):
            style13 += ')'
        else:
            style13 += '.'
        indexCount += 1

    # A Stem (right)
    # T Stem (right)
    nucleo14 = ''
    style14 = ''
    for i in range(indexCount, indexCount + len(a_stem)):
        nucleo14 += sequence[i]
        if(a_stem[i - indexCount][2] == 0):
            style14 += ')'
        else:
            style14 += '.'
        indexCount += 1

    # Remaining Nucleotides
    nucleo15 = ''
    style15 = ''
    for i in range(indexCount, len(sequence)):
        nucleo15 += sequence[i]
        style15 += '.'
        indexCount += 1

    # Print Key
    print("\033[1;37;40m KEY: \033[0m")
    print("\t\033[1;35;40m" + "A Stem" + "\033[0m")
    print("\t\033[1;32;40m" + "D Stem" + "\033[0m")
    print("\t\033[1;31;40m" + "Anticodon Stem" + "\033[0m")
    print("\t\033[1;34;40m" + "T Stem" + "\033[0m")
    print("\tLoops & Etc.\n\n")
    # Print color-coded nucleotides
    print("\033[1;37;40m SEQUENCE: \033[0m")
    print("  \033[1;35;40m" + nucleo1 + "\033[0m" +
    nucleo2 + 
    "\033[1;32;40m" + nucleo3 + "\033[0m"  + nucleo4 + "\033[1;32;40m" + nucleo5 + "\033[0m" +
    nucleo6 +
    "\033[1;31;40m" + nucleo7 + "\033[0m" + nucleo8 + "\033[1;31;40m" + nucleo9 + "\033[0m" +
    nucleo10 +
    "\033[1;34;40m" + nucleo11 + "\033[0m" + nucleo12 + "\033[1;34;40m" + nucleo13 + "\033[0m" +
    "\033[1;35;40m" + nucleo14 + "\033[0m" +
    nucleo15)
    # Print color-coded match/mismatch stylization
    print("\033[1;37;40m FOLDING: \033[0m")
    print("  \033[1;35;40m" + style1 + "\033[0m" +
    style2 + 
    "\033[1;32;40m" + style3 + "\033[0m"  + style4 + "\033[1;32;40m" + style5 + "\033[0m" +
    style6 +
    "\033[1;31;40m" + style7 + "\033[0m" + style8 + "\033[1;31;40m" + style9 + "\033[0m" +
    style10 +
    "\033[1;34;40m" + style11 + "\033[0m" + style12 + "\033[1;34;40m" + style13 + "\033[0m" +
    "\033[1;35;40m" + style14 + "\033[0m" +
    style15)
    
#################################################################################

if __name__ == '__main__':
    main()
