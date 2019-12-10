# This program reads, stores, and determines
# the secondary structure of tRNA sequences.
# @authors Andrew Flynn and Zoe Moore
# @lastModified December 8, 2019

#################################################################################

def main():
    gene_list = readSeq()
    results = []
    for i in gene_list:
        i = i.replace('T', 'U')
        results.append(find_secondary_structure(i))
    #write_alignments(results)
   
#################################################################################
 
'''
readSeq - a function that reads in the sequences from a FASTA formatted file
and stores them in a list
'''
def readSeq():
    # File that contains FASTA formatted sequences
    database = "FinalProject.txt"

    # Declare variables
    gene_count = 0  # NOTE FOR ANDREW: is this necessary?
    gene_list = []
    current_gene = ""

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
        gene_list.append("".join(current_gene))	
    return gene_list

'''
find_secondary_structure - a function that follows biological conventions
to predict the secondary structure of tRNA sequences
Parameter: a given sequence that we are trying to predict the secondary
structure for
'''
def find_secondary_structure(sequence):
    # The following 8 variables are our region parameters
    a_stem_length = 7
    d_stem_length = [3, 4]
    anticodon_stem_length = 5
    t_stem = 5
    t_loop_length = [5, 7]
    v_loop_length = [3, 24]
    anticodon_loop_length = 7
    d_loop_length = [4, 11]

    total_length = a_stem_length * 2 + d_stem_length[0] * 2 + anticodon_stem_length * 2 + t_stem * 2 + t_loop_length[0] + v_loop_length[0] + anticodon_loop_length + d_loop_length[0]

    structure = True                # Assume it has a secondary structure until it doesn't

#################################### A STEM ####################################
   
 # TODO: the output for the stems is of the form [low index, high index, 0/1 for no mismatch/yes mismatch]
    # but we need to update where we assign these values for the input in each of these segments in order to
    # properly have the third 0/1 business
    while sequence[-1] == '\n':
        sequence = sequence[:-1]
    a_stem = []
    i = 0
    j = len(sequence) - 1
    reset_count = 1
    allowed_mismatches = 0
    found_stem = False
    previous_mismatch = False

    while(allowed_mismatches < 3 and found_stem == False):
        current_mismatches = 0
        while(i < 7 and j > total_length):
            if match(sequence[i], sequence[j]):
                pair = [i, j, 0]                # No mismatch = 0
                a_stem.append(pair)
                i += 1
                j -= 1
                previous_mismatch = False
                if i == 7:
                    found_stem = True
            elif current_mismatches < allowed_mismatches and previous_mismatch == False:
                pair = [i, j, 1]                # Mismatch = 1
                a_stem.append(pair)
                i += 1
                j -= 1
                current_mismatches += 1
                previous_mismatch = True
                if i == 7:
                    found_stem = True
            else:
                i = 0
                j = len(sequence) - 1 - reset_count
                a_stem = []
                reset_count += 1
                previous_mismatch = False
                current_mismatches = 0
        if found_stem == False:
            j = len(sequence) - 1
            allowed_mismatches += 1
            reset_count = 1

############################### T STEM & T LOOP ################################

    t_stem = []
    temp_i = j
    i = temp_i
    temp_j = i - 17
    j = temp_j
    reset_count = 1
    allowed_mismatches = 0
    found_stem = 0
    previous_mismatch = False

    while(allowed_mismatches < 3 and found_stem == False):
        current_mismatches = 0
        while(len(t_stem) < 5):
            if match(sequence[i], sequence[j]):
                pair = [j, i, 0]               # j, i because swapped for indices; no mismatch = 0
                t_stem.append(pair)
                i -= 1
                j += 1
                previous_mismatch = False
                if len(t_stem) == 5:
                    found_stem = True
            elif current_mismatches < allowed_mismatches and previous_mismatch == False:
                pair = [j, i, 1]               # j, i because swapped for indices; mismatch = 1
                t_stem.append(pair)
                i -= 1
                j += 1
                current_mismatches += 1
                previous_mismatch = False
                if len(t_stem) == 5:
                    found_stem = True
            else:
                previous_mismatch = False
                i = temp_i
                j = temp_j + reset_count
                t_stem = []
                reset_count += 1
        if found_stem == False:
            j = temp_j
            allowed_mismatches += 1

    j = t_stem[-1][0]
    t_loop = []
    v_max = j - 5
    while (j <= i):
        t_loop.append(j)
        j += 1

############################### D STEM & D LOOP ################################

    d_stem = []
    temp_i = 9
    i = temp_i
    temp_j = i + 19
    j = temp_j
    reset_count = 1
    allowed_mismatches = 0
    found_stem = 0
    previous_mismatch = False

    while(allowed_mismatches < 3 and found_stem == False):
        current_mismatches = 0
        while(len(d_stem) < 4):
            if match(sequence[i], sequence[j]):
                pair = [i, j, 0]                # No mismatch = 0
                d_stem.append(pair)
                i += 1
                j -= 1
                previous_mismatch = False
                if len(d_stem) == 4:
                    found_stem = True
            elif current_mismatches < allowed_mismatches and previous_mismatch == False and len(d_stem) != 3:
                pair = [i, j, 1]                # Mismatch = 1
                d_stem.append(pair)
                i += 1
                j -= 1
                current_mismatches += 1
                previous_mismatch = True
                if len(d_stem) == 4:
                    found_stem = True
            else:
                previous_mismatch = False
                if (len(d_stem) == 3):
                    break
                i = temp_i
                j = temp_j - reset_count
                d_stem = []
                reset_count += 1
        if found_stem == False:
            j = temp_j
            allowed_mismatches += 1
    j = d_stem[-1][1] - 1
    d_loop = []
    while (i <= j):
        d_loop.append(i)
        i += 1

####################### ANTICODON STEM & ANTICODON LOOP ########################

    anticodon_stem = []
    temp_i = len(a_stem) + 2 + len(d_stem) * 2 + len(d_loop) + 1
    i = temp_i
    temp_j = i + 17
    j = temp_j
    reset_count = 1
    allowed_mismatches = 0
    found_stem = 0
    previous_mismatch = False

    while(allowed_mismatches < 3 and found_stem == False):
        current_mismatches = 0
        while(len(anticodon_stem) < 5):
            if match(sequence[i], sequence[j]):
                pair = [i, j, 0]                # No mismatch = 0
                anticodon_stem.append(pair)
                i += 1
                j -= 1
                previous_mismatch = False
                if len(anticodon_stem) == 5:
                    found_stem = True
            elif current_mismatches < allowed_mismatches and previous_mismatch == False:
                pair = [i, j, 1]                # No mismatch = 0
                anticodon_stem.append(pair)
                i += 1
                j -= 1
                current_mismatches += 1
                previous_mismatch = True
                if len(anticodon_stem) == 5:
                    found_stem = True
            else:
                previous_mismatch = False
                i = temp_i
                j = temp_j - reset_count
                anticodon_stem = []
                reset_count += 1
        if found_stem == False:
            j = temp_j
            allowed_mismatches += 1

    anticodon_loop = []
    while (i <= j):
        anticodon_loop.append(i)
        i += 1

#################################### V LOOP ####################################

    v_loop = []
    j += 6
    while (j <= v_max):
        v_loop.append(j)
        j += 1

    # Call a function to output results
    if(structure == True):
        outputResults(sequence, a_stem, t_stem, d_stem, anticodon_stem, t_loop, d_loop, anticodon_loop, v_loop)
    else:
        print("NO SECONDARY STRUCTURE")     # FIX THIS SO IT SAYS SOMETHING BETTER

    return t_stem # we'll just return this for now

#################################################################################

'''
match - a function that checks two nucleotides against one another to determine
whether or not they are a match
Parameter: first nucleotide
Parameter: second nucleotide
'''
def match(i, j):
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
    if i == 'U':
        if j == 'G':
            return True
    if i == 'G':
        if j == 'U':
            return True
    return False # if not a match, return false

#################################################################################

'''
outputResults - a function that, if a valid secondary structure is found,
takes the results of the pairing and outputs them in a color-coded, labeled
fashion
Parameter: the sequence that we are trying to predict the secondary structure
Parameters: the predicted stems, in the form [first index, second index, (0/1)match/mismatch]
Parameters: the predicted loops, listing the indices of the loop's nucleotides
'''
def outputResults(sequence, a_stem, t_stem, d_stem, anticodon_stem, t_loop, d_loop, anticodon_loop, v_loop):
    
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
    temp = indexCount
    for i in range(indexCount, indexCount + len(d_stem)):
        nucleo3 += sequence[i]
        if(d_stem[i - temp][2] == 0):
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
    temp = indexCount
    for i in range(indexCount, indexCount + len(d_stem)):
        nucleo5 += sequence[i]
        if(d_stem[i - temp][2] == 0):
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
    temp = indexCount
    for i in range(indexCount, indexCount + len(anticodon_stem)):
        nucleo7 += sequence[i]
        if(anticodon_stem[i - temp][2] == 0):
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
    temp = indexCount
    for i in range(indexCount, indexCount + len(anticodon_stem)):
        nucleo9 += sequence[i]
        if(anticodon_stem[i - temp][2] == 0):
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
    temp = indexCount
    for i in range(indexCount, indexCount + len(t_stem)):
        nucleo11 += sequence[i]
        if(t_stem[i - temp][2] == 0):
            style11 += '('
        else:
            style11 += '.'
        indexCount += 1

    # T Loop
    nucleo12 = ''
    style12 = ''
    for i in range(indexCount, indexCount + len(t_loop) - 1):       # IN ORDER FOR ANDREW'S VERSION TO CORRECTLY PRINT COLORS, ADDED A -1 TO THIS   
        nucleo12 += sequence[i]
        style12 += '.'
        indexCount += 1

    # T Stem (right)
    nucleo13 = ''
    style13 = ''
    temp = indexCount
    for i in range(indexCount, indexCount + len(t_stem)):
        nucleo13 += sequence[i]
        if(t_stem[i - temp][2] == 0):
            style13 += ')'
        else:
            style13 += '.'
        indexCount += 1

    # A Stem (right)
    nucleo14 = ''
    style14 = ''
    temp = indexCount
    for i in range(indexCount, indexCount + len(a_stem)):
        nucleo14 += sequence[i]
        if(a_stem[i - temp][2] == 0):
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
    print("\n\n\033[1;37;40m KEY: \033[0m")
    print("\t\033[1;35;40m" + "A Stem" + "\033[0m")
    print("\t\033[1;32;40m" + "D Stem" + "\033[0m")
    print("\t\033[1;31;40m" + "Anticodon Stem" + "\033[0m")
    print("\t\033[1;34;40m" + "T Stem" + "\033[0m")
    print("\tLoops & Etc.\n\n")
    # Print color-coded nucleotides
    print("\033[1;37;40m SEQUENCE: \033[0m")
    whitespace = "  1"
    for i in range(1, len(sequence) - 1):
        whitespace += " "
    print(whitespace + str(len(sequence)))
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
