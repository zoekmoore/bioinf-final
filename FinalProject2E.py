# This program reads, stores, and determines
# the secondary structure of tRNA sequences.
# @authors Andrew Flynn and Zoe Moore
# @lastModified December 9, 2019

#################################################################################

def main():
    seq_list = readSeq() # populate our sequence list via "FinalProject.txt"
    for i in seq_list: # for each sequence in our sequence list...
        i = i.replace('T', 'U') # replace all instances of 'T' with 'U'
        find_secondary_structure(i) # find its secondary structure
   
#################################################################################
 
'''
readSeq - a function that reads in the sequences from a FASTA-formatted file and stores them in a list
'''
def readSeq():
    # File that contains the FASTA-formatted sequences
    database = "FinalProject.txt"

    # Declare variables
    seq_count = 0
    seq_list = []
    current_seq = ""

    # Iterate through each line of the file and store the sequences
    with open(database) as fp:
        line = fp.readline()
        while line:
            # If the line is a sequence name...
            if line[0] == '>':
                # If it's the first sequence, read the next line
                if seq_count == 0:
                    seq_count += 1
                    line = fp.readline()
                # If not, store our current sequence, "reset" the current
                # sequence, increment the sequence count, and read the next line
                else:
                    seq_list.append(current_seq)
                    current_seq = ""
                    seq_count += 1
                    line = fp.readline()
            else:
                current_seq += line
                line = fp.readline()
        # Change all Ts in the sequences to Us
        seq_list.append("".join(current_seq))	
    return seq_list

'''
find_secondary_structure - a function that follows biological conventions
to predict the secondary structure of tRNA sequences
Parameter: a given sequence that we are trying to predict the secondary
structure for
'''
def find_secondary_structure(sequence):
    # The following 8 variables are our region parameters
    a_stem_length = 7
    t_stem_length = 5
    t_loop_length = [5, 7]
    v_loop_length = [3, 24]
    d_stem_length = 4
    d_loop_length = [6, 11]
    anticodon_stem_length = 5
    anticodon_loop_length = 7

    # The minimum total length of a tRNA secondary structure
    total_length = a_stem_length * 2 + d_stem_length * 2 + anticodon_stem_length * 2 + t_stem_length * 2 + t_loop_length[0] + v_loop_length[0] + anticodon_loop_length + d_loop_length[0]

    # Remove all instances of '\n' from our sequence
    while sequence[-1] == '\n':
        sequence = sequence[:-1]

    structure = True                # Assume it has a secondary structure until it doesn't

#################################### A STEM ####################################

    # Declare / assign our variables
    a_stem = []
    i = 0
    j = len(sequence) - 1
    reset_count = 1
    allowed_mismatches = 0
    found_stem = False
    previous_mismatch = False

    # While we have less than 3 mismatches and haven't found the stem...
    while(allowed_mismatches < 3 and found_stem == False):
        # Reset our current mismatches
        current_mismatches = 0
        # While i and j are within parameter bounds...
        while(i < a_stem_length and j > total_length):
            # If they match...
            if match(sequence[i], sequence[j]):
                pair = [i, j, 0]                # No mismatch = 0
                # Append the matched pair, move to the next pair, and reset previous_mismatch to False
                a_stem.append(pair)             
                i += 1
                j -= 1
                previous_mismatch = False
                # If i = 6, we've found the a_stem
                if i == a_stem_length - 1:
                    found_stem = True
            # If they don't match, we haven't used our allowed mismatches yet, and previous_mismatch is False...
            elif current_mismatches < allowed_mismatches and previous_mismatch == False:
                # Append the mismatched pair, move to the next pair, increment our current mismatches, and set previous_mismatch to True
                pair = [i, j, 1]                # Mismatch = 1
                a_stem.append(pair)
                i += 1
                j -= 1
                current_mismatches += 1
                previous_mismatch = True
                # If i = 6, we've found the a_stem
                if i == a_stem_length - 1:
                    found_stem = True
            # If we encountered too many mismatches, or two consecutive mismatches, reset our search with the next lowest value for j 
            # This includes incrementing our reset_count, resetting previous_mismatch to False, and resetting our current mismatches to 0
            else:
                i = 0
                j = len(sequence) - 1 - reset_count
                a_stem = []
                reset_count += 1
                previous_mismatch = False
                current_mismatches = 0
        # If we've exhausted our search and failed to find an acceptable stem, COMPLETELY reset our search with an additional allowed mismatch
        if found_stem == False:
            j = len(sequence) - 1
            allowed_mismatches += 1
            reset_count = 1

############################### T STEM & T LOOP ################################

    # Declare / assign our variables
    t_stem = []
    temp_i = j
    i = temp_i
    temp_j = i - 2 * t_stem_length - t_loop_length[1]
    j = temp_j
    reset_count = 1
    allowed_mismatches = 0
    found_stem = False
    previous_mismatch = False

    # While we have less than 3 mismatches and haven't found the stem...
    while(allowed_mismatches < 3 and found_stem == False):
        # Reset our current mismatches
        current_mismatches = 0
        # While t_stem isn't fully populated...
        while(len(t_stem) < 5):
            # If they match...
            if match(sequence[i], sequence[j]):
                # Append the matched pair, move to the next pair, and reset previous_mismatch to False
                pair = [j, i, 0]               # j, i because swapped for indices; no mismatch = 0
                t_stem.append(pair)
                i -= 1
                j += 1
                previous_mismatch = False
                # We've found the t_stem if its length is 5
                if len(t_stem) == 5:
                    found_stem = True
            # If they don't match, we haven't used our allowed mismatches yet, and previous_mismatch is False...
            elif current_mismatches < allowed_mismatches and previous_mismatch == False:
                # Append the mismatched pair, move to the next pair, increment our current mismatches, and set previous_mismatch to True
                pair = [j, i, 1]               # j, i because swapped for indices; mismatch = 1
                t_stem.append(pair)
                i -= 1
                j += 1
                current_mismatches += 1
                previous_mismatch = False
                # We've found the t_stem if its length is 5
                if len(t_stem) == 5:
                    found_stem = True
            # If we encountered too many mismatches, or two consecutive mismatches, reset our search with the next highest value for j 
            # This includes incrementing our reset_count, resetting previous_mismatch to False, and resetting our current mismatches to 0
            else:
                i = temp_i
                j = temp_j + reset_count
                t_stem = []
                reset_count += 1
                previous_mismatch = False
                current_mismatches = 0
        # If we've exhausted our search and failed to find an acceptable stem, COMPLETELY reset our search with an additional allowed mismatch
        if found_stem == False:
            j = temp_j
            allowed_mismatches += 1
            reset_count = 1

    # For the t_loop...
    j = t_stem[-1][0] + 1 # j = next index after t_stem's last index
    t_loop = []
    v_max = temp_j
    # All nucleotides between t_stem belong to the t_loop
    while (j <= i):
        t_loop.append(j)
        j += 1

############################### D STEM & D LOOP ################################

    # Declare / assign our variables
    d_stem = []
    temp_i = a_stem_length + 2
    i = temp_i
    max_j = i + 2 * d_stem_length + d_loop_length[1]
    temp_j = i + 2 * d_stem_length + d_loop_length[0]
    j = temp_j
    reset_count = 1
    allowed_mismatches = 0
    found_stem = False
    previous_mismatch = False

    # While we have less than 3 mismatches and haven't found the stem...
    while(allowed_mismatches < 3 and found_stem == False):
        # Reset our current mismatches
        current_mismatches = 0
        # While d_stem isn't fully populated and j < max_j...
        while(len(d_stem) < 4 and j <= max_j):
            # If they match...
            if match(sequence[i], sequence[j]):
                # Append the matched pair, move to the next pair, and reset previous_mismatch to False        
                pair = [i, j, 0]                # No mismatch = 0
                d_stem.append(pair)
                i += 1
                j -= 1
                previous_mismatch = False
                # We've found the d_stem if its length is 4
                if len(d_stem) == 4:
                    found_stem = True
            # If they don't match, we haven't used our allowed mismatches yet, and previous_mismatch is False...
            elif current_mismatches < allowed_mismatches and previous_mismatch == False:
                # Append the mismatched pair, move to the next pair, increment our current mismatches, and set previous_mismatch to True
                pair = [i, j, 1]                # Mismatch = 1
                d_stem.append(pair)
                i += 1
                j -= 1
                current_mismatches += 1
                previous_mismatch = True
                # We've found the d_stem if its length is 4
                if len(d_stem) == 4:
                    found_stem = True
            # If we encountered too many mismatches, or two consecutive mismatches, reset our search with the next highest value for j 
            # This includes incrementing our reset_count, resetting previous_mismatch to False, and resetting our current mismatches to 0
            else:
                i = temp_i
                j = temp_j + reset_count
                d_stem = []
                reset_count += 1
                previous_mismatch = False
                current_mismatches = 0
        # If we've exhausted our search and failed to find an acceptable stem, COMPLETELY reset our search with an additional allowed mismatch
        if found_stem == False:
            j = temp_j
            allowed_mismatches += 1
            reset_count = 1

    # For the d_loop...
    j = d_stem[-1][1] # j = first index after d_loop
    d_loop = []
    # All nucleotides between d_stem belong to the d_loop
    while (i < j):
        d_loop.append(i)
        i += 1

####################### ANTICODON STEM & ANTICODON LOOP ########################

    # Declare / assign our variables
    anticodon_stem = []
    temp_i = a_stem_length + 2 + d_stem_length * 2 + len(d_loop) + 1
    i = temp_i
    temp_j = i + 2 * anticodon_stem_length + anticodon_loop_length
    j = temp_j
    reset_count = 1
    allowed_mismatches = 0
    found_stem = 0
    previous_mismatch = False

    # While we have less than 3 mismatches and haven't found the stem...
    while(allowed_mismatches < 3 and found_stem == False):
        # Reset our current mismatches
        current_mismatches = 0
        # While d_stem isn't fully populated...
        while(len(anticodon_stem) < 5 and j >= (i + anticodon_loop_length)):
            # If they match...
            if match(sequence[i], sequence[j]):
                # Append the matched pair, move to the next pair, and reset previous_mismatch to False   
                pair = [i, j, 0]                # No mismatch = 0
                anticodon_stem.append(pair)
                i += 1
                j -= 1
                previous_mismatch = False
                # We've found the anticodon_stem if its length is 5
                if len(anticodon_stem) == 5:
                    found_stem = True
            # If they don't match, we haven't used our allowed mismatches yet, and previous_mismatch is False...
            elif current_mismatches < allowed_mismatches and previous_mismatch == False:
                # Append the mismatched pair, move to the next pair, increment our current mismatches, and set previous_mismatch to True
                pair = [i, j, 1]                # No mismatch = 0
                anticodon_stem.append(pair)
                i += 1
                j -= 1
                current_mismatches += 1
                previous_mismatch = True
                # We've found the anticodon_stem if its length is 5
                if len(anticodon_stem) == 5:
                    found_stem = True
            # If we encountered too many mismatches, or two consecutive mismatches, reset our search with the next highest value for j 
            # This includes incrementing our reset_count, resetting previous_mismatch to False, and resetting our current mismatches to 0
            else:
                i = temp_i
                j = temp_j - reset_count
                anticodon_stem = []
                reset_count += 1
                previous_mismatch = False
                current_mismatches = 0
        # If we've exhausted our search and failed to find an acceptable stem, COMPLETELY reset our search with an additional allowed mismatch
        if found_stem == False:
            j = temp_j
            allowed_mismatches += 1
            reset_count = 1

    # For the anti_codon loop...
    anticodon_loop = []
    while (i <= j):
        anticodon_loop.append(i)
        i += 1

#################################### V LOOP ####################################

    # Declare / assign our variables
    v_loop = []
    j += anticodon_stem_length + 1 # j = first index after anticodon_stem
    while (j <= v_max):
        v_loop.append(j)
        j += 1

    # Call a function to output results
    if(structure == True):
        outputResults(sequence, a_stem, t_stem, d_stem, anticodon_stem, t_loop, d_loop, anticodon_loop, v_loop)
    else:
        print("NO SECONDARY STRUCTURE")     # FIX THIS SO IT SAYS SOMETHING BETTER

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
    temp = 3
    for i in range(indexCount, indexCount + len(d_stem)):
        nucleo5 += sequence[i]
        if(d_stem[temp][2] == 0):
            style5 += ')'
        else:
            style5 += '.'
        indexCount += 1
        temp -= 1

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
    temp = 4
    for i in range(indexCount, indexCount + len(anticodon_stem)):
        nucleo9 += sequence[i]
        if(anticodon_stem[temp][2] == 0):
            style9 += ')'
        else:
            style9 += '.'
        indexCount += 1
        temp -= 1
    
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
    for i in range(indexCount, indexCount + len(t_loop)):       # IN ORDER FOR ANDREW'S VERSION TO CORRECTLY PRINT COLORS, ADDED A -1 TO THIS   
        nucleo12 += sequence[i]
        style12 += '.'
        indexCount += 1

    # T Stem (right)
    nucleo13 = ''
    style13 = ''
    temp = 4
    for i in range(indexCount, indexCount + len(t_stem)):
        nucleo13 += sequence[i]
        if(t_stem[temp][2] == 0):
            style13 += ')'
        else:
            style13 += '.'
        indexCount += 1
        temp -= 1

    # A Stem (right)
    nucleo14 = ''
    style14 = ''
    temp = 6
    for i in range(indexCount, indexCount + len(a_stem)):
        nucleo14 += sequence[i]
        if(a_stem[temp][2] == 0):
            style14 += ')'
        else:
            style14 += '.'
        indexCount += 1
        temp -= 1

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
    whitespace = " 1"
    for i in range(1, len(sequence) - 1):
        whitespace += " "
    print(whitespace + str(len(sequence)))
    print(" \033[1;35;40m" + nucleo1 + "\033[0m" +
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
    print(" \033[1;35;40m" + style1 + "\033[0m" +
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

