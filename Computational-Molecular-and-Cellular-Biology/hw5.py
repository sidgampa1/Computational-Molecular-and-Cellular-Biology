import sys
import os.path
import itertools

## check if corect number of args are given
if len(sys.argv) != 3:
    print("Incorrect number of arguments, try again")
    sys.exit()

## block to open and read scoring and sequence files
try:
    file_matrix = open(sys.argv[1], 'r')
    scoring_matrix = {} ## dict. of letter pairs: score
    row_count = 0
    allowed_letters = {} ## dict. of letter_index: letter
    for line in file_matrix:
        row = line.replace(" ", "")
        row = row.rstrip('\n')
        if row_count == 0:
            char_count = 1
            for char in row:
                allowed_letters[char_count] = char
                char_count += 1

        elif row_count <= len(allowed_letters):
            col_count = 1
            negative = False
            score = 0
            for char in row:
                first_letter = allowed_letters[row_count]
                second_letter = allowed_letters[col_count]
                pair = first_letter + second_letter
                if char != '-' and not negative:
                    score = float(char)
                elif char != '-' and negative:
                    score = float(char) * -1
                    negative = False
                elif char == '-':
                    negative = True

                if not negative:
                    scoring_matrix[pair] = score
                    col_count += 1

        else:
            scoring_matrix['GAP'] = float(row)

        row_count += 1

    file_matrix.close()

## read in FASTA file sequences
    sequence_file = open(sys.argv[2], 'r')
    sequences = []
    sequence = ""
    for line in sequence_file:
        if (line[0] == '>'):
            if (sequence != ""):
                sequences.append(sequence)
                sequence = ""
        else:
            sequence += line.rstrip('\n')

    sequences.append(sequence)
    sequence_file.close()

except IOError:
    sys.stderr.write("Incorrect type of file, please try again")
    sys.exit()

## fill in dynamic programming matrix

dyn_matrix = [[0 for x in range(len(sequences[0]) + 1)] for y in range(len(sequences[1]) + 1)]

score = 0
gap = scoring_matrix['GAP']
count = 0
for row in dyn_matrix:
    if (count == 0):
        for i in range(len(row)):
            row[i] = score
            score += gap
        score = gap
        count += 1
    else:
        row[0] = score
        score += gap

for row in range(1, len(dyn_matrix)): ## fill in dyn_matrix greedily
    for col in range(1, len(dyn_matrix[0])):
        letter1 = sequences[0][col - 1]
        letter2 = sequences[1][row - 1]
        pair = letter1 + letter2
        matching = dyn_matrix[row - 1][col - 1] + scoring_matrix[pair]
        upper_gap = dyn_matrix[row - 1][col] + gap
        left_gap = dyn_matrix[row][col - 1] + gap
        dyn_matrix[row][col] = max(matching, upper_gap, left_gap)

row = len(dyn_matrix) - 1
col = len(dyn_matrix[0]) - 1
seq1 = ""
seq2 = ""
while row != 0 and col != 0: ## find best alignment backwards
    letter1 = sequences[0][col - 1]
    letter2 = sequences[1][row - 1]
    pair = letter1 + letter2
    matching = dyn_matrix[row - 1][col - 1] + scoring_matrix[pair]
    left_gap = dyn_matrix[row][col - 1] + gap
    upper_gap = dyn_matrix[row - 1][col] + gap
    current = dyn_matrix[row][col]
    if matching == current:
        row -= 1
        col -= 1
        seq1 += sequences[0][col]
        seq2 += sequences[1][row]
    elif left_gap == current:
        col -= 1
        seq1 += sequences[0][col]
        seq2 += '-'
    elif upper_gap == current:
        row -= 1
        seq1 += '-'
        seq2 += sequences[1][row]

## print dynamic matrix
for row in dyn_matrix:
    string = ""
    for element in row:
        string += str(int(element))
        string += " "
    print(string)

## print seq alignments as FASTA output
print(">AlignmentM")
print(seq1[::-1])
print(">AlignmentN")
print(seq2[::-1])
