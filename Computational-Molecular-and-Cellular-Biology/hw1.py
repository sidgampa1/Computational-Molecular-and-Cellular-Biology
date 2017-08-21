import sys

##helper method
##takes a sequence and creates reverse complement
def complement(seq, name, rna):
    complement_table = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if (rna):
        complement_table["A"] = "U"
    new_seq = ""
    for letter in seq:
        complement_letter = complement_table[letter]
        new_seq += complement_letter
    new_seq = new_seq[::-1]
    length = len(new_seq)

    ##print sequence in FASTA format
    print(name)
    print_seq(new_seq, 0, 70, length)
    print

##helper method
##print sequence; start inclusive, end exclusive
def print_seq(new_seq, start, end, length):
    if end > length:
        print(new_seq[start: length])
    elif (start < length and end <= length):
        print(new_seq[start: end])
        print_seq(new_seq, start + 70, end + 70, length)


rna = False ##reverse complement in rna bases or dna bases

##main block
##open file
try:
    file = open(sys.argv[1], 'r')
    if (len(sys.argv) == 3 and sys.argv[2] == "rna"):
        rna = True
except IOError:
    print("File cannot be opened.")
except IndexError:
    print("You must provide a file name!")

line = file.readline()
sequence = ""

##read file to create reverse complement
while (line != ''):
    if (line[0] == '>'):
        if (sequence != ""):
            complement(sequence, name, rna)
            sequence = ""
        name = line[1:len(line) - 1]
    else:
        sequence += line[:len(line) - 1]
    line = file.readline()

complement(sequence, name, rna)
