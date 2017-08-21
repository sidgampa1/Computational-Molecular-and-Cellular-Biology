import sys
import os.path
import random

'''
***params file format example:***

name: gi|6599362|emb|AJ131352.1| Trema virgata gene encoding hemoglobin, isolate T5
length: 1104
A 0.29981884058
C 0.173007246377
T 0.329710144928
G 0.197463768116

***with motif = nonempty motif with no internal stop codons:***

name: gi|6599362|emb|AJ131352.1| Trema virgata gene encoding hemoglobin, isolate T5
length: 1104
motif: CATCATCATCATCAT
A 0.29981884058
C 0.173007246377
T 0.329710144928
G 0.197463768116
'''

args = sys.argv
args_length = len(args)

## booleans to parse command line args
sim = (args_length == 1 and args[0] == "simulator.py")
save = (args_length == 5 and args[3] == "--save")
calc = (not sim) and (args[1] == "--calc") and ((args_length == 3) or save)
load = (not sim) and (args[1] == "--load") and ((args_length == 3) or save)
usage_help = (args_length == 2 and (args[1] == "-h" or args[1] == "--help"))
subseq_count = {} ## counts for subseqs if necessary
k = 0

## default nucleotide count and output seq name
    ## can be modified to fit commands
Nuc_Count = {"A": 249, "T": 249, "C": 249, "G": 249}
seq_name = "random seq|equal distribution of A, T, G, C"

## for use when searching for internal stop codons
stop_codons = {"TAA": "ATA", "TAG": "AGT", "TGA": "GTA"}

## helper method
## print sequence; start inclusive, end exclusive
def print_seq(new_seq, start, end, length):
    if end > length:
        print(new_seq[start: length])
    elif (start < length and end <= length):
        print(new_seq[start: end])
        print_seq(new_seq, start + 60, end + 60, length)


## helper method
## counts subsequence
def subsequence_count(seq):
    remainder = len(seq) % k
    if (remainder != 0):
        start = 0
        end = len(seq) - remainder
        while end < len(seq):
            subsequences = [seq[i:i+k] for i in range(start, end, k)]
            subsequence_count_helper(subsequences)
            start += 1
            end += 1
    else:
        subsequences = [seq[i:i+k] for i in range(0, len(seq), k)]
        subsequence_count_helper(subsequences)


## helper method
## adds subsequence count to dictionary
def subsequence_count_helper(subsequences):
    for subseq in subsequences:
        if (subseq in subseq_count):
            # print("counting subseq")
            subseq_count[subseq] += 1
        else:
            # print("adding new subseq")
            subseq_count[subseq] = 1


## uses Nuc_Count to create random seq with correct composition
def random_seq(motif = ""):
    seq = []
    for key, value in Nuc_Count.iteritems():
        while (value > 0):
            seq.append(key)
            value -= 1

    seq = randomize_and_clean(seq, motif)

    ## output in FASTA format
    print(">" + seq_name)
    print_seq(seq, 0, 60, len(seq))

    if (k != 0):
        subsequence_count(seq)
    ##save params if needed
    if (save):
        try:
            file = open(args[4], 'a')
            file.write("name: " + seq_name + '\n')
            file.write("length: " + str(len(seq)) + '\n')
            if (motif != ""):
                file.write("motif: " + motif + '\n')
            for key, value in Nuc_Count.iteritems():
                file.write(key + " " + str(float(value) / len(seq)) + '\n')
        except ValueError:
            print("wrong file or incorrect formatting")
        except IOError:
            print("Wrong file type, please try again")


def shuffle_stop_codons(seq):
    list_of_codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    seq = [stop_codons[codon] if codon in stop_codons else codon for codon in list_of_codons]
    seq = "".join(seq)
    return seq

def insert_motif(seq, motif):
    ## throw and error if motif is invalid
    for key in stop_codons.keys():
        if key in motif:
            print("Motif contains internal stop codon... it will not be inserted. Please fix and run the program again.")
            sys.exit()
    if (len(motif) > len(seq)):
        print("Motif is longer than sequence... it will not be inserted. Please fix and run the program again")
        sys.exit()

    ## randomly delete chars in motif from seq
    indices = list(enumerate(seq))
    A_inds = [i for i,letter in indices if letter == "A"]
    T_inds = [i for i,letter in indices if letter == "T"]
    G_inds = [i for i,letter in indices if letter == "G"]
    C_inds = [i for i,letter in indices if letter == "C"]
    dic = {"A": A_inds, "T": T_inds, "G": G_inds, "C": C_inds}
    lst = list(seq)
    for letter in motif:
        index = random.choice(dic[letter])
        lst[index] = ""
    new = "".join(lst)

    ## if length is still not short enough, remove some more nucs
    while (len(new) + len(motif) > len(seq)):
        new = new.replace(random.choice(dic.keys()), "", 1)

    tries = 0
    while (motif != ""):
        insert_point = random.choice(range(0, len(new)))
        new = new[:insert_point] + motif + new[insert_point:]

        list_of_codons = [new[i:i+3] for i in range(0, len(new), 3)]
        stops = [codon for codon in list_of_codons if codon in stop_codons]

        if (len(stops) == 0):
            return new
        elif (tries == 50): ## if re-inserting does not work, shuffle surrounding nucs
            shuffle_stop_codons(new)
            return new
        else:
            new = new[:insert_point] + new[insert_point + len(motif):]
            tries += 1
    return new


## helper function
## changes premature stop codons, shuffles sequence, fixes length
def randomize_and_clean(seq, motif = ""):
    remainder = len(seq) % 3
    ## fix length ~ must be multiple of 3 for codons
    if (remainder):
        for i in range(0, 3 - remainder):
            seq += random.choice(Nuc_Count.keys())
    ## shuffle sequence, change to String obj for cleanup
    random.shuffle(seq)
    seq = "".join(seq)
    ## pre-fix distribution by replacing to be deleted single codons [ATG, TAA]
    replaced_bases = {"A": [seq[0], seq[len(seq) - 2]], "T": [seq[1], seq[len(seq) - 3]], "G": [seq[2], seq[len(seq) - 1]]}
    for key, value in replaced_bases.iteritems():
        for base in value:
            seq = seq.replace(key, base, 1)
    ## search for and shuffle premature stop codons
    seq = shuffle_stop_codons(seq)
    ## insert motif (checks if necessary or not)
    seq = insert_motif(seq, motif)

    ## write in start and stop codons at beginning and end
    seq = "ATG" + seq[3:]
    seq = seq[:(len(seq) - 3)] + "TAG"
    return seq

## helper function
## counts nuc comp, saves to Nuc_Count
def count_nucs(seq):
    ## reset dic before use
    for key in Nuc_Count.keys():
        Nuc_Count[key] = 0

    ## count number of each nuc
    for nuc in seq:
        Nuc_Count[nuc] += 1

def load_file(params):
    if (not os.path.exists(params)):
        print("file '" + str(params) + "' does not exist, please try again")
        sys.exit()

    ##open file and read in data
    global seq_name
    length = 0
    motif = ""
    try:
        file = open(params, 'r')
        for line in file:
            tag = line.split()[0]
            if (tag == "name:"):
                seq_name = line[6:len(line) - 1]
            elif (tag == "length:"):
                length = int(line[8:len(line) - 1])
            elif (tag == "motif:"):
                motif = line[7: len(line) - 1]
            else:
                count = float(line[2:len(line) - 1]) * length
                Nuc_Count[tag] = round(count)
                if (tag == "G"):
                    random_seq(motif) ## generate seq at end of each seq data
                    motif = ""
    except ValueError:
        print("wrong file or incorrect formatting")
    except IOError:
        print("Wrong file type, please try again")



def calc_stats(seq_file):

    ## calc stats
    seq = ""
    global seq_name
    ##read file to get sequence
    try:
        file = open(seq_file, 'r')
        ## ask user if subseq count desired
        answer = raw_input("Would you like frequencies of subsequences? [Y/N]")
        if (answer == "y" or answer == "Y"):
            global k
            k = input("What length subsequence would you like to analyze?")
            if (not isinstance(k, int) or k < 0):
                k = input("k is not a valid positive integer, please enter a number again or enter a non number to continue")
            if (not isinstance(k, int) or k < 0):
                k = 0
        ## read file
        for line in file:
            if (line[0] == '>'):
                if (seq != ""):
                    count_nucs(seq)
                    random_seq()
                    seq = ""
                seq_name = line[1:len(line) - 1]
            else:
                seq += line[:len(line) - 1]
        count_nucs(seq)
        random_seq()
    except ValueError:
        print("Incorrect formatting in file, please try again")
        sys.exit()
    except IOError:
        print("File does not exist")
        sys.exit()




def print_help():
    print("SIMULATOR COMMANDS:\n")

    print("(no arguments given)")
    print("generates random sequence of length 996 with equal numbers of all nucleotides\n")

    print("--calc seq_file")
    print("  calculates nuc comp and length from seq_file and outputs random sequence with same stats\n")

    print("--calc seq_file --save my_params")
    print("  does calculations as per --calc and also saves seq statistics to my_params file\n")

    print("--load my_params")
    print("  loads seq statistics from my_params and generates random sequence with same stats\n")

    print("--load my_params --save my_params_2")
    print("  loads my_params as per --load and also saves same statistics to my_params_2\n")

    print("--help (or -h)")
    print("  outputs this help message\n")

    print("*** Any other combinations, such as --calc seq_file --load my_params, will throw an error ***")

    print("To insert a motif:\nAfter the length line in the params file, include a line as follows: 'motif: ' followed by your motif\n\nYour motif should be shorter or equal to length of sequence and also not contain any possible internal stop codons. The program will terminate if there are stop codons in the motif, even if there are other sequences to simulate.\n")

    print("To get kmer frequencies:\nFollow prompts when running --calc command. Please provide a positive whole number. Simulated sequences and kmers will be outputted after prompts")


## main block
## reads and processes commands
if (sim): ##create random seq with equal nuc composition
    random_seq()
elif (calc): ## calc stats from seq on file, generate seq
    seq_file = args[2]
    calc_stats(seq_file)
    if (k != 0):
        print("kmer frequencies:")
        total = sum(subseq_count.values())
        for key, value in subseq_count.iteritems():
            print(key + ": " + str(float(value) / total))
elif (load): ## load stats from file, generate seq
    param_file = args[2]
    load_file(param_file)
elif (usage_help): ## output help message with instructions
    print_help()
else:
    print("Unknown commands given, run with '-h' or '--help' for help")
    sys.exit()
