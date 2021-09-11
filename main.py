

def input_and_replace(inputfile):
    file: str = inputfile
    f = open(file, "r")
    seq = f.read()

    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    seq = seq.replace("U", "T")
    #print(seq)
    return seq

def triplet_to_aminoacid(seq):
    aa_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    protein = []
    for i in range(len(seq)):
        codon = (seq[i+1 : i+4])
        if len(codon) != 3:
          break
        #print(codon)
        protein.append(aa_table[codon])
        i += 3

    #print(protein)
    return protein

def out_file(protein, filename):
    a_list = protein
    with open(filename, "w") as f: 
      for item in a_list:
        f.write("%s" % item)

def index_gene(gene):
    max = len(gene) + 1

    temp = []
    temp = [0 for i in range(len(gene))]
    j = 0
    i = 1
    while i < len(gene):
        if gene[i] == gene[j]:
            temp[i] = j + 1
            j += 1
            i += 1
        elif i != j and j > 0:
            j = temp[j - 1]
        elif i == j and j == 0:
            temp[i] = j + 1
            i += 1
        elif i != j and j == 0:
            temp[i] = 0
            i += 1
    return temp

def gene_search(seq, gene, temp):
    i = 0
    j = 0
    while j < len(seq):
        if i == len(gene):
            break
        j += 1
        if gene[i] == seq[j - 1]:
            i += 1
        elif gene[i] != seq[j - 1]:
            i = temp[i - 1]

    if i == len(gene):
        position = j - len(gene)
        print("Gene found at position: " + str(position + 1))
        return


def main():
    inputfile = "DNA_seq_original.txt"
    seq = input_and_replace(inputfile)
    out_file(seq, "DNA processed seq.txt")
    x = "-1"
    while x != "0":
        x = input("For Amino acid conversion enter 1,\n"
                  "For Gene searching enter 2,\n"
                  "To end enter 0: ")
        if x == "1":
            protein = triplet_to_aminoacid(seq)
            outfile = input("Please enter the out file name: ")
            out_file(protein, "" + outfile + ".txt")
            print("Amino Acid conversion complete.")

        if x == "2":
            gene = input("Enter the gene sequence to search for: ")
            temp = index_gene(gene)
            gene_search(seq, gene, temp)

        if int(x) < 0 or int(x) > 2:
            print("Value entered must be 1, 2 or 0 for exit. Please try again.")






if __name__ == '__main__':
    main()

