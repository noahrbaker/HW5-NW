# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    species_files = {
        "HS": "./data/Homo_sapiens_BRD2.fa",
        "GG": "./data/Gallus_gallus_BRD2.fa",
        "MM": "./data/Mus_musculus_BRD2.fa",
        "BR": "./data/Balaeniceps_rex_BRD2.fa",
        "TT": "./data/tursiops_truncatus_BRD2.fa"
    }

    sequences = {species: read_fasta(file)[0] for species, file in species_files.items()}

    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    scores = {}
    for species, seq in sequences.items():
        if species != "HS":
            scores[species] = NW.align(sequences["HS"], seq)[0]

    # Sort the scores dictionary by value in descending order
    sorted_scores = dict(sorted(scores.items(), key=lambda x: x[1], reverse=True))

    print("Alignment scores of species to human BRD2:")
    for species, score in sorted_scores.items():
        print(f"{species}: {score}")

if __name__ == "__main__":
    main()