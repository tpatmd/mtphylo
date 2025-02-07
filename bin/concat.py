#!/usr/bin/env python3

import sys

if sys.argv[2] == "13":
    genes = [
        "nad2",
        "cox1",
        "cox2",
        "atp8",
        "atp6",
        "cox3",
        "nad3",
        "nad5",
        "nad4",
        "nad4l",
        "nad6",
        "cytb",
        "nad1",
    ]
elif sys.argv[2] == "12":
    genes = [
        "nad2",
        "cox1",
        "cox2",
        "atp6",
        "cox3",
        "nad3",
        "nad5",
        "nad4",
        "nad4l",
        "nad6",
        "cytb",
        "nad1",
    ]

sequences = {}
length = {}
tmp = 0

if sys.argev[4] == "nuc":
    with open(sys.argv[1], "r") as file:
        for line in file:
            if line.strip():
                if line.startswith(">"):
                    species = "_".join(line.strip().split("_")[:3])[1:]
                    gene = line.strip().split("_")[3]
                    sequences[species] = sequences.get(species, {})
                else:
                    sequences[species][gene] = (
                        sequences[species].get(gene, "")
                        + "".join(line.strip().split(" ")).upper()
                    )
elif sys.argv[4] == "aa":
    with open(sys.argv[1], "r") as file:
        for line in file:
            if line.strip():
                if line.startswith(">"):
                    species = line.strip().split(" ")[0][1:]
                    gene = line.strip().split(" ")[1]
                    sequences[species] = sequences.get(species, {})
                else:
                    sequences[species][gene] = (
                        sequences[species].get(gene, "")
                        + "".join(line.strip().split(" ")).upper()
                    )

with open(f"concatenated_dataset.fasta", "w") as file:
    for species, sequence in sequences.items():
        file.write(f">{species}\n")
        for gene in genes:
            file.write(f"{sequences[species][gene]}")
            if length.get(gene, 0) == 0:
                length[gene] = len(sequences[species][gene])
        file.write(f"\n")

if sys.argv[3] == "cds":
    with open(f"concatenated_dataset.partition", "w") as file:
        start = 1
        file.write(f"#nexus\n")
        file.write(f"begin sets;\n")
        for gene in genes:
            file.write(f"    charset {gene} = {start}-{start + length[gene] - 1};\n")
            start += length[gene]
        file.write(f"end;")
elif sys.argv[3] == "codon":
    with open(f"concatenated_dataset.partition", "w") as file:
        start = 1
        file.write(f"#nexus\n")
        file.write(f"begin sets;\n")
        for gene in genes:
            file.write(
                f"    charset {gene}-pt1 = {start}-{start + length[gene] - 1}\\3;\n"
            )
            file.write(
                f"    charset {gene}-pt2 = {start + 1}-{start + length[gene] - 1}\\3;\n"
            )
            file.write(
                f"    charset {gene}-pt3 = {start + 2}-{start + length[gene] - 1}\\3;\n"
            )
            start += length[gene]
        file.write(f"end;")
