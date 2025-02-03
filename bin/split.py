#!/usr/bin/env python3

import sys

cds = {}
genes = [
    "atp6",
    "atp8",
    "cox1",
    "cox2",
    "cox3",
    "cytb",
    "nad1",
    "nad2",
    "nad3",
    "nad4",
    "nad4l",
    "nad5",
    "nad6",
]

with open(sys.argv[1], "r") as file:
    for line in file:
        if line.strip():
            if line.startswith(">"):
                species = "_".join(line.strip().split(" ")[:3])[1:]
                gene = line.strip().split(" ")[3]
                cds[species] = cds.get(species, {})
            else:
                cds[species][gene] = cds[species].get(gene, "") + line.strip()

for gene in genes:
    with open(f"{gene}.fasta", "w") as file:
        for species, seq in cds.items():
            if cds[species].get(gene):
                file.write(f">{species} {gene}\n{cds[species][gene]}\n")

with open("split-report.tsv", "w") as file:
    file.write("accession no.\tspecies\t" + "\t".join(genes) + "\n")
    for species, seq in cds.items():
        file.write(f"{species.split("_")[0]}\t{" ".join(species.split("_")[1:3])}")
        for gene in genes:
            if cds[species].get(gene):
                file.write(f"\t{len(cds[species][gene])}")
            else:
                file.write(f"\t0")
        file.write("\n")

print(len(cds.items()))