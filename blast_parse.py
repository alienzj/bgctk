#!/usr/bin/env python3

from Bio import SearchIO
from pprint import pprint
import pandas as pd
import sys
import math


def parse_fai(fai):
    """
   see https://samtools.github.io/hts-specs/tabix.pdf
   """
    df = pd.read_csv(
        fai,
        sep="\t",
        header=None,
        names=["gene_name", "gene_len", "gene_beg", "gene_end", "total_bytes"],
    )
    ave = sum(df["gene_len"]) / len(df)
    df["ratio"] = df["gene_len"] / ave
    df["bgc"] = df.apply(lambda x: x["gene_name"].split("|")[0], axis=1)
    return df.set_index("gene_name").loc[:, ["bgc", "gene_len", "ratio"]]


def weight(evalue):
    return math.exp(-evalue)


# BLAST tab format
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# qseqid sseqid pident length evalue

BGC_GENE = parse_fai(sys.argv[1])
i = 0
blast_df = pd.DataFrame(
    columns=["read_id", "gene_name", "ident_pct", "aln_span", "evalue"]
)
epsilon_dict = {}

for qr in SearchIO.parse(sys.stdin, format="blast-tab"):
    # print("Search %s has %i hits" % (qr.id, len(qr)))
    i += 1
    epsilon_dict[qr.id] = 1e-20
    evalue_weight_sum = 0
    for hit in qr.hits:
        for hsp in hit.hsps:
            # print("%s %s %s %s %s %s" % (qr.id, hit.id, hsp.ident_pct, hsp.evalue, hsp.aln_span, cur))
            evalue_weight = weight(hsp.evalue)
            evalue_weight_sum += evalue_weight
            blast_df = blast_df.append(
                {
                    "read_id": qr.id,
                    "gene_name": hit.id,
                    "ident_pct": hsp.ident_pct,
                    "aln_span": hsp.aln_span,
                    "evalue": hsp.evalue,
                    "evalue_weight": evalue_weight,
                },
                ignore_index=True,
            )
    epsilon_dict[qr.id] = max(epsilon_dict[qr.id], evalue_weight_sum)
    if i > 4:
        break

# blast_df = pd.merge(blast_df, BGC_GENE.reset_index(), on="gene_name")
# pprint(blast_df)
# print("\n")
# pprint(epsilon_dict)
# print("\n")
# sum_df = blast_df.groupby("read_id")["evalue_weight"].agg(max=max)
# pprint(sum_df)
# print(sum_df.loc["CL100103977L2C001R001_2429/2", "max"])

blast_df["gene_abun"] = blast_df.apply(
    lambda x: x["evalue_weight"]
    / BGC_GENE.loc[x["gene_name"], "ratio"]
    / epsilon_dict[x["read_id"]],
    axis=1
)
pprint(blast_df)
