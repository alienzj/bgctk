# /usr/bin/env python3

import argparse
import concurrent.futures
import math
import os
import sys
import time
from pprint import pprint

import pandas as pd
from Bio import SearchIO
from scipy import stats


def count_no_zero(x):
    return x[x > 0.0].count()


def abun(x):
    s = x[x > 0.0]
    if not s.empty:
        return stats.hmean(s)
    else:
        return 0.0


def abun_global(x):
    if x["count"] / 2 <= x["count_no_zero"]:
        return x["abun"]
    else:
        return 0.0


def abun_strict(x):
    if BGC_COUNT.loc[x["bgc"], "count"] / 2 <= x["count_no_zero"]:
        return x["abun"]
    else:
        return 0.0


def abun_strict_v2(x):
    # if BGC_COUNT.loc[x["bgc"], "count"] / 2 <= BGC_COUNT2[x["bgc"]]["count"]:
    if BGC_COUNT.loc[x["bgc"], "count"] / 2 <= x["count_no_zero_v2"]:
        return x["abun"]
    else:
        return 0.0


def parse(profile_file):
    sample_id = os.path.basename(profile_file).split(".")[0]
    try:
        if os.path.exists(profile_file):
            abun_df = pd.read_csv(profile_file, sep="\t")
            if abun_df.empty:
                print("%s is empty" % profile_file)
                return None, None, None, None, None
        else:
            print("%s is not exists" % profile_file)
            return None, None, None, None, None
    except pd.io.common.EmptyDataError:
        print("%s is empty" % profile_file)
        return None, None, None, None, None

    count_df1 = (
        abun_df.loc[:, ["bgc", "count_no_zero"]]
        .rename(columns={"count_no_zero": sample_id})
        .set_index("bgc")
    )
    count_df2 = (
        abun_df.loc[:, ["bgc", "count_no_zero_v2"]]
        .rename(columns={"count_no_zero_v2": sample_id})
        .set_index("bgc")
    )
    abun_df1 = (
        abun_df.loc[:, ["bgc", "abun_norm"]]
        .rename(columns={"abun_norm": sample_id})
        .set_index("bgc")
    )
    abun_df2 = (
        abun_df.loc[:, ["bgc", "abun_norm_strict"]]
        .rename(columns={"abun_norm_strict": sample_id})
        .set_index("bgc")
    )
    abun_df3 = (
        abun_df.loc[:, ["bgc", "abun_norm_strict_v2"]]
        .rename(columns={"abun_norm_strict_v2": sample_id})
        .set_index("bgc")
    )
    return count_df1, count_df2, abun_df1, abun_df2, abun_df3


def parse_fai(fai):
    """
   see https://samtools.github.io/hts-specs/tabix.pdf
   """
    df = pd.read_csv(
        fai,
        sep="\t",
        header=None,
        names=["gene_name", "gene_len", "gene_beg", "gene_end", "total_bytes"],
        encoding="utf8",
    )
    ave = sum(df["gene_len"]) / len(df)
    df["ratio"] = df["gene_len"] / ave
    df["bgc"] = df.apply(lambda x: x["gene_name"].split("|")[0], axis=1)
    return (
        df.set_index("gene_name").loc[:, ["bgc", "gene_len", "ratio"]],
        df.groupby("bgc")["gene_name"].agg(count="count"),
    )


def weight(evalue):
    return math.exp(-evalue)


def profiler(args):
    start_time = time.time()
    print("parse begin")

    global BGC_COUNT
    BGC_GENE, BGC_COUNT = parse_fai(args.bgc)

    BGC_COUNT2 = {}

    blast_dict = {
        "read_id": [],
        "gene_name": [],
        "ident_pct": [],
        "aln_span": [],
        "evalue": [],
        "evalue_weight": [],
    }
    epsilon_dict = {}
    # i = 0

    for qr in SearchIO.parse(args.files, format="blast-tab"):
        # print("Search %s has %i hits" % (qr.id, len(qr)))
        # i += 1
        epsilon_dict[qr.id] = 1e-20
        evalue_weight_sum = 0
        for hit in qr.hits:
            for hsp in hit.hsps:
                if (hsp.evalue <= 1e-5) and (hsp.bitscore >= 10):
                    evalue_weight = weight(hsp.evalue)
                    evalue_weight_sum += evalue_weight

                    # Pandas/numpy structures are fundamenally not suited for efficiently growing
                    blast_dict["read_id"].append(qr.id)
                    blast_dict["gene_name"].append(hit.id)
                    blast_dict["ident_pct"].append(hsp.ident_pct)
                    blast_dict["aln_span"].append(hsp.aln_span)
                    blast_dict["evalue"].append(hsp.evalue)
                    blast_dict["evalue_weight"].append(evalue_weight)

        epsilon_dict[qr.id] = max(epsilon_dict[qr.id], evalue_weight_sum)
        # if i > 10000:
        #     break

    blast_df = pd.DataFrame(blast_dict)
    blast_df["gene_abun"] = blast_df.apply(
        lambda x: x["evalue_weight"]
        / BGC_GENE.loc[x["gene_name"], "ratio"]
        / epsilon_dict[x["read_id"]],
        axis=1,
    )
    bgc_df = pd.merge(blast_df, BGC_GENE, on="gene_name")

    for i in range(len(blast_df)):
        gene_abun = blast_df.at[i, "gene_abun"]
        if gene_abun > 0.0:
            gene_name = blast_df.at[i, "gene_name"]
            bgc_id = gene_name.split("|")[0]
            if not bgc_id in BGC_COUNT2:
                BGC_COUNT2[bgc_id] = {}
                BGC_COUNT2[bgc_id][gene_name] = 1
                BGC_COUNT2[bgc_id]["count"] = 1
            else:
                if not gene_name in BGC_COUNT2[bgc_id]:
                    BGC_COUNT2[bgc_id][gene_name] = 1
                    BGC_COUNT2[bgc_id]["count"] += 1
                else:
                    BGC_COUNT2[bgc_id][gene_name] += 1

    abun_df = (
        bgc_df.groupby("bgc")["gene_abun"]
        .agg(count="count", count_no_zero=count_no_zero, abun=abun)
        .reset_index()
    )

    abun_df["count_no_zero_v2"] = abun_df.apply(
        lambda x: BGC_COUNT2[x["bgc"]]["count"], axis=1
    )

    abun_df["abun_strict"] = abun_df.apply(lambda x: abun_strict(x), axis=1)
    abun_df["abun_strict_v2"] = abun_df.apply(lambda x: abun_strict_v2(x), axis=1)

    if sum(abun_df["abun"]) != 0.0:
        abun_df["abun_norm"] = abun_df["abun"] / sum(abun_df["abun"])
    else:
        abun_df["abun_norm"] = 0.0
    if sum(abun_df["abun_strict"]) != 0.0:
        abun_df["abun_norm_strict"] = abun_df["abun_strict"] / sum(
            abun_df["abun_strict"]
        )
    else:
        abun_df["abun_norm_strict"] = 0.0
    if sum(abun_df["abun_strict_v2"]) != 0.0:
        abun_df["abun_norm_strict_v2"] = abun_df["abun_strict_v2"] / sum(
            abun_df["abun_strict_v2"]
        )
    else:
        abun_df["abun_norm_strict_v2"] = 0.0

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    abun_df.to_csv(
        os.path.join(args.outdir, args.sample_id + ".profile"), sep="\t", index=False
    )
    print("parse end")
    print("%f seconds" % (time.time() - start_time))


def merger(args):
    count1_list = []
    count2_list = []
    abun1_list = []
    abun2_list = []
    abun3_list = []

    profile_list = [line.strip() for line in open(args.files, "r").readlines()]

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        for count1_df, count2_df, abun1_df, abun2_df, abun3_df in executor.map(
            parse, profile_list
        ):
            if (
                (count1_df is not None)
                and (count2_df is not None)
                and (abun1_df is not None)
                and (abun2_df is not None)
                and (abun3_df is not None)
            ):
                count1_list.append(count1_df)
                count2_list.append(count2_df)
                abun1_list.append(abun1_df)
                abun2_list.append(abun2_df)
                abun3_list.append(abun3_df)

    count1_df_ = pd.concat(count1_list, axis=1, sort=True).fillna(0).reset_index()
    count2_df_ = pd.concat(count2_list, axis=1, sort=True).fillna(0).reset_index()

    abun1_df_ = pd.concat(abun1_list, axis=1, sort=True).fillna(0).reset_index()
    abun2_df_ = pd.concat(abun2_list, axis=1, sort=True).fillna(0).reset_index()
    abun3_df_ = pd.concat(abun3_list, axis=1, sort=True).fillna(0).reset_index()

    count1_df_.to_csv(os.path.join(args.outdir, "count.profile"), sep="\t", index=False)
    count2_df_.to_csv(
        os.path.join(args.outdir, "count_v2.profile"), sep="\t", index=False
    )

    abun1_df_.to_csv(
        os.path.join(args.outdir, "abundance.profile"), sep="\t", index=False
    )
    abun2_df_.to_csv(
        os.path.join(args.outdir, "abundance_strict.profile"), sep="\t", index=False
    )
    abun3_df_.to_csv(
        os.path.join(args.outdir, "abundance_strict_v2.profile"), sep="\t", index=False
    )


def main():
    parser = argparse.ArgumentParser(
        prog="bgctk",
        usage="bgctk [subcommand] [options]",
        description="bgctk, a toolkit to construct BGC profile",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print version and exit",
    )

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "--files",
        default=sys.stdin,
        help="files to read, blast-out or profile-out, if empty, stdin is used",
    )

    subparsers = parser.add_subparsers(title="available subcommands", metavar="")

    parser_profiler = subparsers.add_parser(
        "profiler",
        parents=[parent_parser],
        prog="bgctk profiler",
        description="blastout praser",
        help="parse blastout to compute BGC profile",
    )
    parser_profiler.add_argument(
        "--bgc", required=True, type=str, help="bgc gene metadata"
    )
    parser_profiler.add_argument(
        "--sample_id", required=True, type=str, help="sample id"
    )
    parser_profiler.add_argument(
        "--outdir", required=True, type=str, help="bgc profile output dir"
    )

    parser_merger = subparsers.add_parser(
        "merger",
        parents=[parent_parser],
        prog="bgctk merger",
        description="bgc profile merger",
        help="merge bgc profile",
    )
    parser_merger.add_argument(
        "--outdir", required=True, type=str, help="merged bgc profile output dir"
    )
    parser_merger.add_argument(
        "--threads", type=int, default=8, help="threads, default: 8"
    )

    parser_profiler._optionals.title = "arguments"
    parser_profiler.set_defaults(func=profiler)
    parser_merger._optionals.title = "arguments"
    parser_merger.set_defaults(func=merger)

    args = parser.parse_args()
    try:
        if args.version:
            print("bgctk version 0.1.0")
            sys.exit(0)
        args.func(args)
    except AttributeError as e:
        print(e)
        parser.print_help()


if __name__ == "__main__":
    main()
