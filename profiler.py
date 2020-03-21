#/usr/bin/env python3

import pandas as pd
from scipy import stats
import os
import sys
import concurrent.futures


def parse(depth_file):
    sample_id = os.path.basename(depth_file).split(".")[0]
    try:
        if os.path.exists(depth_file):
            df = pd.read_csv(depth_file, sep='\t')
            if df.empty:
                print("%s is empty" % depth_file)
                return None, None
        else:
            print("%s is not exists" % depth_file)
            return None, None
    except pd.io.common.EmptyDataError:
        print("%s is empty" % depth_file)
        return None, None

    df["BGC_id"] = df.apply(lambda x: x["contigName"].split("|")[0], axis=1)
    df["totalAvgDepth_norm"] = df["totalAvgDepth"] / df["contigLen"]

    def count_no_zero(x):
        return x[x > 0.0].count()

    def abun(x):
         s = x[x > 0.0]
         if not s.empty:
             return stats.hmean(s)
         else:
             return 0.0

    def abun2(x):
        if x["count"] / 2 <= x["count_no_zero"]:
            return x["abun"]
        else:
            return 0.0

    result = df.groupby(["BGC_id"])["totalAvgDepth_norm"].agg(['count', count_no_zero, abun])
    result["abun2"] = result.apply(lambda x: abun2(x), axis=1)
    result["abun_norm"] = result["abun"] / sum(result["abun"])
    result["abun2_norm"] = result["abun2"] / sum(result["abun2"])
    result = result.reset_index()

    return result.loc[:, ["BGC_id", "count_no_zero"]].rename(columns={"count_no_zero": sample_id}).set_index("BGC_id"), \
           result.loc[:, ["BGC_id", "abun_norm"]].rename(columns={"abun_norm": sample_id}).set_index("BGC_id"), \
           result.loc[:, ["BGC_id", "abun2_norm"]].rename(columns={"abun2_norm": sample_id}).set_index("BGC_id")


def merge(depth_files, workers):
    count_list = []
    abun_list = []
    abun2_list = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for count_df, abun_df, abun2_df in executor.map(parse, depth_files):
            if (count_df is not None) and (abun_df is not None) and (abun2_df is not None):
                count_list.append(count_df)
                abun_list.append(abun_df)
                abun2_list.append(abun2_df)

    count_df_ = pd.concat(count_list, axis=1).reset_index()
    abun_df_ = pd.concat(abun_list, axis=1).reset_index()
    abun2_df_ = pd.concat(abun2_list, axis=1).reset_index()

    return count_df_, abun_df_, abun2_df_


def main():
    depth_files = [x.strip() for x in open(sys.argv[1], 'r').readlines()]
    count_df, abun_df, abun2_df = merge(depth_files, int(sys.argv[2]))

    count_df.to_csv("count.profile", sep='\t', index=False)
    abun_df.to_csv("abundance.profile", sep='\t', index=False)
    abun2_df.to_csv("abundance_strict.profile", sep='\t', index=False)


if __name__ == '__main__':
    main()
