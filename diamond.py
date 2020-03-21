#!/usr/bin/env python3

import os
import sys
import pandas as pd
from metapi import sample


MIBIG_DB = "/ldfssz1/ST_META/share/User/zhujie/database/mibig/mibig_prot_seqs_2.0.dmnd"
MIBIG_FA = "/ldfssz1/ST_META/share/User/zhujie/database/mibig/mibig_prot_seqs_2.0.fasta"
BGCTK = (
    "/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/igo/assay/function/bgc_profile/bgctk.py"
)
BGC_FAI = (
    "/ldfssz1/ST_META/share/User/zhujie/database/mibig/mibig_prot_seqs_2.0.fasta.fai"
)
VFDB = "/ldfssz1/ST_META/share/User/zhujie/database/VFDB/VFDB_setB_pro.fas.gz.dmnd"
VFDB_FAI = "/ldfssz1/ST_META/share/User/zhujie/database/VFDB/VFDB_setB_pro.fas.fai"
BACTIBASE_DB = "/ldfssz1/ST_META/share/User/zhujie/database/bactibase/bactibase.fa.dmnd"
BACTIBASE_DB_FAI = (
    "/ldfssz1/ST_META/share/User/zhujie/database/bactibase/bactibase.fa.fai"
)
PROFILING_TK = "/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/igo/assay/function/bacteriocins_vfs/profiling.py"


if len(sys.argv) != 3:
    print(
        "python %s [samples.tsv] [jgi|blast|bgctk|profiling_bactibase|profiling_vfdb]"
        % os.path.abspath(__file__)
    )
    sys.exit(1)
elif (
    (sys.argv[2] != "jgi")
    and (sys.argv[2] != "blast")
    and (sys.argv[2] != "bgctk")
    and (sys.argv[2] != "profiling_bactibase")
    and (sys.argv[2] != "profiling_vfdb")
):
    print(
        "python %s [samples.tsv] [jgi|blast|bgctk|profiling_bactibase|profiling_vfdb]"
        % os.path.abspath(__file__)
    )
    sys.exit(1)

SAMPLES = sample.parse_samples(sys.argv[1], "fastq", True, False)

blastx_tab = "blastx_tab"
blastx_logs = "blastx_logs"

blastx_bam = "blastx_bam"
jgi_out = "jgi_out"

bgctk_out = "bgctk_out"
bgctk_logs = "bgctk_logs"

profiling_out = "profiling_out"
profiling_logs = "profiling_logs"

if sys.argv[2] == "jgi":
    os.makedirs(blastx_bam, exist_ok=True)
    os.makedirs(jgi_out, exist_ok=True)
    os.makedirs(blastx_logs, exist_ok=True)

if sys.argv[2] == "blast":
    os.makedirs(blastx_tab, exist_ok=True)
    os.makedirs(blastx_logs, exist_ok=True)

if sys.argv[2] == "bgctk":
    os.makedirs(bgctk_out, exist_ok=True)
    os.makedirs(bgctk_logs, exist_ok=True)

if (sys.argv[2] == "profiling_bactibase") or (sys.argv[2] == "profiling_vfdb"):
    os.makedirs(profiling_out, exist_ok=True)
    os.makedirs(profiling_logs, exist_ok=True)


for id in SAMPLES.index.unique():
    fq1 = SAMPLES.loc[id, "fq1"]
    fq2 = SAMPLES.loc[id, "fq2"]

    if sys.argv[2] == "jgi":
        if not pd.isna(fq2):
            cmd_1 = (
                "seqtk mergepe %s %s | diamond blastx -d %s -o %s/%s.sam -f 101 -p 4 -e 1e-05 > %s/%s.diamond.log"
                % (fq1, fq2, MIBIG_DB, blastx_bam, id, blastx_logs, id)
            )
        else:
            cmd_1 = (
                "diamond blastx -q %s -d %s -o %s/%s.sam -f 101 -p 4 -e 1e-05 > %s/%s.diamond.log"
                % (fq1, MIBIG_DB, blastx_bam, id, blastx_logs, id)
            )
        print(cmd_1)

        cmd_2 = (
            "samtools view -bhT %s %s/%s.sam | samtools sort -@4 -T %s/%s.temp -O BAM - | jgi_summarize_bam_contig_depths --outputDepth %s/%s.depth -"
            % (MIBIG_FA, blastx_bam, id, blastx_bam, id, jgi_out, id)
        )
        print(cmd_2)

        cmd_3 = "gzip -f %s/%s.depth" % (jgi_out, id)
        cmd_4 = "rm -rf %s/%s.sam %s/%s.temp.*.bam" % (blastx_bam, id, blastx_bam, id)
        print(cmd_3)
        print(cmd_4)

    elif sys.argv[2] == "blast":
        if not pd.isna(fq2):
            cmd_1 = (
                "seqtk mergepe %s %s | diamond blastx -d %s -o %s/%s.blast.m6 -f 6 -p 4 -e 1e-05 > %s/%s.diamond.log"
                % (fq1, fq2, MIBIG_DB, blastx_tab, id, blastx_logs, id)
            )
        else:
            cmd_1 = (
                "diamond blastx -q %s -d %s -o %s/%s.blast.m6 -f 6 -p 4 -e 1e-05 > %s/%s.diamond.log"
                % (fq1, MIBIG_DB, blastx_tab, id, blastx_logs, id)
            )
        print(cmd_1)

    elif sys.argv[2] == "bgctk":
        if not pd.isna(fq2):
            cmd_1 = (
                "seqtk mergepe %s %s | diamond blastx -d %s -f 6 -p 4 -e 1e-05 | python %s profiler --bgc %s --sample_id %s --outdir %s > %s/%s.bgctk.log"
                % (fq1, fq2, MIBIG_DB, BGCTK, BGC_FAI, id, bgctk_out, bgctk_logs, id)
            )
        else:
            cmd_1 = (
                "diamond blastx -q %s -d %s -f 6 -p 4 -e 1e-05 | python %s profiler --bgc %s --sample_id %s --outdir %s > %s/%s.bgctk.log"
                % (fq1, MIBIG_DB, BGCTK, BGC_FAI, id, bgctk_out, bgctk_logs, id)
            )
        print(cmd_1)

    elif sys.argv[2] == "profiling_bactibase":
        if not pd.isna(fq2):
            cmd_1 = (
                "seqtk mergepe %s %s | diamond blastx -d %s -f 6 -p 4 -e 1e-05 | python %s profiler --metadata %s --sample_id %s --outdir %s > %s/%s.profiling.log"
                % (
                    fq1,
                    fq2,
                    BACTIBASE_DB,
                    PROFILING_TK,
                    BACTIBASE_DB_FAI,
                    id,
                    profiling_out,
                    profiling_logs,
                    id,
                )
            )
        else:
            cmd_1 = (
                "diamond blastx -q %s -d %s -f 6 -p 4 -e 1e-05 | python %s profiler --metadata %s --sample_id %s --outdir %s > %s/%s.profiling.log"
                % (
                    fq1,
                    BACTIBASE_DB,
                    PROFILING_TK,
                    BACTIBASE_DB_FAI,
                    id,
                    profiling_out,
                    profiling_logs,
                    id,
                )
            )
        print(cmd_1)

    elif sys.argv[2] == "profiling_vfdb":
        if not pd.isna(fq2):
            cmd_1 = (
                "seqtk mergepe %s %s | diamond blastx -d %s -f 6 -p 4 -e 1e-05 | python %s profiler --metadata %s --sample_id %s --outdir %s > %s/%s.profiling.log"
                % (
                    fq1,
                    fq2,
                    VFDB,
                    PROFILING_TK,
                    VFDB_FAI,
                    id,
                    profiling_out,
                    profiling_logs,
                    id,
                )
            )
        else:
            cmd_1 = (
                "diamond blastx -q %s -d %s -f 6 -p 4 -e 1e-05 | python %s profiler --metadata %s --sample_id %s --outdir %s > %s/%s.profiling.log"
                % (
                    fq1,
                    VFDB,
                    PROFILING_TK,
                    VFDB_FAI,
                    id,
                    profiling_out,
                    profiling_logs,
                    id,
                )
            )
        print(cmd_1)

    else:
        sys.exit(1)
