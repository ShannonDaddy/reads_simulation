#!/usr/bin/python
import argparse
import textwrap
import pandas as pd
import gzip
import os
import sys
import random
import re

# initialize
amplicon_list = ""
variant_table = ""
art_path = "/home/lcbio/Softwares/art_bin_MountRainier/art_illumina"
art_opts = "-amp --paired -qL 15 -ss MSv3 -l 151 -m 250 -s 50 -f 20"
molecule_num = 6000
# molecule_pcr_fold = 20
read1_out = "sim_read1.fastq.gz"
read2_out = "sim_read2.fastq.gz"

SEQ1_BEYOND = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCAAGCTTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA"
SEQ2_BEYOND = "GATCGTCGGACTGTAGAAATCTGAACGTGTACATCAGCGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAA"

def parse_arguments():
    global amplicon_list
    global variant_table
    global art_path
    global art_opts
    global molecule_num
    # global molecule_pcr_fold
    global read1_out
    global read2_out

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='get simulated reads',
                                     epilog=textwrap.dedent('''
                                         version:    1.0
                                         contact:    gqian@lc-bio.com
                                         date:       2021/02/16'''))

    parser.add_argument("--amplicon_list", metavar="amplicon_list", help="the amplicon list", required=True)
    parser.add_argument("--variant_table", metavar="var_table", help="the variant table")
    parser.add_argument("--art_path", metavar="ART_PATH", help="path to the art simulator, default: "
                                                               "/home/lcbio/Softwares/art_bin_MountRainier/art_illumina")
    parser.add_argument("--art_opts", metavar="ART_OPTS",
                        help="the art options, default: '-amp --paired -qL 15 -ss MSv3 -l 151 -m 250 -s 50 -f 20'")
    parser.add_argument("--molecule_num", metavar="molecule_num", help="raw molecule number, default: 6000")
    # parser.add_argument("--molecule_pcr_fold", metavar="molecule_pcr_fold", help="PCR fold coverage of one
    # molecule, default: 20")
    parser.add_argument("--read1_out", metavar="read1_out", help="read1 fastq file, default: ./sim_read1.fastq.gz")
    parser.add_argument("--read2_out", metavar="read2_out", help="read2 fastq file, default: ./sim_read2.fastq.gz")

    args = parser.parse_args()
    argsDict = args.__dict__

    amplicon_list = argsDict['amplicon_list']

    if argsDict['variant_table']:
        variant_table = argsDict['variant_table']

    if argsDict['art_path']:
        art_path = argsDict['art_path']

    if argsDict['art_opts']:
        art_opts = argsDict['art_opts']

    if argsDict['molecule_num']:
        molecule_num = int(argsDict['molecule_num'])

    # if argsDict['molecule_pcr_fold']:
    #     molecule_pcr_fold = argsDict['molecule_pcr_fold']

    if argsDict['read1_out']:
        read1_out = argsDict['read1_out']

    if argsDict['read2_out']:
        read2_out = argsDict['read2_out']


def fill_rand_base_in_mtag(mtag):

    nuc = ['A', 'T', 'C', 'G']

    bases = []
    for i in mtag:
        if i == 'N':
            bases.append(nuc[random.randint(0,3)])
        else:
            bases.append(i)

    mtag_filled = ''.join(bases)

    return mtag_filled


def fill_rand_base_in_mtag_pair(mtag1, mtag2, mtag_recorded):

    mtag1_filled = fill_rand_base_in_mtag(mtag1)
    mtag2_filled = fill_rand_base_in_mtag(mtag2)
    mtag_pair = mtag1_filled + mtag2_filled

    while mtag_recorded.has_key(mtag_pair):
        mtag1_filled = fill_rand_base_in_mtag(mtag1)
        mtag2_filled = fill_rand_base_in_mtag(mtag2)
        mtag_pair = mtag1_filled + mtag2_filled

        if not mtag_recorded.has_key(mtag_pair):
            mtag_recorded[mtag_pair] = 1
            break

    if not mtag_recorded.has_key(mtag_pair):
            mtag_recorded[mtag_pair] = 1

    return mtag1_filled, mtag2_filled


def get_var_seq(prbSeq, start_pos, ref, alt):

    bases = list(prbSeq)
    start = start_pos - 1
    end = start + len(ref)
    bases[start:end] = alt
    var_seq = ''.join(bases)

    return var_seq


def get_var_prbSeq(chr, prbStart, prbEnd, prbStrand, prbSeq, var):
    var_prbSeq = ''
    var_chr = var["chr"]
    var_pos = var["pos"]
    var_ref = var["ref"]
    var_alt = var["alt"]
    var_id = "%s|%d|%s|%s" % (var_chr, var_pos, var_ref, var_alt)
    var_ref_len = len(var_ref)

    if (chr == var_chr) and (prbStart <= var_pos) and ((var_pos + var_ref_len - 1) <= prbEnd):
        if prbStrand == "+":
            var_prbSeq = get_var_seq(prbSeq, var_pos-prbStart+1, var_ref, var_alt)
        else:
            prbSeq_rc = reverse_comp(prbSeq)
            var_prbSeq = reverse_comp(get_var_seq(prbSeq_rc, var_pos-prbStart+1, var_ref, var_alt))
    else:
        sys.stderr.write(var_id + " not covered by the specified amplicon, please check you variant table settings!\n")
        sys.exit(-1)

    return var_prbSeq


def reverse_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    return "".join(complement.get(base, base) for base in reversed(seq))


def get_amp_seq(df_amp, ampIndex, molecule_num, mtag_recorded, out_prefix, read_len, var=None):
    amp_seq_file1 = out_prefix + "1.fasta"
    amp_seq_file2 = out_prefix + "2.fasta"
    fout1 = open_file(amp_seq_file1, 'w')
    fout2 = open_file(amp_seq_file2, 'w')

    ampID = df_amp.loc[ampIndex, "ampID"]
    chr = df_amp.loc[ampIndex, "chr"]
    prbStart = df_amp.loc[ampIndex, "prbStart"]
    prbEnd = df_amp.loc[ampIndex, "prbEnd"]
    prbStrand = df_amp.loc[ampIndex, "prbStrand"]
    prbSeq = df_amp.loc[ampIndex, "prbSeq"]

    mtag1 = ''
    mtag2 = ''
    has_mtag = True
    if 'seqMTag1' in df_amp.columns:
        mtag1 = df_amp.loc[ampIndex, "seqMTag1"]
    if 'seqMTag2' in df_amp.columns:
        mtag2 = df_amp.loc[ampIndex, "seqMTag2"]
    if (mtag1 == '') or (mtag2 == ''):
        has_mtag = False

    if var is not None:
        prbSeq = get_var_prbSeq(chr, prbStart, prbEnd, prbStrand, prbSeq, var)

    prbSeq_rc = reverse_comp(prbSeq)
    SEQ2_BEYOND_RC = reverse_comp(SEQ2_BEYOND)
    for i in range(0, molecule_num):
        if has_mtag:
            mtag1_filled, mtag2_filled = fill_rand_base_in_mtag_pair(mtag1, mtag2, mtag_recorded)
            mtag1_filled_rc = reverse_comp(mtag1_filled)
            name = '_'.join([ampID, "MTagR1", "MTagR2", mtag2_filled, mtag1_filled])
            # seq1 = SEQ1_BEYOND_RC + mtag1_filled + prbSeq + mtag2_filled_rc
            seq1 = mtag2_filled + prbSeq_rc + mtag1_filled_rc + SEQ1_BEYOND
            # seq2 = mtag1_filled + prbSeq + mtag2_filled_rc + SEQ2_BEYOND
            seq2 = SEQ2_BEYOND_RC + mtag2_filled + prbSeq_rc + mtag1_filled_rc
        else:
            name = '_'.join([ampID, str(i+1)])
            seq1 = prbSeq_rc + SEQ1_BEYOND
            seq2 = SEQ2_BEYOND_RC + prbSeq_rc
        seq1 += 'A' * (read_len - len(seq1))
        seq2 = 'T' * (read_len - len(seq2)) + seq2
        fout1.write(">" + name + "\n" + seq1 + "\n")
        fout2.write(">" + name + "\n" + seq2 + "\n")

    fout1.close()
    fout2.close()

    return amp_seq_file1, amp_seq_file2


def get_sim_reads(amp_seq1, amp_seq2, out_prefix):

    out_prefix1 = out_prefix + "1."
    os.system("{art_path} -i {amp_seq} {art_opts} -o {out_prefix}".format(art_path=art_path,
                                                                          amp_seq=amp_seq1,
                                                                          art_opts=art_opts,
                                                                          out_prefix=out_prefix1))
    sim1_r1 = out_prefix1 + "1.fq"
    sim1_r2 = out_prefix1 + "2.fq"

    out_prefix2 = out_prefix + "2."
    os.system("{art_path} -i {amp_seq} {art_opts} -o {out_prefix}".format(art_path=art_path,
                                                                          amp_seq=amp_seq2,
                                                                          art_opts=art_opts,
                                                                          out_prefix=out_prefix2))
    sim2_r1 = out_prefix2 + "1.fq"
    sim2_r2 = out_prefix2 + "2.fq"

    return sim1_r1, sim2_r2


def open_file(input_file, mode):
    assert isinstance(input_file, str)
    if input_file.endswith(".gz"):
        fin = gzip.open(input_file, mode)
    else:
        fin = open(input_file, mode)
    return fin


def write_sim_reads(sim_reads, fw):
    fin = open_file(sim_reads, 'r')

    for line in fin:
        fw.write(line)

    fin.close()


def get_read_length(art_opts):
    m = re.search(r'-l\s+(\d+)', art_opts)
    if m:
        read_len = int(m.group(1))
    else:
        sys.stderr.write("[ERROR] can't parse from art options for read length, please check your art options!\n")
        sys.exit(-1)

    return read_len


def simulate_reads(df_amp, df_var):

    read_len = get_read_length(art_opts)
    ampIndex_list = sorted(list(set(df_amp.index)))
    out_dir = os.path.dirname(os.path.realpath(read1_out))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    tmp_dir = out_dir + "/sim_tmp"
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    fout1 = open_file(read1_out, 'w')
    fout2 = open_file(read2_out, 'w')
    for ind, ampIndex in enumerate(ampIndex_list):
        mtag_recorded = {}
        vars = pd.DataFrame()
        if ampIndex in df_var.index:
            vars = df_var.loc[ampIndex, :]

        var_seq_num = 0
        var_num = 1
        if isinstance(vars, pd.DataFrame):
            var_num = len(vars)

        for i in range(0, var_num):
            var = vars
            if var_num > 1:
                var = vars.iloc[i]
            seq_num = int(molecule_num * float(var["vaf"]))
            var_seq_num += seq_num
            out_prefix = "%s/amp%d_var%d" % (tmp_dir, ampIndex, i)
            var_amp_seq1, var_amp_seq2 = get_amp_seq(df_amp, ampIndex, seq_num, mtag_recorded, out_prefix, read_len, var)
            sim_reads1, sim_reads2 = get_sim_reads(var_amp_seq1, var_amp_seq2, out_prefix)
            write_sim_reads(sim_reads1, fout1)
            write_sim_reads(sim_reads2, fout2)

        raw_seq_num = molecule_num - var_seq_num
        out_prefix = "%s/amp%d_wt" % (tmp_dir, ampIndex)
        raw_amp_seq1, raw_amp_seq2 = get_amp_seq(df_amp, ampIndex, raw_seq_num, mtag_recorded, out_prefix, read_len)
        sim_reads1, sim_reads2 = get_sim_reads(raw_amp_seq1, raw_amp_seq2, out_prefix)
        write_sim_reads(sim_reads1, fout1)
        write_sim_reads(sim_reads2, fout2)
    fout2.close()
    fout1.close()


def main():
    # parse arguments
    parse_arguments()

    # read amplicon list
    df_amp = pd.read_csv(amplicon_list, header=0, sep="\t", index_col="#index")

    # read variant file
    df_var = pd.read_csv(variant_table, header=0, sep="\t", index_col="ampIndex")

    # simulate reads
    simulate_reads(df_amp, df_var)


if __name__ == '__main__':
    main()
