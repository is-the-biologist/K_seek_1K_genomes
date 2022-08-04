import pandas as pd
import multiprocessing
import os
import gzip
import time
import sys


def dupliRemoval(cramName):
    """
    Method for removing duplicate read counts from k-seek .total tables. Uses a txt file containing the query names of reads
    and uses them to remove the counts by checking them agains the .rep.gz reads.

    :param cramName:
    :return:
    """

    start = time.time()
    dupe_handles = set( pd.read_csv(f"/fs/cbsuclarkfs1/storage/30x_1K_genomes/CRAM/{cramName}.dupes.txt", header=None)[0].values)

    #iterate through each of the .rep files

    for rep_file in os.listdir("/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/kmer_splits/reads"):
        if rep_file.split("_")[0] == cramName:
            read_set = {} #read in the reads
            with gzip.open(os.path.join("/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/kmer_splits/reads",rep_file), "rt") as myRep:
                n = 0
                for read in myRep:
                    if n == 0:#query name
                        try:
                            header = read.strip("\n").split(" ")[1].split("/")[0]
                        except IndexError:
                            sys.stderr.write(f"{rep_file} malformed header:\n{read}\n")
                            sys.exit()
                        read_set[header] = 0
                    elif n == 4: #tandem counts
                        read_set[header] = read.strip("\n")
                        n = -1

                    n += 1
            #uniq_reads = read_set.keys() - dupe_handles
            #duplicate_reads = read_set.keys() - uniq_reads

            duplicate_reads = read_set.keys() & dupe_handles #get intersection

            #read in the total file and modify it
            total_file = os.path.join("/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/kmer_splits/total", rep_file.replace(".gz", ".total"))
            temp_total = pd.read_csv(total_file, sep="\t", header=None, index_col=0)

            for read in duplicate_reads:
                #remove the counts from the tandem repeat when they are derived from a pcr/optical duplicate
                tandem, count = read_set[read].split("=")
                try:#exception when repeat not found in file -- tandem could be out of register
                    temp_total.loc[tandem] = temp_total.loc[tandem] - int(count)
                except KeyError:
                    sys.stderr.write(f"{tandem} in {rep_file} not found in .total\n")
            temp_total = temp_total.loc[~(temp_total==0).all(axis=1)] #remove zeros
            final_handle = os.path.join("/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/kmer_splits/temp", rep_file.replace(".gz", ".total"))
            temp_total.to_csv(final_handle,sep="\t", header=None)

    end = time.time()

    print(end-start)

if __name__ == '__main__':
    sample_sheet = pd.read_csv('1k_CRAM_sample_sheet.final.txt', header=None)[0].values
    #sample_sheet = ["HG03646"]
    threads = 30
    myPool = multiprocessing.Pool(threads)
    myPool.map(dupliRemoval, sample_sheet)
    myPool.close()