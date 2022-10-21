"""
Program to calculate degree of interspersion between tandem k-mer sequences


"""
import os
import gzip
import sys
from multiprocessing import Pool
import numpy as np
import pandas as pd

class Quantify:
    def __init__(self):
        self.tandem_offset_table = {}
        self.sorted_tandems = []

    def extractReads(self, tandems, genome_list):
        self.sorted_tandems = tandems
        for tandem_seq in tandems:
            f, r = tandem_seq.split("/")

            for k in range(len(f)):
                offset_f = f[k + 1::] + f[:k + 1]
                offset_r = r[k + 1::] + r[:k + 1]

                #create a hash table where each offset is a key and the cannonical sequence is the value for easy look-up
                #this is probably fine for small table, but may be too memory intensive
                self.tandem_offset_table[offset_f] = tandem_seq
                self.tandem_offset_table[offset_r] = tandem_seq
        #this will probably be multithreaded
        myPool = Pool(processes=30)
        myPool.map(self.findInterspersed, genome_list)


    def findInterspersed(self, genome):
        directory = '/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/reads'
        #directory = "/Users/iskander/Documents/Barbash_lab/HumanGenome_TE/interspersed/dummy_files"
        #generate a table where we store the sequence headers and the tandems they hold:
        seq_header_storage = {}

        ###

        #stores the # of links between tandems btwn PE fastqs
        interspersion_matrix = np.full(fill_value=0, shape=(len(self.sorted_tandems), len(self.sorted_tandems)))

        #all_offsets = set(self.tandem_offset_table.keys())

        #iter thru all the files and store the sequence headers for each tandem if in our set of tandems we care about

        for file in os.listdir(directory):
            if file.startswith(genome):
                #iterate thru the fastq files:
                with gzip.open(os.path.join(directory, file), mode='rt') as fq:

                    line_num = 0
                    for line in fq:
                        if line_num == 4:
                            kmer, count = line.split('=')

                            #if (kmer in all_offsets) == True: #set operation to check O(1)
                            try:#try/except is faster
                                full_tandem = self.tandem_offset_table[kmer]
                            except KeyError:
                            #else:
                                full_tandem = "F" #when tandem is not in keys we want to ignore

                            #headers = set(seq_header_storage.keys())


                            try:#if (seq_header in headers) == True:# check if header is in set of headers quickly
                            #try/except should be faster
                                if full_tandem != "F":

                                    tandem_pe = seq_header_storage[seq_header] #grab tandem from other paired end read

                                    #add +1 link to the tandems on the linker matrix
                                    interspersion_matrix[self.sorted_tandems == tandem_pe, self.sorted_tandems == full_tandem] += 1
                                    interspersion_matrix[self.sorted_tandems == full_tandem, self.sorted_tandems == tandem_pe] += 1
                                    seq_header_storage.pop(seq_header)
                            except KeyError:#else:

                                if full_tandem != "F":
                                    seq_header_storage[seq_header] = full_tandem
                            line_num = 0

                        elif line_num == 0: #grab sequence header
                            seq_header = line.split("/")[0]
                            line_num += 1
                        else:
                            line_num += 1
            else:
                pass

        np.save(f'/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/interspersion/counts/{genome}.interspersed.counts.npy', interspersion_matrix)
    def combine(self, dir='/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/interspersion/counts'):
        sumMatrix = np.zeros(shape=(125,125))
        for f in os.listdir(dir):
            intr = np.load(os.path.join(dir, f))
            sumMatrix = sumMatrix + intr

        np.save('/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/interspersion/1K_genomes.kmer.interspersion.npy', sumMatrix)


if __name__ == '__main__':

    tandem_list = pd.read_csv('2022_08_15_tandem_list.txt', header=None)[0].values
    genome_list = pd.read_csv('1k_CRAM_sample_sheet.final.txt', header=None)[0].values

    myQ = Quantify()
    myQ.extractReads(tandem_list, genome_list)
    myQ.combine()
    #