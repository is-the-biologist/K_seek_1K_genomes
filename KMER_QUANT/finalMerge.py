import pandas as pd
import numpy as np

def merge(fh):
    kmer_table = pd.read_csv(fh, sep="\t", index_col=0)
    sample_sheet = pd.read_csv("1k_CRAM_sample_sheet.final.txt", header=None)[0].values
   # sample_sheet = ["HG03646", "HG00353"]
    #assign a hash table with keys as the sample name and values as the index of each subset:

    combined_inds = {sample:[] for sample in sample_sheet}


    for si in range(kmer_table.shape[0]):#iterate through index labels of df and add index number according to sample name
        combined_inds[ kmer_table.index[si].split("_")[0] ].append(si)


    final_sums = []
    for sample in sample_sheet:

        new_df = np.sum(kmer_table.iloc[combined_inds[sample]], axis=0)
        new_df["sample"] = sample
        final_sums.append(new_df)
        print(sample)
    completed_df = pd.concat(final_sums, axis=1).T
    completed_df.to_csv("1k_genomes.final.merged.rep.compiled",sep="\t", index=None)
#merge(fh="temp.comp.df.rep.compiled")
merge(fh="1k_genomes.final.rep.compiled")
