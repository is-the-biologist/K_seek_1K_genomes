import os.path
import sys
import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering
from sklearn.metrics import rand_score, adjusted_rand_score
from scipy.stats import chi2
from operator import itemgetter
import networkx as nx
from community import community_louvain
from multiprocessing import Pool
from statsmodels.stats.multitest import fdrcorrection
from itertools import combinations


class ClusterTandems:
    """
    This module implements several different methods for clustering although the input data is the same:

    Pre-processing:
    The k-seek .rep files are used to generate an N x N symmetrical matrix, where each cell represent the number of reads
    that have monomer A in read 1 and monomer B in read 2 from the paired-end reads. We compute this for each sample
    individually and also sum these values across all samples to get values for the complete data set. These matrices will
    be the primary data used to generate the graphs for the spectral clustering.

    1. Divisive Spectral:
        This module uses a "divisive hierarchical spectral clustering" approach to cluster together tandem monomers that form larger,
        more complex higher order repeats. This approach is highly inspired by Altemose, 2014 where a similar strategy was
        employed on Hsat2,3 monomers to detect higher order repeats. The novel feature of this being our use of affinity matrix
        and diversity of monomers.



        Step 1.
        Use the matrices of the paired-end reads to generate and Observed / Expected probability of finding two monomers on the
        same read fragment as P(A&B) / P(A)*P(B). We calculate the P(A&B) as the frequency of the 'chimeric' read with A and B.
        Then P(A) as frequency of all reads with A monomer and P(B) as frequency of all reads with monomer B. Thus generating
        values that represent how more likely two monomers are to be found on the same fragment relative to random chance.

        Step 2.
        Use a divisive spectral clustering approach to generate clusters of monomers that are often found in the same fragment.
        This approach uses the Obsv/Exp probability matrices from step 1 as input to a Spectral Clustering algorithm. We use
        this matrix as an affinity matrix from which to construct the graph for spectral clustering where nodes represent
        monomers and edges represent their relative probability of being found in the same fragment, such that edges in the
        graph that connect nodes that are more likely to be found in the same fragment by random chance will have higher weights
        (O/E > 1). We perform recursive clustering wherein we perform sequential binary splits on the graph until the edges of
        each cluster have a median O/E > 1. This means that in each cluster the median enrichment of monomers being found in the
        same fragment is positive. The result is a collection of monomers that on the whole are more likely to be found together
        in the same sets of fragments than by chance, suggesting co-localization and higher order structures.

    2. Louvain:
        Louvain clustering uses modularity optimization on, in this case, a weighted network to create a set of
        communities that partition the graph optimally. In this case the edge-weights are the E/O ratios and the
        communities are clusters that represent tandems that are found to form higher-order interspersed structures
        more often.

    """
    def __init__(self):
        pass

    def upper_tri_indexing(self, A):
        m = A.shape[0]
        r, c = np.triu_indices(m, 1)
        return A[r, c]

    def recurClustering(self, affinity_matrix):

        # function to perform the recursion:
        clustering = SpectralClustering(n_clusters=2, assign_labels='discretize', affinity='precomputed',
                                        random_state=0).fit_predict(affinity_matrix)
        NTs, cluster_labels = zip(*sorted(zip(affinity_matrix.columns, clustering), key=itemgetter(1)))
        NTs = np.asarray(NTs)
        cluster_labels = np.asarray(cluster_labels)

        for c in range(2):

            new_matrix = affinity_matrix.loc[NTs[cluster_labels == c], NTs[cluster_labels == c]]

            if new_matrix.shape[0] > 1: #exit condition for singleton clusters:
                med_enrichment = np.min( #we can use min value here to ensure all links are > 1
                    self.upper_tri_indexing(new_matrix.values))

                if med_enrichment <= 1:
                    if new_matrix.shape[0] == 2: #to prevent warnings in the spectral clustering we can just manually split the cluster here
                        for c in NTs[cluster_labels == c]:
                            print(c)
                    else:
                        self.recurClustering(new_matrix)


                else: #guarantee exit condition
                    cluster = NTs[cluster_labels == c]
                    print(",".join(cluster))

            else:
                cluster = NTs[cluster_labels == c]
                print(",".join(cluster))

    def louv(self, affinity_matrix, pval_matrix, sample, correction_type='bh'):

        # simple correction: if p-value > alpha then we set value of odds to 1
        #logic is simple H0 fails to reject therefore b=c
        if correction_type == 'bonf':
            #bonferoni correction:
            alpha = 0.05 / len(np.triu_indices(pval_matrix.shape[0])[0])
            affinity_matrix[pval_matrix > alpha] = 1

        elif correction_type == 'bh':
            #benjamini-hochberg
            alpha_pass, q = fdrcorrection(pval_matrix.values[np.triu_indices(pval_matrix.shape[0])], alpha=0.05)
            # fdr corrected  p-values:
            Q = np.zeros(shape=pval_matrix.shape)
            Q[np.triu_indices(pval_matrix.shape[0])] = q
            Q = np.where(Q, Q, Q.T)
            Q = pd.DataFrame(data=Q, columns=pval_matrix.columns, index=pval_matrix.columns)
            affinity_matrix[Q > 0.05] = 1
        else:
            sys.stderr.write("Incorrect correction type selected. Try: 'bh' or 'bonf'\n")
            sys.exit()

        partition = community_louvain.best_partition(nx.from_pandas_adjacency(affinity_matrix), random_state=42)
        NTs, cluster_labels = zip(*sorted(list(partition.items()), key=itemgetter(1)))

        cluster_df = pd.DataFrame(data={'tandems': NTs, 'label': cluster_labels})
        output_name = os.path.split(sample)[-1].split('.')[0]
        cluster_df.to_csv(f'/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/interspersion/LOUVAIN_CLUSTERS_2022_08_17/{output_name}.louvain.clusters.csv', index=None)

    def calcIndependence(self, interspersed, regularize=False):
        #calculates the Observed / Expected ratio of the number of tandem monomers
        interspersed[np.diag_indices(interspersed.shape[0])] = np.diagonal(interspersed) / 2

        #remove all tandems absent from the dataset:
        #pruned_matrix = interspersed[np.sum(interspersed, axis=0) > 0][:, np.sum(interspersed, axis=0) > 0]

        #necessary for the rare tandems
        tandem_list = np.asarray([k.split("/")[0] for k in pd.read_csv('/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/scripts/2022_08_15_tandem_list.txt', header=None)[0]])

        JS = np.zeros(shape=interspersed.shape)
        PVAL = np.zeros(shape=interspersed.shape)
        total_links = np.sum(interspersed, axis=0)
        for i in range(len(total_links)):
            for j in range(len(total_links)):
                # P(A&B) / P(A)*P(B)
                O = (interspersed[i, j] / (np.sum(total_links)))
                E = (total_links[j] / np.sum(total_links)) * (total_links[i] / np.sum(total_links))
                if E == 0:# handle P(A)=0 or P(B)=0
                    JS[i, j] = 0
                else:
                    JS[i, j] = O / E

                #calculate the p_val via mcnemar test:
                #chi2 = (b-c)**2 / (b+c)
                b = O*np.sum(total_links)
                c = E*np.sum(total_links)
                chi_test = (b-c)**2 / (b+c)
                pval = chi2.sf(chi_test, df=1)
                PVAL[i, j] = pval

        #regularization alternative method is to fully connect the graph by adding a pseudocount to the affinity matrix
        #we can try replacing zeros with the minimum float in python:
        if regularize:
            JS[JS == 0] = sys.float_info.min

        JS_df = pd.DataFrame(data=JS, columns=tandem_list, index=tandem_list)

        PVAL_df = pd.DataFrame(data=PVAL, columns=tandem_list, index=tandem_list)

        return JS_df, PVAL_df

    def invokeCluster(self, sample, clustering_algo='louvain'):

        #load sample
        link_matrix = np.load(sample)
        #calculate indepdence as P(A&B) / P(A)*P(B)
        independence_prob, pval = self.calcIndependence(link_matrix, regularize=False)

        if clustering_algo == 'spectral':
            #compute clusters and print
            self.recurClustering(affinity_matrix=independence_prob)
        elif clustering_algo == 'louvain':
            self.louv(affinity_matrix=independence_prob, sample=sample, pval_matrix=pval)
        else:
            sys.stderr.write('Incorrect argument of clustering type, try: "spectral" or "louvain"\n')

class ClusterAnalysis:

    """

    This module contains the information for performing analysis on the clusters of tandem interspersed structures from
    the 1k genomes dataset.

    Functions to calculate Rand Index and Adjusted Rand Index as pairwise comparisons between individuals.

    """

    def __init__(self):
        self.tandem_list = np.asarray([k.split("/")[0] for k in pd.read_csv('/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/scripts/2022_08_15_tandem_list.txt', header=None)[0]])
        self.sample_names = pd.read_csv("/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/scripts/1k_CRAM_sample_sheet.final.txt", header=None)[0]

    def readCluster(self, sample):
        """
        Function to read in the cluster labels

        :param sample:
        :return:
        """
        fname = self.sample_names[sample]
        cluster_df = pd.read_csv(f"/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/interspersion/LOUVAIN_CLUSTERS_2022_08_17/{fname}.louvain.clusters.csv")

        order_fnc = lambda x: np.where(self.tandem_list == x[0])[0]
        cluster_labels = np.asarray(sorted(cluster_df.values, key=order_fnc))[:,1]

        return cluster_labels

    def computeScore(self, pair):
        #function wrapping the score computation
        A = self.readCluster(pair[0])
        B = self.readCluster(pair[1])

        score = rand_score(A, B)

        return score

    def pairwiseComparisons(self, threads=30):
        """
        Function to perform all PW comparisons using itertools and multiprocessing for speed.

        :param threads:
        :return:
        """

        sample_index = [i for i in range(2504)]
        RI_matrix = np.full(fill_value=1.0, shape=(len(sample_index), len(sample_index)))
        pw_combos = list(combinations(sample_index, 2))

        myPool = Pool(threads)
        scores = myPool.map(self.computeScore, pw_combos)
        for sc in range(len(scores)):
            i = pw_combos[sc][0]
            j = pw_combos[sc][1]
            RI_matrix[i,j] = scores[sc]
            RI_matrix[j,i] = scores[sc]

        np.save('POPULATION.LOUVAIN.RAND_SCORE.npy', RI_matrix)

    def estCN(self, counts_matrix):
        """
        Function to estimate the copy-number from read-depth of the various k-mer clusters. We will just use the counts
        matrices from the PE data as the read depth of each junction between two kmers. We can attempt to correct for
        GC bias by using the combined GC of both k-mers?

        :return:
        """

        #get the GC bias bin info
        f = counts_matrix.split(".")[0]
        depth_summary = pd.read_csv( os.path.join("/fs/cbsuclarkfs1/storage/30x_1K_genomes/gc/GC_EST_UNIQ_2022_05_18", f+".gc.depth.tsv"), sep="\t" , header=None)
        # calculate depth | %GC using bins
        med_value = np.median(depth_summary[3].values)
        # adding conditional argument such that any bins that had < 10000 reads compose them are to be discarded and instead use the median value
        # this is because the edges of the estimates become spikey and unstable as the # of reads that make up a bin
        # is quite small
        gc_depth = {d[0]: (d[3] if d[1] >= 10000 else med_value) for d in
                    depth_summary.values}  # assign bins to dict
        counts = np.load(os.path.join("/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/interspersion/intrsp_counts",counts_matrix) )

        gc_bins = np.asarray([((b + 1) / 100) for b in range(100)])
        get_gc_bins = lambda gc: min(gc_bins[float(gc) <= gc_bins])
        compute_gc = lambda kmer: ( np.sum(np.asarray(list(kmer)) == "G") + np.sum(np.asarray(list(kmer)) == "C") )/ len(kmer)

        correction_factor = np.full(fill_value=np.nan, shape=counts.shape) #collect all correction factors first
        for i in range(counts.shape[0]):
            for j in range(counts.shape[0]):

                #we are going to concatenate the two k-mers as a rough estimate of what the GC bias would be at the junction
                concat = self.tandem_list[i] + self.tandem_list[j]
                gc = compute_gc(concat) #compute the GC of concat
                gc_bin = get_gc_bins(gc) #get gc bin
                correction_factor[i,j] = gc_depth[gc_bin]

        corrected_counts = counts / correction_factor

        np.save(f"/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/interspersion/GC_corrected_counts/{f}.GC.corrected.interspersed.counts.npy", corrected_counts)

    def GC_correction(self, threads):
        myPool = Pool(threads)
        matrix_names = [k+".interspersed.counts.npy" for k in self.sample_names.values]
        myPool.map(self.estCN, matrix_names)

    def computeCN_means(self):
        """
        Compute the mean copy-number of each interspersion junction using the GC corrected read-depth.
        :return:
        """
        mean_CN = np.zeros(shape=(125,125)).astype(float)
        for sample in self.sample_names.values:
            A = np.load(os.path.join("/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/interspersion/GC_corrected_counts", sample+".GC.corrected.interspersed.counts.npy"))
            mean_CN = mean_CN + A

        mean_CN = mean_CN / 2504

        np.save('1k_genomes.meanCN.global.interspersion.npy', mean_CN)

if __name__ == '__main__':

    #cluster
    harch = ClusterTandems()
    sample_list = ["/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/interspersion/counts/"+sample+".interspersed.counts.npy" for sample in pd.read_csv('/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/scripts/1k_CRAM_sample_sheet.final.txt', header=None)[0].values]
    myPool = Pool(30)
    myPool.map(harch.invokeCluster, sample_list)

    #analyze patterns
    cAnal = ClusterAnalysis()
    cAnal.pairwiseComparisons()


