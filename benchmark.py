from skbio.diversity import alpha_diversity
from skbio.tree import TreeNode
import numpy as np
import time
from biom import load_table

def run_faith(counts, otu_ids, tree):
    return alpha_diversity('faith_pd', counts, tree=tree, otu_ids=otu_ids)

def run_fast_faith(counts, otu_ids, tree, shear=False):
    return alpha_diversity('fast_faith_pd', counts, tree=tree, otu_ids=otu_ids, shear=shear)

if __name__ == '__main__':

    #with open('data/moving_pictures_tree.nwk') as fp:
    #    tree_line = fp.readline().rstrip()
    tree = TreeNode.read('data/gg_13_8_otus/trees/97_otus.tree')

    tips = list(tree.tips())
    tree_otus = {tip.name for tip in tips}
    
    bt = load_table('data/moving-pictures/67_otu_table.biom')

    counts = np.asarray([bt.data(i) for i in bt.ids()])
    otu_ids = bt.ids('observation')

    #with open('data/moving_pictures_table.txt', 'r') as fp:
    #    lines = fp.readlines()

    #lines = [line.rstrip().split('\t') for line in lines]

    ## only keep otu's that we have in our tree
    #otu_ids = [line[0] for line in lines if line[0] in tree_otus]
    ##otu_ids = otu_ids[1:]
    #counts = [[entry for entry in line[1:]] for line in lines \
    #        if line[0] in tree_otus]
    #counts = np.array(counts, dtype=np.float32).T
    #counts = counts.T

    print("Adding zero branch lengths")
    for node in tree.traverse():
        if node.length is None:
            node.length = 0

    print("Starting...")
    t0 = time.time()
    f1 = run_faith(counts, otu_ids, tree)
    t1 = time.time()
    print("Faith():          {}".format(t1-t0))
    t2 = time.time()
    f2 = run_fast_faith(counts, otu_ids, tree)
    t3 = time.time()
    print("FastFaith():      {}".format(t3 - t2))
    t4 = time.time()
    f3 = run_fast_faith(counts, otu_ids, tree, shear=True)
    t5 = time.time()
    print("FastFaithShear(): {}".format(t5 - t4))
    if np.allclose(f1.values, f2.values):
        print("Success!")
    else:
        print(":(")

