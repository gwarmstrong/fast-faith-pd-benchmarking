from skbio.diversity import alpha_diversity
from skbio.tree import TreeNode
import numpy as np
import time

def run_faith(counts, otu_ids, tree):
    alpha_diversity('faith_pd', counts, tree=tree, otu_ids=otu_ids)

def run_fast_faith(counts, otu_ids, tree):
    alpha_diversity('fast_faith_pd', counts, tree=tree, otu_ids=otu_ids)

if __name__ == '__main__':

    with open('data/moving_pictures_tree.nwk') as fp:
        tree_line = fp.readline().rstrip()
    tree = TreeNode.read([tree_line])

    tips = list(tree.tips())
    tree_otus = {tip.name for tip in tips}
    
    with open('data/moving_pictures_table.txt', 'r') as fp:
        lines = fp.readlines()

    lines = [line.rstrip().split('\t') for line in lines]

    # only keep otu's that we have in our tree
    otu_ids = [line[0] for line in lines if line[0] in tree_otus]
    #otu_ids = otu_ids[1:]
    counts = [[entry for entry in line[1:]] for line in lines \
            if line[0] in tree_otus]
    counts = np.array(counts, dtype=np.float32).T


    t0 = time.time()
    run_faith(counts, otu_ids, tree)
    t1 = time.time()
    print("Faith(): {}".format(t1-t0))
    t2 = time.time()
    run_fast_faith(counts, otu_ids, tree)
    t3 = time.time()
    print("FastFaith(): {}".format(t3 - t2))


