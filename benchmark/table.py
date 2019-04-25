import itertools
import os
from numpy import random
from util import temporary_seed, safe_dir
from biom.util import biom_open

def random_subset(table, otu_size=None, sample_size=None, seed=None):
    otu_ids = table.ids(axis='observation')
    sample_ids = table.ids(axis='sample')

    # temporary_seed helps make sampling thread-safe
    with temporary_seed(seed):
        filters = []
        if otu_size is not None:
            otu_subset = random.choice(otu_ids, size=otu_size, replace=False)
        else:
            otu_subset = otu_ids

        if sample_size is not None:
            sample_subset = random.choice(sample_ids, size=sample_size, replace=False)
        else:
            sample_subset = sample_ids

    return otu_subset, sample_subset

def filter_table(table, otu_subset=None, sample_subset=None):
    new_table = table.copy()
    if otu_subset is not None:
        filter_otus = lambda values, id_, md: \
                id_ in otu_subset
        new_table.filter(filter_otus, axis='observation', inplace=True)
    if sample_subset is not None:
        filter_samples = lambda values, id_, md: \
                id_ in sample_subset
        new_table.filter(filter_samples, axis='sample', inplace=True)
    return new_table

def get_random_subtable(table, otu_size=None, sample_size=None, seed=None):
    """
    Parameters
    ----------

    table (biom.Table): 
        table to subset

    otu_size (int or None):
        number of otus in final table (must be equal or 
        smaller than `table`). All otus are used if None

    sample_size (int or None):
        number of samples in final table (must be equal or 
        smaller than `table`). All otus are used if None

    seed (int):
        seed for the random generator


    Returns
    -------

    biom.Table: 
        subset of original table with `otu_size` otus and 
        `sample_size` samples.

    """
    otu_subset, sample_subset = random_subset(table, otu_size=otu_size, sample_size=sample_size, seed=seed)
    return filter_table(table, otu_subset=otu_subset, sample_subset=sample_subset)

def prep_args(otu_sizes=None, sample_sizes=None, reps=1):
    """
    Generates (otu_size, sample_size, rep_number, seed) pairs
    """
    if otu_sizes is None:
        otu_sizes = [None]
    if sample_sizes is None:
        sample_sizes = [None]

    total_subsets = len(otu_sizes)*len(sample_sizes)*reps
    sizes = total_subsets * reps
    size_combs = list(itertools.product(otu_sizes, sample_sizes))
    size_reps = [(*size_comb, i)  for i in range(reps) for size_comb in size_combs]
    seeds = range(total_subsets)

    sizes_with_seeds = [(*size_reps[seed], seed) for seed in seeds]
    return sizes_with_seeds

arg_names = ['otu_size', 'sample_size', 'rep', 'seed']
def generate_random_tables(table, tree, output_dir, otu_sizes=None, sample_sizes=None, reps=1):

    # TODO use safe_dir to create output_dir if not created
    safe_dir(output_dir)

    args = prep_args(otu_sizes=otu_sizes, sample_sizes=sample_sizes, reps=reps)

    for arg in args:
        # prepare output info
        file_start = 'otu_size-{}--sample_size-{}--'\
                     'rep-{}--seed-{}'.format(*arg)
        full_file_start = os.path.join(output_dir, file_start)

        # get subsetted table
        table_subset = get_random_subtable(table, otu_size=arg[0],
                           sample_size=arg[1], seed=arg[3])
        table_subset.table_id = file_start

        # write out biom table
        with biom_open(full_file_start + '.biom', 'w') as fp:
            table_subset.to_hdf5(fp, "subset: " + file_start)

        # create a sheared tree based off the table
        otu_ids = table_subset.ids('observation')
        tree_subset = tree.shear(otu_ids)
        for node in tree.traverse():
            if node.length is None:
                node.length = 0
        tree_subset.write(full_file_start + '.newick')

    # write out arguments to file
    args_file = os.path.join(output_dir, 'args.txt')

    with open(args_file, 'w') as fp:
        for arg in args:
            for item in arg:
                fp.write(str(item) + '\n')
            fp.write('\n')

    return output_dir

if __name__ == '__main__':
    from biom import load_table
    from skbio import TreeNode
    table = load_table('../../trial-unifrac/moving-pictures-downloaded/table-dir/feature-table.biom')
    tree = TreeNode.read('../../trial-unifrac/moving-pictures-downloaded/tree-dir/tree.nwk')

    generate_random_tables(table, tree, '../data/demo2', otu_sizes=[1,10], sample_sizes=[1,10], reps=3)
