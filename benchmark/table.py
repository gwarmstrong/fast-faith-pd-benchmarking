import itertools
import os
from numpy import random
from util import temporary_seed, safe_dir
import biom
from biom.util import biom_open
import bp
from typing import Optional


def random_subset(table, otu_size=None, sample_size=None, seed=None):
    """Obtains random otu and sample ids

    Parameters
    ----------

    table : biom.Table
        table to be subset

    otu_size : Optional[int]
        Number of otu's to be in subset. If None, includes all

    sample_size : Optional[int]
        Number of samples ot be in subset. If None, inludes all

    seed : Optional[int]
        Seed for random number generator

    Returns
    -------

    (list[str], list[str])
        The (list of otu ids, list of sample ids) for this subset of the
         table

    """

    otu_ids = table.ids(axis='observation')
    sample_ids = table.ids(axis='sample')

    # temporary_seed helps make sampling thread-safe
    with temporary_seed(seed):
        if otu_size is not None:
            otu_subset = random.choice(otu_ids, size=otu_size, replace=False)
        else:
            otu_subset = otu_ids

        if sample_size is not None:
            sample_subset = random.choice(sample_ids, size=sample_size,
                                          replace=False)
        else:
            sample_subset = sample_ids

    return otu_subset, sample_subset


def filter_table(table, otu_subset=None, sample_subset=None):
    """Filters a table given otu and sample ids

    Parameters
    ----------

    table : biom.Table
        Table to be filtered

    otu_subset : Optional[list[str]]
        List of otu's to keep in table. If None, keep all otu's.

    sample_subset : Optional[list[str]]
        List of samples to keep in table. If None keep all samples.

    Returns
    -------

    biom.Table
        A filtered version of the original table, based on otu and
         sample subsets


    """
    # create a new table then filter this table in place
    new_table = table.copy()
    if otu_subset is not None:
        new_table.filter(ids_to_keep=otu_subset, axis='observation',
                         inplace=True)
    if sample_subset is not None:
        new_table.filter(ids_to_keep=sample_subset, axis='sample',
                         inplace=True)
    return new_table


def get_random_subtable(table: biom.Table,
                        otu_size: Optional[int],
                        sample_size: Optional[int],
                        seed: Optional[int]):
    """Returns a random sub-table with specified otu and sample sizes


    Parameters
    ----------

    table : biom.Table
        table to subset

    otu_size : int, optional
        number of otus in final table (must be equal or 
        smaller than `table`). All otus are used if None.

    sample_size : int, optional
        number of samples in final table (must be equal or 
        smaller than `table`). All otus are used if None.

    seed : int, optional
        seed for the random generator


    Returns
    -------

    biom.Table
        subset of original table with `otu_size` otus and 
         `sample_size` samples.

    """
    # get ids to for subsets
    otu_subset, sample_subset = random_subset(table,
                                              otu_size=otu_size,
                                              sample_size=sample_size,
                                              seed=seed)
    # filter and return table based on otu_subset and sample_subset
    return filter_table(table,
                        otu_subset=otu_subset,
                        sample_subset=sample_subset)


def prep_args(otu_sizes=None, sample_sizes=None, reps=1):
    """Generates (otu_size, sample_size, rep_number, seed) pairs

    Parameters
    ----------

    otu_sizes : Optional[list[int]]
        Specifies the sizes for OTU sets that should be chosen. If None,
         uses all OTUs for all tables.

    sample_sizes : Optional[list[int]]
        Specifies the sizes for sample sets taht should be chose. If
         None, uses all samples for all tables.

    reps : Optional[int]
        Specifies the number of repetitions that each `otu_size` x
         `sample_size` pair should be used to generate a new table.


    Returns
    -------

    list[tuple]
        Each tuple contains (otu_size, sample_size, rep_number, seed)
         for a single table subset

    """
    if otu_sizes is None:
        otu_sizes = [None]
    if sample_sizes is None:
        sample_sizes = [None]

    total_subsets = len(otu_sizes)*len(sample_sizes)*reps
    size_combs = list(itertools.product(otu_sizes, sample_sizes))

    # TODO can probably be done better using zip
    size_reps = [(*size_comb, i) for i in range(reps) for
                 size_comb in size_combs]
    seeds = range(total_subsets)
    sizes_with_seeds = [(*size_reps[seed], seed) for seed in seeds]

    return sizes_with_seeds


arg_names = ['otu_size', 'sample_size', 'rep', 'seed']


def subset_and_write_table_tree(otu_size: int, sample_size: int, rep: int,
                                seed: int, table: biom.Table,
                                tree: bp.BP,
                                output_dir) -> None:
    """Given parameters for a single subset, filter the table and tree
    and write to file

    Parameters
    ----------

    otu_size : int
        Number of OTUs for this specific table

    sample_size : int
        Number of samples for this specific table

    rep : int
        The number indexing the repetition for this table

    seed : int
        The random seed to use when subsetting the table

    table : biom.Table
        The table to subset

    tree : bp.BP
        The tree to be subset for the OTUs in `table`

    output_dir : str
        Location to write the files in


    Returns
    -------

    None

    """
    # prepare output info
    file_start = 'otu_size-{}--sample_size-{}--rep-{}--seed-{}' \
        .format(otu_size, sample_size, rep, seed)
    full_file_start = os.path.join(output_dir, file_start)

    # get subset of table
    table_subset = get_random_subtable(table,
                                       otu_size=otu_size,
                                       sample_size=sample_size,
                                       seed=seed)
    table_subset.table_id = file_start

    # write out biom table
    with biom_open(full_file_start + '.biom', 'w') as fp:
        table_subset.to_hdf5(fp, "subset: " + file_start)

    # create a sheared tree based off the table
    otu_ids = table_subset.ids('observation')

    sheared_tree = tree.shear(set(otu_ids))

    tree_subset = bp.to_skbio_treenode(sheared_tree)

    for node in tree_subset.traverse():
        if node.length is None:
            node.length = 0
    tree_subset.write(full_file_start + '.newick')


def generate_random_tables(table, tree, output_dir, otu_sizes=None,
                           sample_sizes=None, reps=1, job=None):
    """Subset `table` randomly along axes based on {otu,sample}_sizes
    and make the corresponding cuts to tree. Write results to file.

    Parameters
    ----------

    table : biom.Table
        Table to be subset

    tree : skbio.TreeNode
        Tree corresponding to `table`

    output_dir : str
        Directory to write the resulting tables and trees to

    otu_sizes : Optional[list[int]]
        List of numbers of otus to extract in subsets. If None, then all
         OTU's are used

    sample_sizes : Optional[list[int]]
        List of numbers of samples to extract in subsets. If None, then
        all samples are used

    reps : Optional[int]
        Number of tables to generate for each `otu` x `sample` pair. If
         None, do one rep for each table.

    job : Optional[int]
        Index of job. Used to avoid writing args file multiple times. If
         None, the args file is written.


    Returns
    -------

    str
        Directory subsetted trees and tables are written in

    """
    # use safe_dir to create output_dir if not created
    safe_dir(output_dir)

    args = prep_args(otu_sizes=otu_sizes, sample_sizes=sample_sizes, reps=reps)

    # perform table subsetting for all args
    if job is None:  
            for arg in args:
                subset_and_write_table_tree(*arg, table, tree, output_dir)

    # perform table subsetting for specific arg
    else:
        arg = args[job-1]
        subset_and_write_table_tree(*arg, table, tree, output_dir)

    # write file that keeps track of args
    if (job == 1) or (job is None):
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
    table = load_table('../../trial-unifrac/moving-pictures-downloaded'
                       '/table-dir/feature-table.biom')
    tree = TreeNode.read('../../trial-unifrac/moving-pictures-downloaded'
                         '/tree-dir/tree.nwk')

    generate_random_tables(table, tree, '../data/demo2', otu_sizes=[1, 10],
                           sample_sizes=[1, 10], reps=3)
