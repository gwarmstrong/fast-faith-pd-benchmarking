from numpy import random
from util import temporary_seed

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
    otu_subset, sample_subset = random_subset(table, otu_size=otu_size, sample_size=sample_size, seed=seed)
    return filter_table(table, otu_subset=otu_subset, sample_subset=sample_subset)
