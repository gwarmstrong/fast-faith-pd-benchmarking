from table import generate_random_tables
from biom import load_table
from skbio import TreeNode
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--tree', '-t', dest='tree',
                    help='File name for tree')
parser.add_argument('--table' '-T', dest='table',
                    help='File name for table')
parser.add_argument('--otu-sizes', '-o', dest='otu_sizes',
                    nargs='*', type=int, default=None,
                    help='List of number of OTUS to '\
                                'sample.')
parser.add_argument('--sample-sizes', '-s', dest='sample_sizes',
                    nargs='*', type=int, default=None,
                    help='List of number of samples to '\
                                'sample.')
parser.add_argument('--reps', '-r', dest='reps', default=1,
                    type=int,
                    help='Number of times to sample each')
parser.add_argument('--output-dir', '-d', dest='output_dir',
                    help="directory to write results to")

args = parser.parse_args()

table = load_table(args.table)
tree = TreeNode.read(args.tree)
generate_random_tables(table, tree, args.output_dir, 
                       otu_sizes=args.otu_sizes,
                       sample_sizes=args.sample_sizes,
                       reps=args.reps)

# '../../trial-unifrac/moving-pictures-downloaded/table-dir/feature-table.biom'
# '../../trial-unifrac/moving-pictures-downloaded/tree-dir/tree.nwk'
