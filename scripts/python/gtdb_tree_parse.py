#!/usr/bin/python

from ete3 import Tree
import pandas as pd

bac120_r95_tree = Tree('db/gtdb_release202/bac120_r202.tree', format=1, quoted_node_names=True)

for node in bac120_r95_tree.traverse():
    if node.is_leaf() or node.is_root(): continue
    
    node.support = float(node.name.split(':', 1)[0]) / 100
    node.name = (node.name.split(':', 1)[1]).replace('; ', ';') if ':' in node.name else ''

bac120_r95_tree.write(outfile='raw/gtdb/bac120_r202_ete3fmt1.tree', format=1)
bac120_r95_tree.write(outfile='raw/gtdb/bac120_r202_ete3fmt5.tree', format=5)


ar122_r95_tree = Tree('db/gtdb_release202/ar122_r202.tree', format=1, quoted_node_names=True)

for node in ar122_r95_tree.traverse():
    if node.is_leaf() or node.is_root(): continue
    
    node.support = float(node.name.split(':', 1)[0]) / 100
    node.name = (node.name.split(':', 1)[1]).replace('; ', ';') if ':' in node.name else ''

ar122_r95_tree.write(outfile='raw/gtdb/ar122_r202_ete3fmt1.tree', format=1)
ar122_r95_tree.write(outfile='raw/gtdb/ar122_r202_ete3fmt5.tree', format=5)