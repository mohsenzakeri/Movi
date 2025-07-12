import sys
import pandas as pd

ROOT_TAXON = 1
UNDEFINED_TAXON = 0

def find_parent_node(taxon, parents, taxon_to_index):
    if taxon not in taxon_to_index:
        return ROOT_TAXON
    if taxon_to_index[taxon] not in parents:
        return ROOT_TAXON
    return parents[taxon_to_index[taxon]]

def find_LCA(taxon1, taxon2, parents, taxon_to_index):
    if taxon2 == UNDEFINED_TAXON:
        return taxon1
    if taxon1 == taxon2:
        return taxon1

    taxon1_lineage = {}
    while taxon1 != ROOT_TAXON:
        taxon1_lineage[taxon1] = 1
        taxon1 = find_parent_node(taxon1, parents, taxon_to_index)

    while taxon2 != ROOT_TAXON:
        if taxon2 in taxon1_lineage:
            return taxon2
        taxon2 = find_parent_node(taxon2, parents, taxon_to_index)

    return ROOT_TAXON

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Error: Please provide the tree (nodes.dmp) file and the Movi Color results")
        print("Usage: python find_Pseudomondatoa_taxons.py <tree file> <movi_color_results>")
        sys.exit(1)

    tree_file = sys.argv[1]
    movi_color_results_file = sys.argv[2]

    tree_df = pd.read_csv(tree_file,  sep="\t", header=None)
    children, parents = tree_df[0], tree_df[2]

    taxon_to_index = dict(zip(children, range(len(children))))

    movi_results = pd.read_csv(movi_color_results_file, names = ["read_name", "doc", "doc_secondary"])

    movi_results['lca'] = [
        find_LCA(a, b, parents, taxon_to_index)
        for a, b in zip(movi_results['doc'], movi_results['doc_secondary'])
    ]

    movi_results.to_csv(movi_color_results_file, index=False)