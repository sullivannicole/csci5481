# ------------
# Imports
# ------------

# !pip install fastaparser
# !pip install pandas==1.5.0
# !pip install scikit-bio
import fastaparser
import pandas as pd
import numpy as np
import re
import itertools
from skbio.tree import TreeNode as tn
import argparse

def make_argparser():
    '''
    Parameters
    -----------
        Technically, this function takes no parameters. Note, though, that the parser expects 1 string.

    Return
    --------
        This function returns an argparse parser, which is then used as the input for the nei_saitou() function. See the docs for more info on how the argparse methods work: https://docs.python.org/3/library/argparse.html
    '''
    parser = argparse.ArgumentParser(prog = 'neighbor_joining.py', 
                                     description = 'A program to run the Nei-Saitou neighbor joining algorithm on an input fna file')
    
    parser.add_argument('input_fna_filepath', type=str, help = 'Path to fasta file, including .fna extension')
    return parser
    

def parse_header(header: list) -> dict:
    """
    Parameters
    -----------
        header: a list produced of sequence headers produced by the fastaparser library reading a .fna file

    Return
    --------
        A dictionary containing a key-value pair for each of the items in the header list
    """

    header_split = re.findall("\[(.*?)\]", header)
    header_ls = [{x.split("=")[0]: x.split("=")[1]} for x in header_split]

    return {k: v for d in header_ls for k, v in d.items()}


def read_seqs(fna_file: str, 
              seq_only: bool = True) -> pd.DataFrame:
    """
    Parameters
    -----------
    fna_file: str
        The path to the fna file to read in, including the '.fna' extension. E.g. 'pfizer_mrna.fna' for a local fasta file containing Pfizer mRNA sequence.
    seq_only: Boolean
        Whether or not to return just the sequences contained in the file, or to include all the metadata contained in the fasta file (headers, tags, locii, etc.). Set to True by default, since it's often useful to extract just the sequence.

    Return
    --------
    pd.DataFrame 
        A dataframe with one row per sequence in the input fasta file. At the very least (if seq_only = True), the dataframe will contain the columns "sequence", "lcl" (the sequence headers), and "specimen", which is the file name, excluding the ".fna" extension. If seq_only = False, the output dataframe will contain all additional columns that can be extracted from the input fasta file. These additional columns vary depending on what's contained in the input, but can include, tags, locii of the start/end of the sequence within the entire genome, and other useful metadata.
    """
    with open(fna_file) as fasta_file:
        rparser = fastaparser.Reader(fasta_file)
        seqs = [seq.sequence_as_string() for seq in rparser]

        if seq_only:
            return seqs

        qparser = fastaparser.Reader(fasta_file, parse_method="quick")
        headers = [seq.header for seq in qparser]

    seq_df = pd.DataFrame(
        {
            "sequence": seqs,
            "lcl": [x.split("[")[0] for x in headers],
            "specimen": fna_file[:-4],
        }
    )

    headers_df = pd.DataFrame(list(map(parse_header, headers)))
    genome_df = pd.concat([seq_df, headers_df], axis=1)

    return genome_df

def genetic_distance_matrix(msa):
    """
    Calculates the pairwise genetic distance between all sequences in a multiple sequence alignment.

    Parameters
    ----------
    msa : List[str]
        The multiple sequence alignment for a list of sequences.

    Returns
    -------
    np.ndarray
        A distance matrix containing the genetic distance between all pairs of sequences in the MSA.
    """
    # Get all pairwise combos in MSA (assume each seq is unique)
    seq_pairs = list(itertools.combinations(msa, 2))
    pair_coords = list(itertools.combinations(range(len(msa)), 2))

    # Next, get % dissimilarity for each pair set
    pairwise_distances = list(map(lambda z: 1 - (sum(x == y for x, y in zip(z[0], z[1])))/len(msa[0]), seq_pairs))
    distance_matrix = np.zeros((len(msa), len(msa)))

    # Assign distances to coord in distance matrix
    for i, ele in enumerate(pair_coords):
        distance_matrix[ele] = pairwise_distances[i]

    # Since the matrix is symmetrical, just transpose a mirror image of the entries above the diag to below the diag
    distance_matrix = distance_matrix + distance_matrix.T - np.diag(distance_matrix.diagonal())

    return distance_matrix

def format_distance_matrix(distance_matrix, tip_labels):
    """
    Formats a distance matrix so it can be used with other tools.

    Parameters
    ----------
    distance_matrix : np.ndarray
        The pairwise distance matrix, as produced by genetic_distance_matrix().
    tip_labels : List[str]
        List of sequences labels to include in formatted matrix.

    Returns
    -------
    str
        Formatted distance matrix, with column and row names included.
    """
    dist_matrix_df = pd.DataFrame(dist_matrix)
    dist_matrix_df.columns = tip_labels
    dist_matrix_df.index = tip_labels
    dist_matrix_str = dist_matrix_df.to_csv(header = False, sep = '\t')
    dist_matrix_str_w_header = '\t'.join(tip_labels)+'\n'+dist_matrix_str

    return dist_matrix_str_w_header

def nei_saitou(dist_matrix, return_logs = False):

    """
    Computes the unweighted, neighbor-joining tree given a distance matrix.

    Parameters
    ----------
    dist_matrix : np.ndarray
        A square, symmetric matrix of pairwise distances between sequences.
    return_logs : bool, optional
        Whether to return logging information from the algorithm. The default is False.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the edges, in pre-order traversal, in the unrooted tree, with columns ancestor, descendant, and edge_length.
    """

    n_tips = dist_matrix.shape[0] # use this at the end to convert post-order traversal to pre-order
    tip_names = np.array(range(1, dist_matrix.shape[0]+1))
    node_joins = []
    edge_lengths = []
    logs = ''
    iteration = 1

    while dist_matrix.shape[0] > 2:
        logs += 'Iteration ' + str(iteration) + '.\n\n'
        logs += 'Current node IDs are:\n' + str(tip_names) + '\n\n'
        logs += 'D is: \n' + str(dist_matrix) + '\n\n'

        q_matrix = np.zeros(dist_matrix.shape)
        q_coords = list(itertools.combinations(range(dist_matrix.shape[0]), 2))

        # Calculate Q matrix
        for q in q_coords:
            q_matrix[q] = (dist_matrix.shape[0] - 2) * dist_matrix[q] - sum(dist_matrix[q[0]]) - sum(dist_matrix[q[1]])

        # Since the matrix is symmetrical, fill in diag and below with inf
        # With bootstrapping, we might get values of 0 in non-diagonals, so only fill in diagonal and below with inf
        q_matrix[np.tril_indices_from(q_matrix)] = np.inf
        logs += 'Q is:\n' + str(q_matrix) + '\n\n'
        closest_coords = np.array(np.unravel_index(np.argmin(q_matrix), np.array(q_matrix).shape))
        logs += 'Minimum q is ' + str(np.min(q_matrix)) + ' for coords: ' + str(closest_coords) + '\n\n'

        # Distances to node (branch length estimation)
        pair_dist = dist_matrix[tuple(closest_coords)]
        q0_dist = 0.5*pair_dist+ (1/(2*(dist_matrix.shape[0]-2))) * (sum(dist_matrix[closest_coords[0]]) - sum(dist_matrix[closest_coords[1]]))
        q1_dist = pair_dist - q0_dist

        edge_lengths.append(q0_dist)
        edge_lengths.append(q1_dist)

        # Distance matrix update
        # Delete the second row and col of closest coords, overwrite the col remaining after that (consider this the new, "u" col)
        dist_matrix_updt = np.delete(dist_matrix, obj = closest_coords[1], axis = 0)
        dist_matrix_updt = np.delete(dist_matrix_updt, obj = closest_coords[1], axis = 1)

        # Fill in merged row & col w/ new vals
        def u_dist(u, closest_coords, dist_matrix):
            return 0.5 * (dist_matrix[closest_coords[0], u] + dist_matrix[closest_coords[1], u] - dist_matrix[closest_coords[0], closest_coords[1]])
        
        u_dists = [u_dist(u, closest_coords = closest_coords, dist_matrix = dist_matrix) for u in range(dist_matrix.shape[0]) if u != closest_coords[1]] # don't calc distances for closest row/col

        dist_matrix_updt[closest_coords[0]] = u_dists
        dist_matrix_updt[:, closest_coords[0]] = u_dists

        node_joins.append((max(tip_names)+1, tip_names[closest_coords[0]]))
        node_joins.append((max(tip_names)+1, tip_names[closest_coords[1]]))

        # Overwrite 2 merged nodes with single new #
        tip_names[closest_coords[0]] = max(tip_names) + 1
        tip_names = np.delete(tip_names, closest_coords[1])

        # Overwrite distance matrix
        dist_matrix = dist_matrix_updt

        iteration += 1

    # If both remaining nodes are non-tips, order from smallest to largest so they get reordered correctly
    if tip_names[0] > n_tips and tip_names[1] > n_tips:
        node_joins.append(tuple(sorted(tip_names, reverse = True)))
    else:
        node_joins.append(tuple(tip_names))
    edge_lengths.append(dist_matrix[0][1])

    # In post-order traversal initially
    # from_nodes = [x[0] for x in node_joins]
    # to_nodes = [x[1] for x in node_joins]
    # postorder_edges =  pd.DataFrame({'ancestor': from_nodes,'descendant': to_nodes,'edge_length': edge_lengths})

    from_nodes_postorder = [x[0] for x in node_joins]
    postorder_ls = sorted(list(set(from_nodes_postorder)))

    # Convert to pre-order traversal
    # Basically reverse the order of internal node #s in both the ancestor and desc cols 
    node_xwalk = {}
    for i, number in enumerate(postorder_ls[::-1]):
        node_xwalk[number] = postorder_ls[i]

    from_nodes = [node_xwalk[x[0]] for x in node_joins]
    to_nodes = [x[1] if x[1] in range(1, n_tips+1) else node_xwalk[x[1]] for x in node_joins]

    preorder_edges = pd.DataFrame({'ancestor': from_nodes,'descendant': to_nodes,'edge_length': edge_lengths})
    
    if return_logs:
        return preorder_edges, logs # logs are just for debugging
    
    else:
        return preorder_edges
    

def reorient_bootstraps(bootstrapped_sample):
    """
    Reorients bootstrapped DNA samples so that they are arranged by columns.

    Parameters
    ----------
    bootstrapped_sample : List[str]
        A list of DNA sequences.

    Returns
    -------
    List[str]
        A list of reoriented DNA sequences.
    """

    sample = np.array([np.array(list(x), dtype = str) for x in bootstrapped_sample], dtype=str)
    sample_reoriented = [''.join(list(row)) for row in sample.T]

    return sample_reoriented

def bootstrap_sequences(sequences, n_bootstraps = 100):
    """
    Generates bootstrapped DNA samples by randomly selecting columns from the original sequence.

    Parameters
    ----------
    sequences : List[str]
        A list of DNA sequences.
    n_bootstraps : int, optional
        The number of bootstrapped samples to generate. The default is 100.

    Returns
    -------
    List[List[str]]
        A list of bootstrapped DNA samples.
    """

    x = [np.array(list(x), dtype = str) for x in sequences]
    x_np = np.array(x, dtype=str)
    x_T = x_np.T # transpose so we can sample cols as rows (lots of built-in functions for rows)
    msa_cols = np.array([''.join(list(x)) for x in x_T], dtype = str) # concatenate all the separated chars into one long string
    x_bootstrapped = [msa_cols[np.random.randint(len(msa_cols), size=len(msa_cols))] for x in range(n_bootstraps)] # create n bootstrap samples
    bootstrap_samples = list(map(reorient_bootstraps, x_bootstrapped))

    return bootstrap_samples

def build_newick_tree(edge_df, dist_matrix):
    """
    Constructs a Phylogenetic tree given an edge dataframe and a distance matrix.

    Parameters
    ----------
    edge_df : pd.DataFrame
        DataFrame of edges with ancestor and descendant node ids, as well as the length of each edge.
    dist_matrix : np.ndarray
        A square, symmetric numpy array of pairwise distances between sequences used to generate edge_df.

    Returns
    -------
    anytree.node.nodemixin.NodeMixin
        Root node of the Phylogenetic tree object.
    """
    # Kinda ugly implementation bc i didn't read the docs AT ALL *smh*
    tip_names = np.array(range(1, dist_matrix.shape[0]+1))
    inodes = edge_df.ancestor.unique()
    inode_ls = []

    for i,ele in enumerate(inodes):
        inode_children = edge_df.query(f'ancestor == {ele}').descendant.values

        children_as_tn = []
        # If the children are tips, create tree nodes
        children_as_tn += [tn(str(x)) for x in inode_children if x in tip_names]
        
        # If any of the children are internal nodes, grab those internal nodes from the list
        if any(map(lambda v: v in inodes, inode_children)):
            children_as_tn += [x for x in inode_ls if int(x.name) in inodes and int(x.name) in inode_children]

        inode_ls.append(tn(ele.astype(str), children = children_as_tn))

        # Set parent and edge length for each child
        for j in children_as_tn:
            j.parent = inode_ls[i]
            j.length = edge_df.query(f'ancestor == {ele} and descendant == {int(j.name)}')['edge_length'].values[0]

    tr = inode_ls[-1:][0] # Take arbitrary inode as root
    return tr


def get_partitions(tree):
    """
    Given a Phylogenetic tree, returns all possible tip partitions that exist below each node in the tree.

    Parameters:
    -----------
    tree : anytree.node.nodemixin.NodeMixin
        The Phylogenetic tree object containing the partitions to be computed.

    Return:
    -------
    List[List[Tuple[int]]]
        A list of lists containing nodes and all descendent tips in that node group,
        ordered by node height (descending).
    """

    tip_names = [x.name for x in list(tree.tips())]

    # 1. Get all inodes
    tree_inodes = list(tree.non_tips()) + [tree.root()]

    # 2. Get children of each inode
    inode_children = [x.children for x in tree_inodes]

    # 3. Traverse down left and right side of each child
    inode_tips = []
    for child in inode_children:
        inode_tips.append(list(map(lambda y: sorted([int(x) for x in list(child[y].subset())+[child[y].name] if x in tip_names]), range(len(child)))))

    return inode_tips[::-1] # order children from root downwards

def calculate_support(ref_tree_partitions, bootstrap_partitions, ref_newick):
    """
    Calculates bootstrap support scores for partitions of a Phylogenetic tree.

    Parameters
    ----------
    ref_tree_partitions: List[List[int]]
      A list of all possible tip partitions, ordered by node height, in the reference tree.
    bootstrap_partitions: List[List[int]]
      A list of all possible  tip partitions, ordered by node height, in the bootstrap tree.
    ref_newick: anytree.node.nodemixin.NodeMixin
      The reference Phylogenetic tree object containing the partitions against which bootstrap partitions will be compared.

    Returns
    -------
    Dict[str, list]
      A dictionary containing the inode names of nodes in the reference tree and their associated bootstrap support score.
    """

    boot_support_ls = [0]*len(ref_tree_partitions)

    for i, bstrap in enumerate(bootstrap_partitions):
        for j, node in enumerate(bstrap):
            if node in ref_tree_partitions:
                ref_node = ref_tree_partitions.index(bootstrap_partitions[i][j])
                boot_support_ls[ref_node] += 1


    boot_support_np = np.array(boot_support_ls)
    boot_support = np.divide(boot_support_np, n_bootstraps)
    inode_names = [x.name for x in list(ref_newick.non_tips())] + [ref_newick.root().name]
    boot_support_dict = {'inode': inode_names[::-1], # partitions are ordered from root down
                        'bootstrap_support': boot_support}
    
    return boot_support_dict

# -----------
# Main
# -----------

if __name__ == "__main__":
    
    parser = make_argparser()
    args = parser.parse_args()

    msa_w_labels = read_seqs(args.input_fna_filepath, seq_only = False)
    msa_w_labels['seq_id'] = range(1, len(msa_w_labels)+1)
    tip_labels = msa_w_labels[['seq_id', 'lcl']].copy().rename(columns = {'lcl': 'tip_label'})
    tip_labels['tip_label'] = tip_labels.tip_label.str.replace('>', '')

    # Distance matrix
    dist_matrix = genetic_distance_matrix(list(msa_w_labels.sequence.values))
    dist_matrix_str = format_distance_matrix(dist_matrix, tip_labels = list(tip_labels.tip_label.values))
    with open('outputs/genetic-distances.txt', 'w') as textfile:
        textfile.write(dist_matrix_str)

    # Neighbor-joined tree
    edge_df = nei_saitou(dist_matrix)
    edge_df.to_csv('outputs/edges.txt', index = False, header = False, sep = '\t')

    # Bootstrapping
    n_bootstraps = 100
    bootstrap_samples = bootstrap_sequences(msa_w_labels.sequence.values, n_bootstraps = n_bootstraps)
    bootstrap_matrices = [genetic_distance_matrix(x) for x in bootstrap_samples]
    bootstrap_trees = list(map(lambda x: nei_saitou(x), bootstrap_matrices))
    bootstrap_newicks = [build_newick_tree(bootstrap_trees[x], bootstrap_matrices[x]) for x in range(n_bootstraps)]
    ref_tree = build_newick_tree(edge_df, dist_matrix)
    bootstrap_partitions = list(map(lambda x: get_partitions(x), bootstrap_newicks))
    ref_tree_partitions = get_partitions(ref_tree)

    boot_support = calculate_support(ref_tree_partitions=ref_tree_partitions,
                                    bootstrap_partitions= bootstrap_partitions,
                                    ref_newick=ref_tree)

    pd.DataFrame(boot_support).to_csv('outputs/bootstrap.txt', sep = '\t', index = False, header = False)
    
