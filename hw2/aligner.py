# ------------
# Imports
# ------------

# !pip install fastaparser
import fastaparser
import pandas as pd
import numpy as np
import re
import argparse

# ------------------
# Custom functions
# ------------------

def make_argparser():
    '''
    Parameters
    -----------
        Technically, this function takes no parameters. Note, though, that the parser expects 3 strings and 3 integers. The final flag, --ignore_outer_gaps, is optional

    Return
    --------
        This function returns an argparse parser, which is then used as the input for the needleman_wunsch() function. See the docs for more info on how the argparse methods work: https://docs.python.org/3/library/argparse.html
    '''
    parser = argparse.ArgumentParser(prog = 'aligner.py', 
                                     description = 'A program to align a query sequence to a reference sequence, both contained in fasta (.fna) files')
    
    parser.add_argument('-q', '--query', required=True, type=str, help = 'Path to query fasta file, including .fna extension')
    parser.add_argument('-r', '--reference', required=True, type=str, help = 'Path to reference fasta file, including .fna extension')
    parser.add_argument('-o', '--output', required=True, type=str, help = 'Path to file output by the program, including the .txt. extension')
    parser.add_argument('-g', '--gap_penalty', required=True, type=int, help = 'Fixed gap penalty for all positions, unless the --ignore_outer_gaps flag is set; then, this argument only applies to non-start and non-terminal gaps')
    parser.add_argument('-p', '--mismatch_penalty', required=True, type=int, help = 'Fixed penalty for any mismatched bases or residues between the query and the reference')
    parser.add_argument('-m', '--match_score', required=True, type=int, help = 'Fixed score for any matching bases or residues between the query and the reference')
    parser.add_argument('--ignore_outer_gaps', action='store_true')
    
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

def needleman_wunsch(
    ref_sequence: str,
    query_sequence: str,
    gap_penalty: int,
    mismatch_penalty: int,
    match_score: int,
    free_initial_gaps: bool = False
) -> tuple:
    """
    Calculates the optimal alignment of two sequences using the Needleman-Wunsch algorithm.

    Parameters
    ----------
    ref_sequence : str
        The reference sequence of the alignment.
    query_sequence : str
        The query sequence to align to the reference.
    gap_penalty : int
        The penalty for introducing a gap in either sequence.
    mismatch_penalty : int
        The penalty for mismatching two bases or residues.
    match_score : int
        The score for matching two bases or residues.
    free_initial_gaps : bool
        Whether or not to allow gaps at the beginning or end of the sequences without penalty.

    Returns
    -------
    tuple
        A tuple containing two matrices: 
        1. A matrix with the optimal alignment score at every position for each prefix of the two sequences
        2. A matrix with the coordinates of the previous position that gave the optimal alignment score at every position for each prefix of the two sequences.
        Both matrices are of shape (len(query_sequence) + 1, len(ref_sequence) + 1) and numpy.ndarrays of dtype int.

    Raises
    ------
    None

    Notes
    -----
    The sequences are assumed to be zero-indexed.

    The gap penalty is usually a negative number, while the mismatch penalty and match score are usually positive.

    This function implements the Needleman-Wunsch algorithm [1], also known as the global alignment algorithm, which finds the optimal alignment of two sequences based on a similarity matrix or scoring scheme that assigns scores or penalties to pairs of symbols in the sequences. It does so by searching for the best alignment among all possible alignments of two sequences.

    The algorithm works as follows:

    1. Initialize two matrices: one with the optimal alignment score at each position for each prefix of the two sequences, and another with the coordinates of the previous position that gave the best score at each position for each prefix of the two sequences. Both matrices have dimensions (len(query_sequence) + 1, len(ref_sequence) + 1), where the first row/column corresponds to gaps or prefix/suffix of the sequences.

    2. Fill in the first row/column of the matrices with the gap penalty times the position. For the first row, this corresponds to introducing gaps in the reference sequence. For the first column, this corresponds to introducing gaps in the query sequence.

    3. For each combination of two symbols in the sequences, calculate the score of matching them (if they match) or mismatching them (if they differ) according to the scoring scheme, and update the optimal score and coordinates matrices based on the best of three possibilities: 
        a. Align the symbols by matching them and add the matching score to the previous optimal score at the diagonal position of the matrix
        b. Align the query symbol to a gap in the reference sequence and add the gap penalty to the previous optimal score at the position above in the matrix
        c. Align the reference symbol to a gap in the query sequence and add the gap penalty to the previous optimal score at the position to the left in the matrix.

    4. Trace back the optimal path in the coordinates matrix from the bottom right cell to the top left cell, following the position that gave the best score at each step. This yields the optimal alignment.

    References
    ----------
    [1] Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology, 48(3), 443-453.
    """

    score_matrix = np.zeros((len(query_sequence) + 1, len(ref_sequence) + 1))  # add 1 to account for initial cell in upper left
    coord_matrix = np.zeros((len(query_sequence) + 1, len(ref_sequence) + 1, 2), dtype = int)

    # Create row 1 and col 1 of each matrix, which requires knowing only the gap penalty
    initial_gap_penalty = gap_penalty if not free_initial_gaps else 0
    row1 = np.array([x * initial_gap_penalty for x in range(len(ref_sequence) + 1)])
    col1 = np.array([x * initial_gap_penalty for x in range(len(query_sequence) + 1)])
    coord_row1 = np.array([np.array([0, 0])] + [np.array([0, x]) for x in range(len(ref_sequence))])
    coord_col1 = np.array([np.array([0, 0])] + [np.array([x, 0]) for x in range(len(query_sequence))])  # originally seq1

    # Initialize matrices
    score_matrix[0, :] = row1
    score_matrix[:, 0] = col1.T
    coord_matrix[0, :, :] = coord_row1
    coord_matrix[:, 0, :] = coord_col1

    for row in range(1, len(query_sequence) + 1):
        for col in range(1, len(ref_sequence) + 1):

            # Where did the score originate from? Store the coords of each position (diagonal, above, left)
            candidate_coords = [(row - 1, col - 1), (row - 1, col), (row, col - 1)]

            # Calculate scores for each direction (diagonal, left, above)
            dscore = (score_matrix[candidate_coords[0]] + match_score if ref_sequence[col - 1] == query_sequence[row - 1] else score_matrix[candidate_coords[0]] + mismatch_penalty) # -1 to get us back to indexing from 0
            ascore = score_matrix[candidate_coords[1]] + gap_penalty
            lscore = score_matrix[candidate_coords[2]] + gap_penalty

            candidate_scores = np.array([dscore, ascore, lscore])

            # Get highest score and which direction it came from (diagonal, above, left)
            winning_score = np.max(candidate_scores)
            winning_position = np.argmax(candidate_scores)

            # Fill in matrices with best score and coordinate from which it originated
            score_matrix[row, col] = winning_score
            coord_matrix[row, col, :] = np.array(candidate_coords[winning_position], dtype = int)

    return score_matrix, coord_matrix

def retrieve_optimal_path(ref_sequence, query_sequence, score_matrix, free_terminal_gaps):
    """
    Given the score matrix generated via the Needleman-Wunsch algorithm, return the optimal path of the alignment and the alignment score. 

    Parameters
    ----------
    ref_sequence : str
        The reference sequence of the alignment.
    query_sequence : str
        The query sequence to align to the reference.
    score_matrix : numpy.ndarray
        A matrix that stores the optimal alignment score at every position for each prefix of the two sequences. The matrix is obtained through Needleman-Wunsch algorithm.
    free_terminal_gaps : bool
        Whether to allow gaps at the beginning or end of the sequences without penalty. 

    Returns
    -------
    tuple
        A tuple where the first element is a list containing the coordinates of previous position that gave the optimal alignment score at every position for each prefix of the two sequences, starting from the bottom-right position and going leftward/upward until the upper-left position is reached. The second element of the tuple is the alignment score.

    Notes
    -----
    This function uses the optimal score matrix to trace back the optimal alignment. Starting from the bottom right corner of the matrix, it traverses a series of positions on the matrix until it reaches the upper left corner, collecting the position at each step that gave the optimal score. It also keeps track of the alignment score while doing so. Once it has collected all the positions, it returns them in a list that can be used to construct the optimal alignment.

    Following are the steps performed by the function:
    - If free_terminal_gaps are allowed, find the cell with highest score in the last row and col, which will correspond to the final diagonal path taken to fill them. 
    - Initialize 'optimal_coords' w/ coords of best terminal score + the coords in that position of the matrix (the next step backwards in the path); traverse the matrices and store the position every time a diagonal path is found by comparing current and previous positions
    - Calculate and return the alignment score.
    """

    if free_terminal_gaps:
        # Create a matrix with every score except the final row and col overwritten w/ very neg. #
        # Basically so that the argmax only looks in the final row or col for the best score
        # Then find the best score in the entire matrix
        terminal_matrix = score_matrix.copy()
        terminal_matrix[:-1, :-1] = float('-inf')
        best_terminal_score = np.argmax(terminal_matrix)
        best_terminal_coords = np.array(np.unravel_index(best_terminal_score, np.array(terminal_matrix).shape))
        alignment_score = terminal_matrix[best_terminal_coords[0], best_terminal_coords[1]]

    else:
        # Start from bottom right corner if terminal gaps are penalized
        alignment_score = score_matrix[int(len(query_sequence)), int(len(ref_sequence))]
        best_terminal_coords = np.array([int(len(query_sequence)), int(len(ref_sequence))])

    # Initialize w/ coords of best terminal score + the coords in that position of the matrix (ie the next step backwards in the path)
    optimal_coords = [best_terminal_coords] + [coord_matrix[best_terminal_coords[0], best_terminal_coords[1], :]]
    next_coords = coord_matrix[best_terminal_coords[0], best_terminal_coords[1], :]

    # Traverse the matrices until the upper left corner is reached; then stop
    while not np.array_equal(next_coords, np.array([0, 0])):

        coords = coord_matrix[next_coords[0], next_coords[1], :]

        optimal_coords.append(coords)
        next_coords = coords

    optimal_coords_fwd = optimal_coords[::-1]  # Reverse the coords so they start from upper left

    return optimal_coords_fwd, alignment_score

def retrieve_alignments(optimal_coords: list, 
                        ref_sequence: str, 
                        query_sequence: str,
                        bottom_right_coords = np.array) -> tuple:
    """
    Given a list with the path traversed in the matrix to get the optimal alignment, constructs the 2 strings representing the optimal alignment, inserting gaps at specific positions in either the given reference or query sequence.

    Parameters
    ----------
    optimal_coords : list
        A list containing the coordinates of previous position that gave the optimal alignment score at every position for each prefix of the two sequences, starting from the bottom right position and going leftward/upward until the upper left position is reached.
    ref_sequence : str
        The reference sequence of the alignment.
    query_sequence : str
        The query sequence to align to the reference.
    bottom_right_coords : numpy.ndarray
        The coordinates of the bottom-right cell of the score matrix obtained through Needleman-Wunsch algorithm.

    Returns
    -------
    tuple
        A tuple consisting of the two strings representing the sentences in the optimal path, having gaps (if applicable) inserted at the positions given by optimal_coords. 

    Notes
    -----
    This function takes the list of optimal coordinates obtained from retrieve_optimal_path() function, constructs optimal alignment of the two sequences by inserting gaps at the specific positions in either sequence as dictated by optimal_coords.

    Following the optimal_coords path, compared at each position to the preceding 
    position, where the aim here is to get the optimal alignment of each sentence. Two strings, 
    `query_aligned` and `ref_aligned` are returned, with gaps where appropriate. 

    If `bottom_right_coords` is provided, algebraic steps are done after traversing the optimal path to ensure any missing terminal gaps are added to the output.

    """

    ref_aligned = ""
    query_aligned = ""

    ref_init = " " + ref_sequence
    query_init = " " + query_sequence

    ref_consumed = 0
    query_consumed = 0

    for i in range(1, len(optimal_coords)):

        # Compare current coords to previous coords to get dir the path came from
        current_row = optimal_coords[i][0]
        current_col = optimal_coords[i][1]
        row_dir = optimal_coords[i - 1][0] - current_row
        col_dir = optimal_coords[i - 1][1] - current_col

        # Path went diagonally
        if row_dir == -1 and col_dir == -1:
            ref_aligned += ref_init[current_col]
            query_aligned += query_init[current_row]

            ref_consumed += 1
            query_consumed += 1

        # Path went right (came from left)
        elif col_dir == -1 and row_dir == 0:
            ref_aligned += ref_init[current_col]
            query_aligned += "_"  # insert a gap

            ref_consumed += 1

        # Path went down (came from above)
        elif col_dir == 0 and row_dir == -1:
            ref_aligned += "_"
            query_aligned += query_init[current_row]

            query_consumed += 1

    # If terminal gaps were free, need to add those on
    if not np.array_equal(optimal_path[-1:], bottom_right_coords):
        # If x of final coords = last row, then add gaps to q and consume ref
        if optimal_path[-1:][0][0] == bottom_right_coords[0]:
            outstanding_len = len(ref_sequence) - ref_consumed
            query_aligned += "_"*outstanding_len
            ref_aligned += ref_sequence[ref_consumed:]

        # If y of final coords = last col, then add gaps to ref and consume q
        if optimal_path[-1:][0][1] == bottom_right_coords[1]:
            outstanding_len = len(query_sequence) - query_consumed
            ref_aligned += "_"*outstanding_len
            query_aligned += query_sequence[query_consumed:]

    return ref_aligned, query_aligned


# ------------------
# Main
# ------------------

if __name__ == "__main__":
    parser = make_argparser()
    args = parser.parse_args()

    # If the input is an fna, use fastaparser + custom functions
    # Otherwise, if it's a .aa, assume the first line is the header and the 2nd is the sequence
    if args.reference[-3:] == 'fna':
        seq1_dict = read_seqs(args.reference, seq_only = False)
        seq1 = seq1_dict['sequence'][0]
        seq1_header = seq1_dict['lcl'][0]
    
    elif args.reference[-2:] == 'aa':
        seq1_df = pd.read_csv(args.reference)
        seq1 = seq1_df.iloc[0][0]
        seq1_header = seq1_df.columns[0]

    if args.query[-3:] == 'fna':
        seq2_dict = read_seqs(args.query, seq_only = False)
        seq2 = seq2_dict['sequence'][0]
        seq2_header = seq2_dict['lcl'][0]
    
    elif args.query[-2:] == 'aa':
        seq2_df = pd.read_csv(args.query)
        seq2 = seq2_df.iloc[0][0]
        seq2_header = seq2_df.columns[0]

    # Run the Needleman-Wunsch algo
    score_matrix, coord_matrix = needleman_wunsch(
                                                ref_sequence = seq1, # cols
                                                query_sequence = seq2, # rows
                                                gap_penalty = args.gap_penalty,
                                                mismatch_penalty = args.mismatch_penalty,
                                                match_score = args.match_score,
                                                free_initial_gaps = args.ignore_outer_gaps
                                                )

    # Get the bottom right corner
    last_row = score_matrix.shape[0]-1
    last_col = score_matrix.shape[1]-1

    # Retrieve optimal path and corresponding alignments; generate viz
    optimal_path, alignment_score = retrieve_optimal_path(ref_sequence = seq1, 
                                                          query_sequence = seq2, 
                                                          score_matrix = score_matrix, 
                                                          free_terminal_gaps = args.ignore_outer_gaps)
    
    # Retrieve alignments using the coord_matrix that correspond to the optimal path
    ref_aligned, query_aligned = retrieve_alignments(optimal_coords = optimal_path, 
                                                     ref_sequence = seq1, 
                                                     query_sequence = seq2, 
                                                     bottom_right_coords = np.array([last_row, last_col]))

    align_viz = list(map(lambda x, y: '|' if x == y else ' ' if x == '_' or y == '_' else 'x', ref_aligned, query_aligned))
    # alignment_output = str(int(alignment_score)) + '\n' + seq1_header + '\n' + ref_aligned + '\n' + ''.join(align_viz) + '\n' + query_aligned + '\n' + seq2_header
    alignment_output = str(int(alignment_score)) + '\n' + seq2_header + '\n' + query_aligned + ''.join(align_viz) + '\n' + ref_aligned + '\n' + seq1_header

    # Write out result
    with open(args.output, 'w') as text_file:
        text_file.write(alignment_output)
