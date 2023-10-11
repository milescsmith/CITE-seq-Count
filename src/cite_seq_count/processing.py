import gzip
import os
import sys
import time
from collections import Counter, defaultdict
from itertools import islice
from pathlib import Path

import Levenshtein
import pybktree
from numpy import int32
from rich import print as rprint
from scipy import sparse
from tqdm.rich import tqdm
from umi_tools import network, whitelist_methods

from cite_seq_count import secondsToText


def find_best_match(TAG_seq: str, tags: dict[str, str], maximum_distance: int) -> str:
    """
    Find the best match from the list of tags.

    Compares the Levenshtein distance between tags and the trimmed sequences.
    The tag and the sequence must have the same length.
    If no matches found returns 'unmapped'.
    We add 1
    Args:
        TAG_seq (string): Sequence from R1 already start trimmed
        tags (dict): A dictionary with the TAGs as keys and TAG Names as values.
        maximum_distance (int): Maximum distance given by the user.

    Returns:
        best_match (string): The TAG name that will be used for counting.
    """
    best_score = maximum_distance
    for tag, name in tags.items():
        score = Levenshtein.hamming(tag, TAG_seq[: len(tag)])
        if score == 0:
            # Best possible match
            return name
        elif score <= best_score:
            best_score = score
            return name
    return "unmapped"


def find_best_match_shift(TAG_seq: str, tags: dict[str, str], maximum_distance: int) -> str:
    """
    Find the best match from the list of tags with sliding window.

    Compares the Levenshtein distance between tags and the trimmed sequences.
    The tag and the sequence must have the same length.
    If no matches found returns 'unmapped'.
    We add 1
    Args:
        TAG_seq (string): Sequence from R1 already start trimmed
        tags (dict): A dictionary with the TAGs as keys and TAG Names as values.
        maximum_distance (int): Maximum distance given by the user.

    Returns:
        best_match (string): The TAG name that will be used for counting.
    """
    best_score = maximum_distance
    shifts = range(len(TAG_seq) - len(max(tags, key=len)))

    for shift in shifts:
        for tag, name in tags.items():
            score = Levenshtein.hamming(tag, TAG_seq[shift : len(tag) + shift])
            if score == 0:
                # Best possible match
                return name
            elif score <= best_score:
                best_score = score
                return name
    return "unmapped"


def map_reads(
    read1_path: Path,
    read2_path: Path,
    tags: dict[str, str],
    barcode_slice: slice,
    umi_slice: slice,
    indexes: tuple[int, int],
    debug: bool,
    start_trim: int,
    maximum_distance: int,
    sliding_window: bool,
) -> tuple[dict[str, str], Counter[str]]:
    """Read through R1/R2 files and generate a islice starting at a specific index.

    It reads both Read1 and Read2 files, creating a dict based on cell barcode.

    Args:
        read1_path (string): Path to R1.fastq.gz
        read2_path (string): Path to R2.fastq.gz
        chunk_size (int): The number of lines to process
        tags (dict): A dictionary with the TAGs + TAG Names.
        barcode_slice (slice): A slice for extracting the Barcode portion from the
            sequence.
        umi_slice (slice): A slice for extracting the UMI portion from the
            sequence.
        indexes (list): Pair of first and last index for islice
        whitelist (set): The set of white-listed barcodes.
        debug (bool): Print debug messages. Default is False.
        start_trim (int): Number of bases to trim at the start.
        maximum_distance (int): Maximum distance given by the user.
        sliding_window (bool): A bool enabling a sliding window search

    Returns:
        results (dict): A dict of dict of Counters with the mapping results.
        no_match (Counter): A counter with unmapped sequences.
    """
    # Initiate values
    results = {}
    no_match = Counter()
    n = 1
    t = time.time()
    with gzip.open(read1_path, "rt") as textfile1, gzip.open(read2_path, "rt") as textfile2:
        # Read all 2nd lines from 4 line chunks. If first_n not None read only 4 times the given amount.
        secondlines = islice(zip(textfile1, textfile2, strict=True), indexes[0] * 4 + 1, indexes[1] * 4 + 1, 4)
        for read1, read2 in secondlines:
            read1 = read1.strip()
            read2 = read2.strip()

            # Progress info
            if n % 1000000 == 0:
                rprint(
                    f"Processed 1,000,000 reads in {secondsToText.secondsToText(time.time() - t)}. Total "
                    f"reads: {n:,} in child {os.getpid()}"
                )
                sys.stdout.flush()
                t = time.time()

            # Get cell and umi barcodes.
            cell_barcode = read1[barcode_slice]
            # This change in bytes is required by umi_tools for umi correction
            umi = bytes(read1[umi_slice], "ascii")
            # Trim potential starting sequences
            tag_seq = read2[start_trim:]

            if cell_barcode not in results:
                results[cell_barcode] = defaultdict(Counter)

            if sliding_window:
                best_match = find_best_match_shift(tag_seq, tags, maximum_distance)
            else:
                best_match = find_best_match(tag_seq, tags, maximum_distance)

            results[cell_barcode][best_match][umi] += 1

            if best_match == "unmapped":
                no_match[tag_seq] += 1

            if debug:
                rprint(
                    f"\nline:{read1 + read2}\ncell_barcode:{cell_barcode}\tUMI:{umi}\tTAG_seq:{tag_seq}\nline length:{len(read1 + read2)}\tcell barcode length:{len(cell_barcode)}\tUMI length:{len(umi)}\tTAG sequence length:{len(tag_seq)}\nBest match is: {best_match}"
                )
                sys.stdout.flush()
            n += 1
    rprint(f"Mapping done for process {os.getpid()}. Processed {n - 1:,} reads")
    sys.stdout.flush()
    return (results, no_match)


def merge_results(parallel_results):
    """Merge chunked results from parallel processing.

    Args:
        parallel_results (list): List of dict with mapping results.

    Returns:
        merged_results (dict): Results combined as a dict of dict of Counters
        umis_per_cell (Counter): Total umis per cell as a Counter
        reads_per_cell (Counter): Total reads per cell as a Counter
        merged_no_match (Counter): Unmapped tags as a Counter
    """
    merged_results: defaultdict[str, Counter[str]] = {}
    merged_no_match: Counter[str] = Counter()
    umis_per_cell: Counter[str] = Counter()
    reads_per_cell: Counter[str] = Counter()
    for chunk in parallel_results:
        mapped = chunk[0]
        unmapped = chunk[1]
        for cell_barcode in mapped:
            if cell_barcode not in merged_results:
                merged_results[cell_barcode] = defaultdict(Counter)
            for tag in mapped[cell_barcode]:
                # Test the counter. Returns false if empty
                if mapped[cell_barcode][tag]:
                    for umi in mapped[cell_barcode][tag]:
                        merged_results[cell_barcode][tag][umi] += mapped[cell_barcode][tag][umi]
                        umis_per_cell[cell_barcode] += len(mapped[cell_barcode][tag])
                        reads_per_cell[cell_barcode] += mapped[cell_barcode][tag][umi]
        merged_no_match |= unmapped
    return (merged_results, umis_per_cell, reads_per_cell, merged_no_match)


def correct_umis(final_results, collapsing_threshold, top_cells, max_umis):
    """
    Corrects umi barcodes within same cell/tag groups.

    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        collapsing_threshold (int): Max distance between umis.
        top_cells (set): Set of cells to go through.
        max_umis (int): Maximum UMIs to consider for one cluster.

    Returns:
        final_results (dict): Same as input but with corrected umis.
        corrected_umis (int): How many umis have been corrected.
        aberrant_umi_count_cells (set): Set of uncorrected cells.
    """
    rprint("Correcting umis")
    corrected_umis = 0
    aberrant_umi_count_cells = set()
    for cell_barcode in top_cells:
        for tag in final_results[cell_barcode]:
            n_umis = len(final_results[cell_barcode][tag])
            if n_umis > 1 and n_umis <= max_umis:
                umi_clusters = network.UMIClusterer()
                umiclusters = umi_clusters(final_results[cell_barcode][tag], collapsing_threshold)
                (new_res, temp_corrected_umis) = update_umi_counts(umiclusters, final_results[cell_barcode][tag])
                final_results[cell_barcode][tag] = new_res
                corrected_umis += temp_corrected_umis
            elif n_umis > max_umis:
                aberrant_umi_count_cells.add(cell_barcode)
    return (final_results, corrected_umis, aberrant_umi_count_cells)


def update_umi_counts(UMIclusters, cell_tag_counts):
    """
    Update a dict object with umis corrected.

    Args:
        UMIclusters (list): List of lists with corrected umis
        cell_tag_counts (Counter): Counter of umis

    Returns:
        cell_tag_counts (Counter): Updated Counter of umis
        temp_corrected_umis (int): Number of corrected umis
    """
    temp_corrected_umis = 0
    for umi_cluster in UMIclusters:  # This is a list with the first element the dominant barcode
        if len(umi_cluster) > 1:  # This means we got a correction
            major_umi = umi_cluster[0]
            for minor_umi in umi_cluster[1:]:
                temp_corrected_umis += 1
                temp = cell_tag_counts.pop(minor_umi)
                cell_tag_counts[major_umi] += temp
    return (cell_tag_counts, temp_corrected_umis)


def collapse_cells(true_to_false, umis_per_cell, final_results, ab_map):
    """
    Collapses cell barcodes based on the mapping true_to_false

    Args:
        true_to_false (dict): Mapping between the reference and the "mutated" barcodes.
        umis_per_cell (Counter): Counter of number of umis per cell.
        final_results (dict): Dict of dict of Counters with mapping results.
        ab_map (dict): Dict of the TAGS.

    Returns:
        umis_per_cell (Counter): Counter of number of umis per cell.
        final_results (dict): Same as input but with corrected cell barcodes.
        corrected_barcodes (int): How many cell barcodes have been corrected.
    """
    rprint("Collapsing cell barcodes")
    corrected_barcodes = 0
    for real_barcode in true_to_false:
        # If the cell barcode is not in the results
        if real_barcode not in final_results:
            final_results[real_barcode] = defaultdict()
            for tag in ab_map:
                final_results[real_barcode][tag] = Counter()
        for fake_barcode in true_to_false[real_barcode]:
            temp = final_results.pop(fake_barcode)
            corrected_barcodes += 1
            for tag in temp.keys():
                final_results[real_barcode][tag].update(temp[tag])
            temp_umi_counts = umis_per_cell.pop(fake_barcode)
            # temp_read_counts = reads_per_cell.pop(fake_barcode)

            umis_per_cell[real_barcode] += temp_umi_counts
            # reads_per_cell[real_barcode] += temp_read_counts

    return (umis_per_cell, final_results, corrected_barcodes)


def correct_cells(
    final_results,
    reads_per_cell,
    umis_per_cell,
    collapsing_threshold,
    expected_cells,
    ab_map,
):
    """
    Corrects cell barcodes.

    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        umis_per_cell (Counter): Counter of number of umis per cell.
        collapsing_threshold (int): Max distance between umis.
        expected_cells (int): Number of expected cells.
        ab_map (dict): Dict of the TAGS.

    Returns:
        final_results (dict): Same as input but with corrected umis.
        umis_per_cell (Counter): Counter of umis per cell after cell barcode correction
        corrected_umis (int): How many umis have been corrected.
    """
    rprint("Looking for a whitelist")
    cell_whitelist, true_to_false = whitelist_methods.getCellWhitelist(
        cell_barcode_counts=reads_per_cell,
        expect_cells=expected_cells,
        cell_number=expected_cells,
        error_correct_threshold=collapsing_threshold,
        plotfile_prefix=False,
    )

    (umis_per_cell, final_results, corrected_barcodes) = collapse_cells(
        true_to_false=true_to_false,
        umis_per_cell=umis_per_cell,
        final_results=final_results,
        ab_map=ab_map,
    )
    return (final_results, umis_per_cell, corrected_barcodes)


def correct_cells_whitelist(final_results, umis_per_cell, whitelist, collapsing_threshold, ab_map):
    """
    Corrects cell barcodes.

    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        umis_per_cell (Counter): Counter of UMIs per cell.
        whitelist (set): The whitelist reference given by the user.
        collapsing_threshold (int): Max distance between umis.
        ab_map (OrederedDict): Tags in an ordered dict.


    Returns:
        final_results (dict): Same as input but with corrected umis.
        umis_per_cell (Counter): Updated UMI counts after correction.
        corrected_barcodes (int): How many umis have been corrected.
    """
    barcode_tree = pybktree.BKTree(Levenshtein.hamming, whitelist)
    rprint("Generated barcode tree from whitelist")
    cell_barcodes = list(final_results.keys())
    n_barcodes = len(cell_barcodes)
    rprint("Finding reference candidates")
    rprint(f"Processing {n_barcodes:,} cell barcodes")

    # Run with one process
    true_to_false = find_true_to_false_map(
        barcode_tree=barcode_tree,
        cell_barcodes=cell_barcodes,
        whitelist=whitelist,
        collapsing_threshold=collapsing_threshold,
    )
    (umis_per_cell, final_results, corrected_barcodes) = collapse_cells(
        true_to_false, umis_per_cell, final_results, ab_map
    )
    return (final_results, umis_per_cell, corrected_barcodes)


def find_true_to_false_map(barcode_tree, cell_barcodes, whitelist, collapsing_threshold):
    """
    Creates a mapping between "fake" cell barcodes and their original true barcode.

    Args:
        barcode_tree (BKTree): BKTree of all original cell barcodes.
        cell_barcodes (List): Cell barcodes to go through.
        whitelist (Set): Set of the whitelist, the "true" cell barcodes.
        collasping_threshold (int): How many mistakes to correct.

    Return:
        true_to_false (defaultdict(list)): Contains the mapping between the fake and real barcodes. The key is the real one.
    """
    true_to_false = defaultdict(list)
    for cell_barcode in tqdm(cell_barcodes):
        if cell_barcode in whitelist:
            # if the barcode is already whitelisted, no need to add
            continue
        # get all members of whitelist that are at distance of collapsing_threshold
        candidates = [white_cell for d, white_cell in barcode_tree.find(cell_barcode, collapsing_threshold) if d > 0]
        if len(candidates) == 1:
            white_cell_str = candidates[0]
            true_to_false[white_cell_str].append(cell_barcode)
        else:
            # the cell doesnt match to any whitelisted barcode,
            # hence we have to drop it
            # (as it cannot be asscociated with any frequent barcode)
            continue
    return true_to_false


def generate_sparse_matrices(
    final_results: defaultdict[str, defaultdict[str, Counter[str]]],
    ordered_tags_map: dict[str, int],
    top_cells: set[str],
) -> tuple[sparse.dok_matrix, sparse.dok_matrix]:
    """
    Create two sparse matrices with umi and read counts.

    Args:
        final_results (dict): Results in a dict of dicts of Counters.
        ordered_tags_map (dict): Tags in order with indexes as values.

    Returns:
        umi_results_matrix (scipy.sparse.dok_matrix): UMI counts
        read_results_matrix (scipy.sparse.dok_matrix): Read counts

    """
    umi_results_matrix = sparse.dok_matrix((len(ordered_tags_map), len(top_cells)), dtype=int32)
    read_results_matrix = sparse.dok_matrix((len(ordered_tags_map), len(top_cells)), dtype=int32)
    for i, cell_barcode in enumerate(top_cells):
        for tag in final_results[cell_barcode]:
            if final_results[cell_barcode][tag]:
                umi_results_matrix[ordered_tags_map[tag], i] = len(final_results[cell_barcode][tag])
                read_results_matrix[ordered_tags_map[tag], i] = sum(final_results[cell_barcode][tag].values())
    return umi_results_matrix, read_results_matrix
