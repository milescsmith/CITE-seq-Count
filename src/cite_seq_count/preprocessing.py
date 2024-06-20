import csv
import gzip
import sys
from collections import OrderedDict
from collections.abc import Callable, Generator
from functools import partial
from itertools import combinations, islice
from math import floor
from pathlib import Path
from typing import BinaryIO, IO

import Levenshtein
import regex
from loguru import logger
from rich import print as rprint

STRIP_CHARS = '"0123456789- \t\n'

gzopen = partial(gzip.open, mode="rt")


class MismatchedReadlengthsError(Exception):
    def __init__(self, message: str="Read lengths do not match") -> None:
        self.message = message
        super().__init__(message)


class NotMultipleofFourError(Exception):
    def __init__(self, message: str="Number of FASTQ lines is not multiple of 4") -> None:
        self.message = message
        super().__init__(message)


def get_indexes(start_index: int, chunk_size: int, nth: int) -> list[int]:
    """
    Creates indexes from a reference index, a chunk size an nth number

    Args:
        start_index (int): first position
        chunk_size (int): Chunk size
        nth (int): The nth number

    Returns:
        list: First and last position of indexes
    """
    start_index = nth * chunk_size
    stop_index = chunk_size + nth * chunk_size
    return [start_index, stop_index]


def chunk_reads(n_reads: int, n: int) -> list[list[int]]:
    """
    Creates a list of indexes for the islice iterator from the map_reads function.

    Args:
        n_reads (int): Number of reads to split
        n (int): How many buckets for the split.
    Returns:
        indexes (list(list)): Each entry contains the first and the last index for a read.
    """
    if n_reads % n == 0:
        chunk_size = n_reads // n
        rest = 0
    else:
        chunk_size = floor(n_reads / n)
        rest = n_reads - (n * chunk_size)
    indexes = [get_indexes(i, chunk_size, i) for i in range(n)]
    indexes[-1][1] += rest
    return indexes


def parse_whitelist_csv(filename: Path, barcode_length: int, collapsing_threshold: int) -> tuple[set[str], int]:
    """Reads white-listed barcodes from a CSV file.

    The function accepts plain barcodes or even 10X style barcodes with the
    `-1` at the end of each barcode.

    Args:
        filename (str): Whitelist barcode file.
        barcode_length (int): Length of the expected barcodes.
        collapsing_threshold (int): Maximum distance to collapse cell barcodes.

    Returns:
        set: The set of white-listed barcodes.
        int: Collasping threshold

    """
    cell_pattern = regex.compile(fr"[ATGC]{{{barcode_length}}}")

    file_open: Callable[[str | Path, str], IO] = gzopen if filename.suffix == ".gz" else open
    with file_open(filename) as csv_file:
        whitelist = [
            row[0].strip(STRIP_CHARS)
            for row in csv.reader(csv_file)
            if (len(row[0].strip(STRIP_CHARS)) == barcode_length)
        ]

    for cell_barcode in whitelist:
        if not cell_pattern.match(cell_barcode):
            sys.exit(f"This barcode {cell_barcode} is not only composed of ATGC bases.")
    # collapsing_threshold=test_cell_distances(whitelist, collapsing_threshold)
    if not whitelist:
        sys.exit("Please check cell barcode indexes -cbs, -cbl because none of the given whitelist is valid.")
    return (set(whitelist), collapsing_threshold)


def test_cell_distances(whitelist: set[str], collapsing_threshold: int) -> int:
    """Tests cell barcode distances to validate provided cell barcode collapsing threshold

    Function needs the given whitelist as well as the threshold.
    If the value is too high, it will rerun until an acceptable value is found.

    Args:
        whitelist (set): Whitelist barcode set
        collapsing_threshold (int): Value of threshold

    Returns:
        collapsing_threshold (int): Valid threshold
    """
    ok = False
    while not ok:
        rprint(f"Testing cell barcode collapsing threshold of {collapsing_threshold}")
        all_comb = combinations(whitelist, 2)
        for comb in all_comb:
            if Levenshtein.hamming(comb[0], comb[1]) <= collapsing_threshold:
                collapsing_threshold -= 1
                rprint("Value is too high, reducing it by 1")
                break
        else:
            ok = True
    rprint(f"Using {collapsing_threshold} for cell barcode collapsing threshold")
    return collapsing_threshold


def parse_tags_csv(filename: Path) -> dict[str:str]:
    """Reads the TAGs from a CSV file.

    The expected file format (no header) is: TAG,TAG_NAME.
    e.g. file content
        GTCAACTCTTTAGCG,Hashtag_1
        TGATGGCCTATTGGG,Hashtag_2
        TTCCGCCTCTCTTTG,Hashtag_3

    Args:
        filename (str): TAGs file.

    Returns:
        dict: A dictionary containing the TAGs and their names.

    """
    file_open = gzopen if filename.suffix == ".gz" else open
    with file_open(filename) as csv_file:
        tags = {row[0].strip(): row[1].strip() for row in csv.reader(csv_file)}
    return tags


def check_tags(tags: dict[str, str], maximum_distance: int) -> OrderedDict[str, str]:
    """Evaluates the distance between the TAGs based on the `maximum distance`
    argument provided.

    Additionally, it adds the barcode to the name of the TAG circumventing the
    need of having to share the mapping of the antibody and the barcode.

    The output will have the keys sorted by TAG length (longer first). This
    way, longer barcodes will be evaluated first.

    Args:
        tags (dict): A dictionary with the TAGs + TAG Names.
        maximum_distance (int): The maximum Levenshtein distance allowed
            between two TAGs.

    Returns:
        collections.OrderedDict: An ordered dictionary containing the TAGs and
            their names in descendent order based on the length of the TAGs.

    """
    ordered_tags = OrderedDict()
    for tag in sorted(tags, key=len, reverse=True):
        ordered_tags[tag] = f"{tags[tag]}-{tag}"
    # If only one TAG is provided, then no distances to compare.
    if len(tags) == 1:
        return ordered_tags

    offending_pairs = []
    for a, b in combinations(ordered_tags.keys(), 2):
        distance = Levenshtein.distance(a, b)
        if distance <= (maximum_distance - 1):
            offending_pairs.append([a, b, distance])
    dna_pattern = regex.compile("^[ATGC]*$")
    for tag in tags:
        if not dna_pattern.match(tag):
            rprint(f"This tag {tag} is not only composed of ATGC bases.\nPlease check your tags file")
            sys.exit("Exiting the application.\n")
    # If offending pairs are found, print them all.
    if offending_pairs:
        rprint(
            "[ERROR] Minimum Levenshtein distance of TAGs barcode is less "
            "than given threshold.\n"
            "Please use a smaller distance.\n\n"
            "Offending case(s):\n"
        )
        for pair in offending_pairs:
            rprint(f"\t{pair[0]}\n\t{pair[1]}\n\tDistance = {pair[2]}\n")
        sys.exit("Exiting the application.\n")
    return ordered_tags


def get_read_length(filename: Path) -> int:
    """Check wether SEQUENCE lengths are consistent in a FASTQ file and return
    the length.

    Args:
        filename (str): FASTQ file.

    Returns:
        int: The file's SEQUENCE length.

    """
    with gzip.open(filename, "r") as fastq_file:
        secondlines = islice(fastq_file, 1, 1000, 4)
        temp_length = len(next(secondlines).rstrip())
        for sequence in secondlines:
            read_length = len(sequence.rstrip())
            if temp_length != read_length:
                raise
                logger.exception(
                    f"[ERROR] Sequence length in {filename} is not consistent. Please, trim all "
                    "sequences at the same length.\n"
                    "Exiting the application.\n"
                )
    return read_length


def check_barcodes_lengths(
    read1_length: int, cb_first: int, cb_last: int, umi_first: int, umi_last: int
) -> tuple[slice, slice, int]:
    """Check Read1 length against CELL and UMI barcodes length.

    Args:
        read1_length (int): Read1 length.
        cb_first (int): Barcode first base position for Read1.
        cb_last (int): Barcode last base position for Read1.
        umi_first (int): UMI first base position for Read1.
        umi_last (int): UMI last base position for Read1.

    Returns:
        slice: A `slice` object to extract the Barcode from the sequence string.
        slice: A `slice` object to extract the UMI from the sequence string.
        int: The Barcode + UMI length.

    """
    barcode_length = cb_last - cb_first + 1
    umi_length = umi_last - umi_first + 1
    barcode_umi_length = barcode_length + umi_length
    barcode_slice = slice(cb_first - 1, cb_last)
    umi_slice = slice(umi_first - 1, umi_last)

    if barcode_umi_length > read1_length:
        msg = (
            "[ERROR] Read1 length is shorter than the option you are using for "
            "Cell and UMI barcodes length. Please, check your options and rerun.\n\n"
            "Exiting the application.\n"
        )
        logger.exception(msg)
        raise MismatchedReadlengthsError
    elif barcode_umi_length < read1_length:
        rprint(
            f"[WARNING] Read1 length is {read1_length}bp but you are using {barcode_umi_length}bp for Cell "
            "and UMI barcodes combined.\nThis might lead to wrong cell "
            "attribution and skewed umi counts.\n"
        )

    return (barcode_slice, umi_slice, barcode_umi_length)


def blocks(files: BinaryIO, size: int = 65536) -> Generator[bytes, None, None]:
    """
    A fast way of counting the lines of a large file.
    Ref:
        https://stackoverflow.com/a/9631635/9178565

    Args:
        files (io.handler): A file handler
        size (int): Block size
    Returns:
        A generator
    """
    while True:
        if b := files.read(size):
            yield b
        else:
            break


def get_n_lines(file_path: Path) -> int:
    """
    Determines how many lines have to be processed
    depending on options and number of available lines.
    Checks that the number of lines is a multiple of 4.

    Args:
        file_path (string): Path to a fastq.gz file

    Returns:
        n_lines (int): Number of lines in the file
    """
    rprint("Counting number of reads")
    with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as f:
        n_lines = sum(bl.count("\n") for bl in blocks(f))
    if n_lines % 4 != 0:
        logger.exception(f"{file_path}'s number of lines is not a multiple of 4. The file might be corrupted. Exiting")
        raise NotMultipleofFourError
    return n_lines


def get_read_paths(read1_path: tuple[Path], read2_path: tuple[Path]) -> tuple[tuple[Path], tuple[Path]]:
    """
    Splits up 2 comma-separated strings of input files into list of files
    to process. Ensures both lists are equal in length.

    Args:
        read1_path (string): Comma-separated paths to read1.fq
        read2_path (string): Comma-separated paths to read2.fq
    Returns:
        _read1_path (list(string)): list of paths to read1.fq
        _read2_path (list(string)): list of paths to read2.fq
    """
    # _read1_path = read1_path.split(",")
    # _read2_path = read2_path.split(",")
    if len(read1_path) != len(read2_path):
        logger.exception(
            f"Unequal number of read1 ({len(read1_path)}) and read2({len(read2_path)}) files provided. Exiting"
        )
    return (read1_path, read2_path)
