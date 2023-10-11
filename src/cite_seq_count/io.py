import gzip
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd
from scipy import io, sparse


def write_to_files(
    sparse_matrix: sparse.dok_matrix,
    top_cells: set[str],
    ordered_tags_map: dict[str, int],
    data_type: str,  # should probably be an enum
    outfolder: Path,
) -> None:
    """Write the umi and read sparse matrices to file in gzipped mtx format.

    Args:
        sparse_matrix (dok_matrix): Results in a sparse matrix.
        top_cells (set): Set of cells that are selected for output.
        ordered_tags_map (dict): Tags in order with indexes as values.
        data_type (string): A string definning if the data is umi or read based.
        outfolder (string): Path to the output folder.
    """
    prefix = outfolder.joinpath(f"{data_type}_count").resolve()
    if not prefix.exists():
        prefix.mkdir(exist_ok=True)
    io.mmwrite(prefix.joinpath("matrix.mtx"), sparse_matrix)

    with gzip.open(prefix.joinpath("barcodes.tsv.gz"), "wb") as barcode_file:
        for barcode in top_cells:
            barcode_file.write(f"{barcode}\n".encode())

    with gzip.open(prefix.joinpath("features.tsv.gz"), "wb") as feature_file:
        for feature in ordered_tags_map:
            feature_file.write(f"{feature}\n".encode())

    with open(prefix.joinpath("matrix.mtx"), "rb") as mtx_in:
        with gzip.open(prefix.joinpath("matrix.mtx.gz"), "wb") as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)

    prefix.joinpath("matrix.mtx").unlink()


def write_dense(
    sparse_matrix: sparse.dok_matrix, index: list[str], columns: tuple[str], outfolder: Path, filename: Path
) -> None:
    """
    Writes a dense matrix in a csv format

    Args:
       sparse_matrix (dok_matrix): Results in a sparse matrix.
       index (list): List of TAGS
       columns (set): List of cells
       outfolder (str): Output folder
       filename (str): Filename
    """
    outfolder = outfolder.resolve()
    if not outfolder.exists():
        outfolder.mkdir(exist_ok=True)
    pandas_dense = pd.DataFrame(data=sparse_matrix.todense(), columns=tuple(columns), index=tuple(index))
    pandas_dense.to_csv(path_or_buf=outfolder.joinpath(filename.name), sep="\t")


def write_unmapped(merged_no_match: Counter[str], top_unknowns: int, outfolder: Path, filename: Path) -> None:
    """
    Writes a list of top unmapped sequences

    Args:
        merged_no_match (Counter): Counter of unmapped sequences
        top_unknowns (int): Number of unmapped sequences to output
        outfolder (string): Path of the output folder
        filename (string): Name of the output file
    """

    top_unmapped = merged_no_match.most_common(top_unknowns)
    if not outfolder.exists():
        outfolder.mkdir(exist_ok=True)

    with outfolder.joinpath(filename).open("w") as unknown_file:
        unknown_file.write("tag,count\n")
        for element in top_unmapped:
            unknown_file.write(f"{element[0]},{element[1]}\n")
