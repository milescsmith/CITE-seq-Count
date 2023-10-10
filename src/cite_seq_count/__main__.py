#!/usr/bin/env python3.6
"""
Author: Patrick Roelli
"""
import datetime
import os
import sys
import time
from collections import Counter, OrderedDict, defaultdict

# import logging
from pathlib import Path
from typing import Annotated, Optional

import better_exceptions
import typer
from dateutil import tz
from loguru import logger
from multiprocess import Pool, cpu_count
from rich import print as rprint

from cite_seq_count import (
    __version__,
    app,
    io,
    preprocessing,
    processing,
    secondsToText,
    verbosity_level,
    version_callback,
)
from cite_seq_count.logger import init_logger

better_exceptions.hook()

logger.remove()

NUM_CPUS=cpu_count()
DEFAULT_RESULTS_OUTPUT = Path().cwd().joinpath("Results")
DEFAULT_UNMAPPED_OUTPUT = Path().cwd().joinpath("unmapped.csv")
DEFAULT_MULTIPROCESSING_THRESHOLD = 1000001

def create_report(
    outfolder: Path,
    n_reads,
    reads_per_cell,
    no_match,
    start_time,
    umis_corrected,
    bcs_corrected,
    bad_cells,
    bc_threshold,
    umi_threshold,
    read1_path,
    read2_path,
    cb_first,
    cb_last,
    umi_first,
    umi_last,
    expected_cells,
    max_error,
    start_trim,
) -> None:
    """
    Creates a report with details about the run in a yaml format.

    Args:
        n_reads (int): Number of reads that have been processed.
        reads_matrix (scipy.sparse.dok_matrix): A sparse matrix continining read counts.
        no_match (Counter): Counter of unmapped tags.
        version (string): CITE-seq-Count package version.
        start_time (time): Start time of the run.
        args (arg_parse): Arguments provided by the user.

    """
    total_unmapped = sum(no_match.values())
    total_mapped = sum(reads_per_cell.values()) - total_unmapped
    mapped_perc = round((total_mapped / n_reads) * 100)
    unmapped_perc = round((total_unmapped / n_reads) * 100)

    with outfolder.joinpath("run_report.yaml").open("w") as report_file:
        report_file.write(
            f"Date: {datetime.datetime.now(tz=tz.tzlocal()).strftime('%Y-%m-%d')}"
            f"Running time: {secondsToText.secondsToText(time.time() - start_time)}"
            f"CITE-seq-Count Version: {__version__}"
            f"Reads processed: {n_reads}"
            f"Percentage mapped: {mapped_perc}"
            f"Percentage unmapped: {unmapped_perc}"
            f"Uncorrected cells: {len(bad_cells)}"
            f"Correction:"
            f"\tCell barcodes collapsing threshold: {bc_threshold}"
            f"\tCell barcodes corrected: {bcs_corrected}"
            f"\tUMI collapsing threshold: {umi_threshold}"
            f"\tUMIs corrected: {umis_corrected}"
            f"Run parameters:"
            f"\tRead1_paths: {read1_path}"
            f"\tRead2_paths: {read2_path}"
            f"\tCell barcode:"
            f"\t\tFirst position: {cb_first}"
            f"\t\tLast position: {cb_last}"
            f"\tUMI barcode:"
            f"\t\tFirst position: {umi_first}"
            f"\t\tLast position: {umi_last}"
            f"\tExpected cells: {expected_cells}"
            f"\tTags max errors: {max_error}"
            f"\tStart trim: {start_trim}"
        )


@app.callback(invoke_without_command=True)
@app.command(
    # name="count",
    no_args_is_help=True,
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
    )
def main(
    read1_path: Annotated[
        Path,
        typer.Option(
            "-R1",
            "--read1",
            help=(
                "The path of Read1 in gz format, or a comma-separated list of paths to all Read1 files in "
                "gz format (E.g. A1.fq.gz, B1.fq,gz, ..."
            ),
            rich_help_panel = "Inputs",
            file_okay=True,
            resolve_path=True,
            dir_okay=False,
            readable=True,
        )
    ],
    read2_path: Annotated[
        Path,
        typer.Option(
            "-R2",
            "--read2",
            help=(
                "The path of Read2 in gz format, or a comma-separated list of paths to all Read2 files in "
                "gz format (E.g. A2.fq.gz, B2.fq,gz, ..."
            ),
            rich_help_panel = "Inputs",
            file_okay=True,
            resolve_path=True,
            dir_okay=False,
            readable=True,
        )
    ],
    tags: Annotated[
        Path,
        typer.Option(
            "-t",
            "--tags",
            help=(
                "The path to the csv file containing the antibody\n"
                "barcodes as well as their respective names.\n\n"
                "Example of an antibody barcode file structure:\n\n"
                "\tATGCGA,First_tag_name"
                "\n\n\tGTCATG,Second_tag_name"
            ),
            rich_help_panel = "Inputs",
            file_okay=True,
            resolve_path=True,
            dir_okay=False,
            readable=True,
        )
    ],

    cb_first: Annotated[
        int,
        typer.Option(
            "-cbf",
            "--cell_barcode_first_base",
            help="Postion of the first base of your cell barcodes.",
            rich_help_panel = "Barcodes",
        )
    ],
    cb_last: Annotated[
        int,
        typer.Option(
            "-cbl",
            "--cell_barcode_last_base",
            help=("Postion of the last base of your cell barcodes."),
            rich_help_panel = "Barcodes",
        )
    ],
    umi_first: Annotated[
        int,
        typer.Option(
            "-umif",
            "--umi_first_base",
            help="Postion of the first base of your UMI.",
            rich_help_panel = "Barcodes",
        )
    ],
    umi_last: Annotated[
        int,
        typer.Option(
            "-umil",
            "--umi_last_base",
            help="Postion of the last base of your UMI.",
            rich_help_panel = "Barcodes",
        )
    ],
    umi_threshold: Annotated[
        int,
        typer.Option(
            "--umi_collapsing_dist",
            help="threshold for umi collapsing.",
            rich_help_panel = "Barcodes",
        )
    ] = 2,
    bc_threshold: Annotated[
        int,
        typer.Option(
            "--bc_collapsing_dist",
            help="threshold for cellular barcode collapsing.",
            rich_help_panel = "Barcodes",
        )
    ] = 1,
    # cells = parser.add_argument_group(
    #     "Cells", description=("Expected number of cells and potential whitelist")
    # ),

    expected_cells: Annotated[
        int,
        typer.Option(
            "-cells",
            "--expected_cells",
            help=("Number of expected cells from your run."),
            rich_help_panel = "Cells",
        )
    ] = 0,
    whitelist: Annotated[
        Optional[Path], # noqa UP007
        typer.Option(
            "-wl",
            "--whitelist",
            help=(
                "A csv file containning a whitelist of barcodes produced"
                " by the mRNA data.\n\n"
                "\tExample:\n\n"
                "\t  ATGCTAGTGCTA\n\n"
                "\t  GCTAGTCAGGAT\n\n"
                "\t  CGACTGCTAACG\n\n\n\n"
                "Or 10X-style:\n\n"
                "\t  ATGCTAGTGCTA-1\n\n"
                "\t  GCTAGTCAGGAT-1\n\n"
                "\t  CGACTGCTAACG-1\n"
            ),
            rich_help_panel = "Cells",
            file_okay=True,
            resolve_path=True,
            dir_okay=False,
            readable=True,
        )
    ] = None,

    # FILTERS group.
    # filters = parser.add_argument_group(
    #     "TAG filters", description=("Filtering and trimming for read2.")
    # ),
    max_error: Annotated[
        int,
        typer.Option(
            "--max-errors",
            help=("Maximum Levenshtein distance allowed for antibody barcodes."),
            rich_help_panel = "Filters",
        )
    ] = 2,
    start_trim: Annotated[
        int,
        typer.Option(
            "-trim",
            "--start-trim",
            help=("Number of bases to discard from read2."),
            rich_help_panel = "Filters",
        )
    ] = 0,
    # Remaining arguments.
    n_threads: Annotated[
        int,
        typer.Option(
            "-T",
            "--threads",
            help="How many threads are to be used for running the program",
        )
    ] = NUM_CPUS,
    first_n: Annotated[
        Optional[int], # noqa UP007
        typer.Option(
            "-n",
            "--first_n",
            help="Select N reads to run on instead of all.",
        )
    ] = None,
    outfolder: Annotated[
        Path,
        typer.Option(
            "-o",
            "--output",
            help="Results will be written to this folder",
        )
    ] = DEFAULT_RESULTS_OUTPUT,
    unmapped_file: Annotated[
        Path,
        typer.Option(
            "-u",
            "--unmapped-tags",
            help="Write table of unknown TAGs to file.",
        )
    ] = DEFAULT_UNMAPPED_OUTPUT,
    unknowns_top: Annotated[
        int,
        typer.Option(
            "-ut",
            "--unknown-top-tags",
            help="Top n unmapped TAGs.",
        )
    ] = 100,
    *,
    no_umi_correction: Annotated[
        bool,
        typer.Option(
            "--no_umi_correction",
            help="Deactivate UMI collapsing",
            rich_help_panel = "Barcodes",
        )
    ] = False,
    sliding_window: Annotated[
        bool,
        typer.Option(
            "--sliding-window",
            help=("Allow for a sliding window when aligning."),
            rich_help_panel = "Filters",
        )
    ] = False,
    dense: Annotated[
        bool,
        typer.Option(
            "--dense",
            help="Add a dense output to the results folder",
        )
    ] = False,
    debug: Annotated[
        bool,
        typer.Option(
            "--debug",
            help="Print extra information for debugging."
        )
    ] = False,
    version: Annotated[ # noqa ARG001
        bool,
        typer.Option(
            "--version",
            callback=version_callback,
            help="Print version number.",
        )
    ] = False,
) -> None:
    # Create logger and stream handler
    init_logger(verbosity_level)

    start_time = time.time()


    if whitelist:
        rprint("Loading whitelist")
        (whitelist, bc_threshold) = preprocessing.parse_whitelist_csv(
            filename=whitelist,
            barcode_length=cb_last - cb_first + 1,
            collapsing_threshold=bc_threshold,
        )
    else:
        whitelist = False

    # Load TAGs/ABs.
    ab_map = preprocessing.parse_tags_csv(tags)
    ab_map = preprocessing.check_tags(ab_map, max_error)

    # Identify input file(s)
    read1_paths, read2_paths = preprocessing.get_read_paths(
        read1_path, read2_path
    )

    # preprocessing and processing occur in separate loops so the program can crash earlier if
    # one of the inputs is not valid.
    read1_lengths = []
    read2_lengths = []
    for read1_path, read2_path in zip(read1_paths, read2_paths, strict=True):
        # Get reads length. So far, there is no validation for Read2.
        read1_lengths.append(preprocessing.get_read_length(read1_path))
        read2_lengths.append(preprocessing.get_read_length(read2_path))
        # Check Read1 length against CELL and UMI barcodes length.
        (
            barcode_slice,
            umi_slice,
            barcode_umi_length,
        ) = preprocessing.check_barcodes_lengths(
            read1_lengths[-1],
            cb_first,
            cb_last,
            umi_first,
            umi_last,
        )
    # Ensure all files have the same input length
    # if len(set(read1_lengths)) != 1:
    #    sys.exit('Input barcode fastqs (read1) do not all have same length.\nExiting')

    # Initialize the counts dicts that will be generated from each input fastq pair
    final_results = defaultdict(lambda: defaultdict(Counter))
    umis_per_cell = Counter()
    reads_per_cell = Counter()
    merged_no_match = Counter()
    number_of_samples = len(read1_paths)
    n_reads = 0

    # Print a statement if multiple files are run.
    if number_of_samples != 1:
        rprint(f"Detected {number_of_samples} files to run on.")

    for read1_path, read2_path in zip(read1_paths, read2_paths, strict=True):
        if first_n:
            n_lines = (first_n * 4) / number_of_samples
        else:
            n_lines = preprocessing.get_n_lines(read1_path)
        n_reads += int(n_lines / 4)
        rprint("Started mapping")
        rprint(f"Processing {n_reads:,} reads")
        # Run with one process
        if n_threads <= 1 or n_reads < DEFAULT_MULTIPROCESSING_THRESHOLD:
            rprint("CITE-seq-Count is running with one core.")
            (_final_results, _merged_no_match) = processing.map_reads(
                read1_path=read1_path,
                read2_path=read2_path,
                tags=ab_map,
                barcode_slice=barcode_slice,
                umi_slice=umi_slice,
                indexes=[0, n_reads],
                whitelist=whitelist,
                debug=debug,
                start_trim=start_trim,
                maximum_distance=max_error,
                sliding_window=sliding_window,
            )
            rprint("Mapping done")
            _umis_per_cell = Counter()
            _reads_per_cell = Counter()
            for cell_barcode, counts in _final_results.items():
                _umis_per_cell[cell_barcode] = sum(len(counts[UMI]) for UMI in counts)
                _reads_per_cell[cell_barcode] = sum(
                    sum(counts[UMI].values()) for UMI in counts
                )
        else:
            # Run with multiple processes
            rprint(f"CITE-seq-Count is running with {n_threads} cores.")
            p = Pool(processes=n_threads)
            chunk_indexes = preprocessing.chunk_reads(n_reads, n_threads)
            parallel_results = []

            for indexes in chunk_indexes:
                p.apply_async(
                    processing.map_reads,
                    args=(
                        read1_path,
                        read2_path,
                        ab_map,
                        barcode_slice,
                        umi_slice,
                        indexes,
                        whitelist,
                        debug,
                        start_trim,
                        max_error,
                        sliding_window,
                    ),
                    callback=parallel_results.append,
                    error_callback=sys.stderr,
                )
            p.close()
            p.join()
            rprint("Mapping done")
            rprint("Merging results")

            (
                _final_results,
                _umis_per_cell,
                _reads_per_cell,
                _merged_no_match,
            ) = processing.merge_results(parallel_results=parallel_results)
            del parallel_results

        # Update the overall counts dicts
        umis_per_cell |= _umis_per_cell
        reads_per_cell |= _reads_per_cell
        merged_no_match |= _merged_no_match
        for cell_barcode in _final_results:
            for tag in _final_results[cell_barcode]:
                if tag in final_results[cell_barcode]:
                    # Counter + Counter = Counter
                    final_results[cell_barcode][tag] += _final_results[cell_barcode][
                        tag
                    ]
                else:
                    # Explicitly save the counter to that tag
                    final_results[cell_barcode][tag] = _final_results[cell_barcode][tag]
    ordered_tags_map = OrderedDict()
    for i, tag in enumerate(ab_map.values()):
        ordered_tags_map[tag] = i
    ordered_tags_map["unmapped"] = i + 1

    # Correct cell barcodes
    if bc_threshold > 0:
        if len(umis_per_cell) <= expected_cells:
            rprint(
                f"Number of expected cells, {expected_cells}, is higher "
                f"than number of cells found {len(umis_per_cell)}.\nNot performing"
                "cell barcode correction"
                ""
            )
            bcs_corrected = 0
        else:
            rprint("Correcting cell barcodes")
            if not whitelist:
                (
                    final_results,
                    umis_per_cell,
                    bcs_corrected,
                ) = processing.correct_cells(
                    final_results=final_results,
                    reads_per_cell=reads_per_cell,
                    umis_per_cell=umis_per_cell,
                    expected_cells=expected_cells,
                    collapsing_threshold=bc_threshold,
                    ab_map=ordered_tags_map,
                )
            else:
                (
                    final_results,
                    umis_per_cell,
                    bcs_corrected,
                ) = processing.correct_cells_whitelist(
                    final_results=final_results,
                    umis_per_cell=umis_per_cell,
                    whitelist=whitelist,
                    collapsing_threshold=bc_threshold,
                    ab_map=ordered_tags_map,
                )
    else:
        bcs_corrected = 0

    # If given, use whitelist for top cells
    if whitelist:
        top_cells = whitelist
        for missing_cell in top_cells:
            if missing_cell in final_results:
                continue
            final_results[missing_cell] = {}
            for tag in ordered_tags_map:
                final_results[missing_cell][tag] = Counter()
            top_cells.add(missing_cell)
    else:
        # Select top cells based on total umis per cell
        top_cells_tuple = umis_per_cell.most_common(expected_cells)
        top_cells = {pair[0] for pair in top_cells_tuple}

    # UMI correction

    if no_umi_correction:
        # Don't correct
        umis_corrected = 0
        aberrant_cells = set()
    else:
        # Correct UMIS
        (final_results, umis_corrected, aberrant_cells) = processing.correct_umis(
            final_results=final_results,
            collapsing_threshold=umi_threshold,
            top_cells=top_cells,
            max_umis=20000,
        )

    # Remove aberrant cells from the top cells
    for cell_barcode in aberrant_cells:
        top_cells.remove(cell_barcode)

    # Create sparse aberrant cells matrix
    (umi_aberrant_matrix, read_aberrant_matrix) = processing.generate_sparse_matrices(
        final_results=final_results,
        ordered_tags_map=ordered_tags_map,
        top_cells=aberrant_cells,
    )

    # Write uncorrected cells to dense output
    io.write_dense(
        sparse_matrix=umi_aberrant_matrix,
        index=list(ordered_tags_map.keys()),
        columns=aberrant_cells,
        outfolder=os.path.join(outfolder, "uncorrected_cells"),
        filename="dense_umis.tsv",
    )

    # Create sparse matrices for results
    (umi_results_matrix, read_results_matrix) = processing.generate_sparse_matrices(
        final_results=final_results,
        ordered_tags_map=ordered_tags_map,
        top_cells=top_cells,
    )

    # Write umis to file
    io.write_to_files(
        sparse_matrix=umi_results_matrix,
        top_cells=top_cells,
        ordered_tags_map=ordered_tags_map,
        data_type="umi",
        outfolder=outfolder,
    )

    # Write reads to file
    io.write_to_files(
        sparse_matrix=read_results_matrix,
        top_cells=top_cells,
        ordered_tags_map=ordered_tags_map,
        data_type="read",
        outfolder=outfolder,
    )

    # Write unmapped sequences
    io.write_unmapped(
        merged_no_match=merged_no_match,
        top_unknowns=unknowns_top,
        outfolder=outfolder,
        filename=unmapped_file,
    )

    # Create report and write it to disk
    create_report(
        outfolder=outfolder,
        n_reads=n_reads,
        reads_per_cell=reads_per_cell,
        no_match=merged_no_match,
        start_time=start_time,
        umis_corrected=umis_corrected,
        bcs_corrected=bcs_corrected,
        bad_cells=aberrant_cells,
        bc_threshold=bc_threshold,
        umi_threshold=umi_threshold,
        read1_path=read1_path,
        read2_path=read2_path,
        cb_first=cb_first,
        cb_last=cb_last,
        umi_first=umi_first,
        umi_last=umi_last,
        expected_cells=expected_cells,
        max_error=max_error,
        start_trim=start_trim,
    )

    # Write dense matrix to disk if requested
    if dense:
        rprint("Writing dense format output")
        io.write_dense(
            sparse_matrix=umi_results_matrix,
            index=list(ordered_tags_map.keys()),
            columns=top_cells,
            outfolder=outfolder,
            filename="dense_umis.tsv",
        )
