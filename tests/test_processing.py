import copy
import itertools
import random
from collections import Counter
from pathlib import Path

import pytest

from cite_seq_count import preprocessing, processing


def complete_poly_A(seq: str, final_length: int=40) -> str:
    poly_A_len = final_length - len(seq)
    return(seq + "A" * poly_A_len)

def get_sequences(ref_path) -> list[str]:
    with open(ref_path) as adt_ref:
        lines = adt_ref.readlines()
        entries = len(lines) // 2
        sequences = [
            complete_poly_A(lines[i + 1].strip()) for i in range(0, entries, 2)
        ]
    return(sequences)

def extend_seq_pool(ref_seq: str, distance: int) -> list[str]:
    extended_pool = [complete_poly_A(ref_seq)]
    extended_pool.append(modify(ref_seq, distance, modification_type="mutate"))
    extended_pool.append(modify(ref_seq, distance, modification_type="mutate"))
    extended_pool.append(modify(ref_seq, distance, modification_type="mutate"))
    return(extended_pool)

def modify(seq: str, n: int, modification_type: str) -> str:
    bases=list("ATGCN")
    positions = list(range(len(seq)))
    exploded_seq = list(seq)
    for _i in range(n):
        match modification_type:
            case "mutate":
                position = random.choice(positions) # noqa S311
                positions.remove(position)
                temp_bases = copy.copy(bases)
                del temp_bases[bases.index(exploded_seq[position])]
                exploded_seq[position] = random.choice(temp_bases) # noqa S311
            case "delete":
                del exploded_seq[random.randint(0,len(exploded_seq)-2)] # noqa S311
            case "add":
                position = random.randint(0,len(exploded_seq)-1) # noqa S311
                exploded_seq.insert(position, random.choice(bases)) # noqa S311
    return(complete_poly_A("".join(exploded_seq)))

@pytest.fixture
def data() -> None:
    # import json
    from collections import Counter, OrderedDict #x , defaultdict
    # from itertools import islice
    # Test file paths
    pytest.correct_R1_path = "tests/test_data/fastq/correct_R1.fastq.gz"
    pytest.correct_R2_path = "tests/test_data/fastq/correct_R2.fastq.gz"

    pytest.chunk_size = 800
    pytest.tags = OrderedDict({
        "CGTACGTAGCCTAGC": "test2-CGTACGTAGCCTAGC",
        "CGTAGCTCG": "test1-CGTAGCTCG"
        })
    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.indexes = [0,800]
    pytest.correct_whitelist = {"ACTGTTTTATTGGCCT","TTCATAAGGTAGGGAT"}
    pytest.legacy = False
    pytest.debug = False
    pytest.start_trim = 0
    pytest.maximum_distance = 5
    pytest.results = {
        "ACTGTTTTATTGGCCT":
        {"test1-CGTAGCTCG":
        Counter({b"CATTAGTGGT": 3, b"CATTAGTGGG": 2, b"CATTCGTGGT": 1})},
        "TTCATAAGGTAGGGAT":
        {"test2-CGTACGTAGCCTAGC":
        Counter({b"TAGCTTAGTA": 3, b"TAGCTTAGTC": 2, b"GCGATGCATA": 1})}
    }
    pytest.corrected_results = {
        "ACTGTTTTATTGGCCT":
        {"test1-CGTAGCTCG":
        Counter({b"CATTAGTGGT": 6})},
        "TTCATAAGGTAGGGAT":
        {"test2-CGTACGTAGCCTAGC":
        Counter({b"TAGCTTAGTA": 5, b"GCGATGCATA": 1})}
    }
    pytest.umis_per_cell = Counter({
        "ACTGTTTTATTGGCCT": 1,
        "TTCATAAGGTAGGGAT": 2
    })
    pytest.reads_per_cell = Counter({
        "ACTGTTTTATTGGCCT": 3,
        "TTCATAAGGTAGGGAT": 6
    })
    pytest.expected_cells = 2
    pytest.no_match = Counter()
    pytest.collapsing_threshold = 1
    pytest.sliding_window = False
    pytest.max_umis = 20000

    pytest.sequence_pool = []
    pytest.tags_complete = preprocessing.check_tags(preprocessing.parse_tags_csv(Path("tests/test_data/tags/correct.csv")), 5)


@pytest.mark.dependency()
def test_find_best_match_with_1_distance(data) -> None:
    distance = 1
    for tag,name in pytest.tags_complete.items():
        counts = Counter()
        for seq in extend_seq_pool(tag, distance):
            counts[processing.find_best_match(seq, pytest.tags_complete, distance)] += 1
        assert counts[name] == 4 # noqa S101

@pytest.mark.dependency()
def test_find_best_match_with_2_distance(data) -> None:
    distance = 2
    for tag,name in pytest.tags_complete.items():
        counts = Counter()

        for seq in extend_seq_pool(tag, distance):
            counts[processing.find_best_match(seq, pytest.tags_complete, distance)] += 1
        assert counts[name] == 4 # noqa S101

@pytest.mark.dependency()
def test_find_best_match_with_3_distance(data) -> None:
    distance = 3
    for tag,name in pytest.tags_complete.items():
        counts = Counter()

        for seq in extend_seq_pool(tag, distance):
            counts[processing.find_best_match(seq, pytest.tags_complete, distance)] += 1
        assert counts[name] == 4 # noqa S101

@pytest.mark.dependency()
def test_find_best_match_with_3_distance_reverse(data) -> None:
    distance = 3
    for tag,name in sorted(pytest.tags_complete.items()):
        counts = Counter()
        for seq in extend_seq_pool(tag, distance):
            counts[processing.find_best_match(seq, pytest.tags_complete, distance)] += 1
        assert counts[name] == 4 # noqa S101

@pytest.mark.dependency(depends=[
    "test_find_best_match_with_1_distance",
    "test_find_best_match_with_2_distance",
    "test_find_best_match_with_3_distance",
    "test_find_best_match_with_3_distance_reverse",])
def test_classify_reads_multi_process(data) -> None:
    (results, no_match) = processing.map_reads(
        pytest.correct_R1_path,
        pytest.correct_R2_path,
        pytest.tags,
        pytest.barcode_slice,
        pytest.umi_slice,
        pytest.indexes,
        pytest.debug,
        pytest.start_trim,
        pytest.maximum_distance,
        pytest.sliding_window)
    assert len(results) == 2 # noqa S101


@pytest.mark.dependency(depends=["test_classify_reads_multi_process"])
def test_correct_umis(data) -> None:
    temp = processing.correct_umis(pytest.results, 2, pytest.corrected_results.keys(), pytest.max_umis)
    results = temp[0]
    n_corrected = temp[1]
    for cell_barcode in results.keys():
        for TAG in results[cell_barcode]:
            assert len(results[cell_barcode][TAG]) == len(pytest.corrected_results[cell_barcode][TAG])
            assert sum(results[cell_barcode][TAG].values()) == sum(pytest.corrected_results[cell_barcode][TAG].values())
    assert n_corrected == 3 # noqa S101


@pytest.mark.dependency(depends=["test_correct_umis"])
def test_correct_cells(data) -> None:
    processing.correct_cells(pytest.corrected_results, pytest.reads_per_cell, pytest.umis_per_cell, pytest.expected_cells, pytest.collapsing_threshold, pytest.tags)


@pytest.mark.dependency(depends=["test_correct_umis"])
def test_generate_sparse_matrices(data) -> None:
    (umi_results_matrix, read_results_matrix) = processing.generate_sparse_matrices(
        pytest.corrected_results, pytest.ordered_tags_map,
        {"ACTGTTTTATTGGCCT","TTCATAAGGTAGGGAT"}
        )
    assert umi_results_matrix.shape == (4,2) # noqa S101
    assert read_results_matrix.shape == (4,2) # noqa S101
    read_results_matrix = read_results_matrix.tocsr()
    total_reads = sum(
        read_results_matrix[i, j]
        for i, j in itertools.product(
            range(read_results_matrix.shape[0]),
            range(read_results_matrix.shape[1]),
        )
    )
    assert total_reads == 12 # noqa S101
