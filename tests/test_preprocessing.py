import importlib.resources as ir
import io
from pathlib import Path

import pytest

from cite_seq_count import preprocessing


@pytest.fixture
def data() -> None:
    from collections import OrderedDict
    from itertools import islice

    # Test file paths
    pytest.correct_whitelist_path = ir.files("tests").joinpath("test_data", "whitelists", "correct.csv")
    pytest.correct_tags_path = ir.files("tests").joinpath("test_data", "tags", "correct.csv")
    pytest.correct_R1_path = ir.files("tests").joinpath("test_data", "fastq", "correct_R1.fastq.gz")
    pytest.correct_R2_path = ir.files("tests").joinpath("test_data", "fastq", "correct_R2.fastq.gz")
    pytest.corrupt_R1_path = ir.files("tests").joinpath("test_data", "fastq", "corrupted_R1.fastq.gz")
    pytest.corrupt_R2_path = ir.files("tests").joinpath("test_data", "fastq", "corrupted_R2.fastq.gz")

    pytest.correct_R1_multipath = "path/to/R1_1.fastq.gz,path/to/R1_2.fastq.gz"
    pytest.correct_R2_multipath = "path/to/R2_1.fastq.gz,path/to/R2_2.fastq.gz"
    pytest.incorrect_R2_multipath = "path/to/R2_1.fastq.gz,path/to/R2_2.fastq.gz,path/to/R2_3.fastq.gz"

    pytest.correct_multipath_result = (
        [Path("path/to/R1_1.fastq.gz"), Path("path/to/R1_2.fastq.gz")],
        [Path("path/to/R2_1.fastq.gz"), Path("path/to/R2_2.fastq.gz")]
        )

    # Create some variables to compare to
    pytest.correct_whitelist = {"ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"}
    pytest.correct_tags = {
        "AGGACCATCCAA":"CITE_LEN_12_1",
        "ACATGTTACCGT":"CITE_LEN_12_2",
        "AGCTTACTATCC":"CITE_LEN_12_3",
        "TCGATAATGCGAGTACAA":"CITE_LEN_18_1",
        "GAGGCTGAGCTAGCTAGT":"CITE_LEN_18_2",
        "GGCTGATGCTGACTGCTA":"CITE_LEN_18_3",
        "TGTGACGTATTGCTAGCTAG":"CITE_LEN_20_1",
        "ACTGTCTAACGGGTCAGTGC":"CITE_LEN_20_2",
        "TATCACATCGGTGGATCCAT":"CITE_LEN_20_3"}
    pytest.correct_ordered_tags = OrderedDict({
        "TGTGACGTATTGCTAGCTAG":"CITE_LEN_20_1-TGTGACGTATTGCTAGCTAG",
        "ACTGTCTAACGGGTCAGTGC":"CITE_LEN_20_2-ACTGTCTAACGGGTCAGTGC",
        "TATCACATCGGTGGATCCAT":"CITE_LEN_20_3-TATCACATCGGTGGATCCAT",
        "TCGATAATGCGAGTACAA":"CITE_LEN_18_1-TCGATAATGCGAGTACAA",
        "GAGGCTGAGCTAGCTAGT":"CITE_LEN_18_2-GAGGCTGAGCTAGCTAGT",
        "GGCTGATGCTGACTGCTA":"CITE_LEN_18_3-GGCTGATGCTGACTGCTA",
        "AGGACCATCCAA":"CITE_LEN_12_1-AGGACCATCCAA",
        "ACATGTTACCGT":"CITE_LEN_12_2-ACATGTTACCGT",
        "AGCTTACTATCC":"CITE_LEN_12_3-AGCTTACTATCC"})
    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.barcode_umi_length = 26

@pytest.mark.dependency()
def test_parse_whitelist_csv(data: None) -> None:
    assert preprocessing.parse_whitelist_csv(pytest.correct_whitelist_path, 16, 1) == (pytest.correct_whitelist,1) # noqa S101

@pytest.mark.dependency()
def test_parse_tags_csv(data: None) -> None:
    assert preprocessing.parse_tags_csv(pytest.correct_tags_path) == pytest.correct_tags # noqa S101

@pytest.mark.dependency(depends=["test_parse_tags_csv"])
def test_check_tags(data: None) -> None:
    assert preprocessing.check_tags(pytest.correct_tags, 5) == pytest.correct_ordered_tags # noqa S101

@pytest.mark.dependency(depends=["test_check_tags"])
def test_check_distance_too_big_between_tags(data: None) -> None:
    with pytest.raises(SystemExit):
        preprocessing.check_tags(pytest.correct_tags, 8)

@pytest.mark.dependency(depends=["test_parse_whitelist_csv"])
def test_check_barcodes_lengths(data: None) -> None:
    assert preprocessing.check_barcodes_lengths(26, 1, 16, 17, 26) == (pytest.barcode_slice, pytest.umi_slice, pytest.barcode_umi_length) # noqa S101

@pytest.mark.dependency()
def test_get_n_lines(data: None):
  assert preprocessing.get_n_lines(pytest.correct_R1_path) == (200 * 4) # noqa S101

@pytest.mark.dependency(depends=["test_get_n_lines"])
def test_get_n_lines_not_multiple_of_4(data: None) -> None:
  with pytest.raises(preprocessing.NotMultipleofFourError):
    preprocessing.get_n_lines(pytest.corrupt_R1_path)

# @pytest.mark.dependency()
# def test_corrrect_multipath(data):
#   assert preprocessing.get_read_paths(pytest.correct_R1_multipath, pytest.correct_R2_multipath) == pytest.correct_multipath_result

# @pytest.mark.dependency(depends=["test_get_n_lines"])
# def test_incorrrect_multipath(data):
#   with pytest.raises(SystemExit):
#     preprocessing.get_read_paths(pytest.correct_R1_multipath, pytest.incorrect_R2_multipath)
