from pathlib import Path

import pytest

from cite_seq_count import io


# pytest_plugins = ("pytest_profiling")
@pytest.fixture
def data() -> None:
    from collections import OrderedDict

    from scipy import sparse
    test_matrix = sparse.dok_matrix((4,2))
    test_matrix[1,1] = 1
    pytest.sparse_matrix = test_matrix
    pytest.top_cells = {"ACTGTTTTATTGGCCT","TTCATAAGGTAGGGAT"}
    pytest.ordered_tags_map = OrderedDict({
        "test3-CGTCGTAGCTGATCGTAGCTGAC":0,
        "test2-CGTACGTAGCCTAGC":1,
        "test1-CGTAGCTCG": 3,
        "unmapped": 4
        })
    pytest.data_type = "umi"
    pytest.outfolder = Path("tests/test_data/")

def test_write_to_files(data, tmp_path) -> None:
    import gzip

    import scipy
    io.write_to_files(
        sparse_matrix=pytest.sparse_matrix,
        top_cells=pytest.top_cells,
        ordered_tags_map=pytest.ordered_tags_map,
        data_type=pytest.data_type,
        outfolder=tmp_path
        )
    file = tmp_path.joinpath("umi_count/matrix.mtx.gz")
    with gzip.open(file, "rb") as mtx_file:
        assert isinstance(scipy.io.mmread(mtx_file) ,scipy.sparse.coo_matrix)
