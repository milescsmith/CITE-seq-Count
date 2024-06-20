from pathlib import Path

from cite_seq_count.__main__ import main


def test_whole_kit_and_kaboodle():
    module_folder = Path().home().joinpath("workspace", "CITE-seq-Count")
    test_data_folder = module_folder.joinpath("tests", "test_data", "real_data")
    main(
        read1_path_list=[str(test_data_folder.joinpath("test_3_fixed_prot_R1_subset.fq.gz"))],
        read2_path_list=[str(test_data_folder.joinpath("test_3_fixed_prot_R2_subset.fq.gz"))],
        tags=test_data_folder.joinpath("asapseq_ab_panel_v1_citeseq_count_formatted.csv"),
        cb_first=1,
        cb_last=16,
        umi_first=17,
        umi_last=35,
        expected_cells=10000,
        whitelist_file=test_data_folder.joinpath("737K-arc-v1.txt.gz"),
        outfolder=test_data_folder
    )

if __name__ == "__main__":
    test_whole_kit_and_kaboodle()