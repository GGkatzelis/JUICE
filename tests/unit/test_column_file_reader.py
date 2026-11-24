"""column_file_reader.py: the four file layouts (JUICE result, FACSIMILE result, ENZ, old JUICE CSV), the column
names of FACSIMILE headers, and the oddities of FACSIMILE numbers.  live_forms.py: its own self-test under pytest."""
import numpy as np
import pytest

from juice.core import column_file_reader as cf
from juice.core import enz_builder, live_forms

ENZ = """COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00
COLUMN 2=sz
COLUMN 3=jNO2
NUMBER OF COLUMNS=3
&&&&&&&&&&&&&&&&&&&&&
 773812800,  87.332742,  3.70654322E-003,
 773812860,  87.199824,  0,
"""
FACSIMILE_CS = """ * Facsimile-Output / EASY Ver. 6.00
 cs: concentrations and settings
 COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00
 COLUMN 2=OH<0>
 COLUMN 3=JSCALE<0>
 NUMBER OF COLUMNS=3
 &&&&&&&&
 773812860.0      9.995893E-01  1.540849-169
 773812920.0      2.5E+06  1.0
 773812920.0      9.9E+06  9.0
"""
FACSIMILE_PLAIN = """COLUMN 1=TIME
COLUMN 2=OH<0>
 773812860.0  1.5E+05  7.0
 773812920.0  2.5E+05  8.0
 773812920.0  9.9E+05  9.0
"""
JUICE_NEW = """# Calculated by JUICE
COLUMN 1=TIME_s
COLUMN 2=OH
COLUMN 3=HO2
NUMBER OF COLUMNS=3
TIME_s,OH,HO2
773812800,1.5e5,2.0e8
773812860,2.5e5,3.0e8
"""
JUICE_OLD = """# old JUICE output
TIME_s,OH,HO2
773812800,1.5e5,2.0e8
773812860,2.5e5,3.0e8
"""


def written(tmp_path, text, name="f.txt"):
    path = tmp_path / name
    path.write_text(text, encoding="latin-1")
    return str(path)


@pytest.mark.parametrize("text, layout", [(ENZ, "enz"), (FACSIMILE_CS, "enz"), (FACSIMILE_PLAIN, "facsimile"),
                                          (JUICE_NEW, "juice_new"), (JUICE_OLD, "juice_old")])
def test_layout_is_recognised(tmp_path, text, layout):
    # a FACSIMILE .cs with its '&&&&' line is read by the same branch as an ENZ file
    assert cf.detect_format(written(tmp_path, text)) == layout


@pytest.mark.parametrize("raw, index, name", [
    ("TIME / SEC SINCE 1.1.2000 00:00:00", 0, "TIME_s"),
    ("OH<0>", 1, "OH"), ("JSCALE<0>", 5, "JSCALE"), ("sz", 1, "sz"),
    ("p<0> HONO", 2, "P_HONO"), ("d<0> NO", 3, "D_NO"),
    ("k<5,0> !UCH4+OH-->CH3!N", 4, "k_R00005:CH4+OH-->CH3"),
    ("v<12,0> !U NO + O3 --> NO2 !N", 4, "v_R00012:NO + O3 --> NO2"),
])
def test_facsimile_column_names(raw, index, name):
    assert cf._normalize_col_name(raw, index) == name


def test_column_definitions_in_file_order(tmp_path):
    assert cf._parse_column_defs(written(tmp_path, FACSIMILE_CS)) == ["TIME_s", "OH", "JSCALE"]
    assert cf._parse_column_defs(written(tmp_path, JUICE_OLD)) == []


def test_enz_file_with_commas_is_read(tmp_path):
    """Bug of 24 Nov 2026: values were split at blanks only, so no ENZ file (comma-separated, with a comma at the
    end of the line) could be read: 'No data found'."""
    df = cf.load_column_file(written(tmp_path, ENZ))
    assert list(df.columns) == ["sz", "jNO2"] and df.index.name == "TIME_s"
    assert list(df.index) == [773812800.0, 773812860.0]
    assert df["sz"].tolist() == [87.332742, 87.199824] and df["jNO2"].tolist() == [3.70654322e-3, 0.0]


def test_enz_file_written_by_juice_reads_back(tmp_path):
    path = str(tmp_path / "kammer.enz")
    grid = np.array([773812800.0, 773812860.0, 773812920.0])
    enz_builder.write_enz(path, grid, {"hv": np.array([0.0, 1.0, 1.0]), "jNO2": np.array([0.0, 0.0052, 0.0061])})
    df = cf.load_column_file(path)
    assert list(df.index) == list(grid) and df["hv"].tolist() == [0.0, 1.0, 1.0]
    assert df["jNO2"].tolist() == pytest.approx([0.0, 0.0052, 0.0061])


def test_facsimile_result_with_three_digit_exponent_and_repeated_last_time(tmp_path):
    df = cf.load_column_file(written(tmp_path, FACSIMILE_CS))
    assert list(df.columns) == ["OH", "JSCALE"]
    assert df.loc[773812860.0, "JSCALE"] == 1.540849e-169               # FACSIMILE drops the E: 1.540849-169
    assert len(df) == 2 and df.loc[773812920.0, "OH"] == 9.9e6          # of a repeated time the last row is kept


def test_facsimile_file_without_separator_line(tmp_path):
    df = cf.load_column_file(written(tmp_path, FACSIMILE_PLAIN))
    assert list(df.columns) == ["OH", "COL2"]                           # a column without a definition gets a number
    assert len(df) == 2 and df.loc[773812920.0, "OH"] == 2.5e5          # here the first row of a repeated time


def test_juice_layouts(tmp_path):
    for text in (JUICE_NEW, JUICE_OLD):
        df = cf.load_column_file(written(tmp_path, text))
        assert list(df.columns) == ["OH", "HO2"] and list(df.index) == [773812800, 773812860]
        assert df["HO2"].tolist() == [2.0e8, 3.0e8]


def test_file_without_data_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="No data found"):
        cf.load_column_file(written(tmp_path, "COLUMN 1=TIME\nCOLUMN 2=OH<0>\n"))
    with pytest.raises(ValueError, match="No data found"):
        cf.load_column_file(written(tmp_path, "COLUMN 1=TIME\nNUMBER OF COLUMNS=1\n&&&&\n"))


def test_live_forms_self_test():
    live_forms._selftest()
