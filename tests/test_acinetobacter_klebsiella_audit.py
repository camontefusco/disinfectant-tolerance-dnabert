from __future__ import annotations

from bs4 import BeautifulSoup

from scripts.audit_acinetobacter_klebsiella_phenotype_genome import (
    extract_ab_2017_biocide_rows,
    table_rows,
)


def test_table_rows_extracts_header_and_values() -> None:
    soup = BeautifulSoup("<table><tr><th>A</th><th>B</th></tr><tr><td>x</td><td>y</td></tr></table>", "html.parser")
    rows = table_rows(soup.find("table"))

    assert rows == [["A", "B"], ["x", "y"]]


def test_extract_ab_2017_biocide_rows_from_minimal_table(tmp_path) -> None:
    raw = tmp_path
    html = """
    <html><body>
    <table></table><table></table><table></table>
    <table>
      <tr><th>Isolate</th><th>Antibiotic susceptibility</th><th>TRI</th><th>CLA</th><th>adeB</th></tr>
      <tr><td>AB01</td><td>S</td><td>16</td><td>32</td><td>0.00</td></tr>
    </table>
    </body></html>
    """
    (raw / "PMC5622949.html").write_text(html)

    rows = extract_ab_2017_biocide_rows(raw)

    assert len(rows) == 2
    assert rows["isolate_id"].nunique() == 1
    assert set(rows["compound"]) == {"triclosan", "chlorhexidine acetate"}
