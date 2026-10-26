# Source-data availability audit

Read `outputs/data_audit_20261026/AUDIT_REPORT.md` in RALSA for conclusions. This audit preserves all source and current analysis files. It does not apply candidate corrections or fit models.

Run these Python scripts in order, with arguments `original_repository RALSA_repository output_directory` for the first two, and `output_directory original_repository RALSA_repository` for the report:

1. `inspect_sources.py`: source inventory, workbook and text-table inspection.
2. `audit_details.py`: site matching, quantitative fields, sample grouping and omitted-field checks.
3. `write_report.py`: evidence-based availability table and audit narrative.

Dependencies: pandas and openpyxl, used only to read original workbooks. No workbook authoring, R execution or source mutation occurs. RDS files are inventoried but not deserialized; their parent CSV and formatter are inspected. Published protocols are documented separately from locally held data.

The report includes expert interpretation of the inspected fields; it is not automatically updated to reflect future source changes without review. In particular, prospective corrections are exported as audit candidates rather than merged into the working inputs.
