# Data

The 30 CSV files in `raw/` are frozen **normalised analysis inputs**, not untouched
field records. They were exported from the InteractionExtinctionDebt working
repository (`All/outputs/43 main figures/inputs`) and imported into RALSA on
2026-09-29. SHA-256 checksums in `checksums.csv` identify the exact snapshots.

The dataset collection is associated with Galiana et al. (2024),
*Power laws in species’ biotic interaction networks can be inferred from
co-occurrence data*, https://doi.org/10.1038/s41559-023-02254-y.
The published data/code record is https://doi.org/10.5281/zenodo.8402455.
RALSA uses the local working collection, with the corrections below; it has not
asserted byte-for-byte equivalence to that archive.

## Schema

Each dataset supplies three tables. All identifiers are strings.

| Suffix | Columns | One row represents |
| --- | --- | --- |
| `_interactions.csv` | site, consumer, resource | A realised pair at one site |
| `_cooccurrences.csv` | site, consumer, resource | A pair occurring together at one site |
| `_occupancy.csv` | site, species, trophic_level | One species recorded in one guild at a site |

Counts are binary per site. Duplicate rows are forbidden. Guild values are
`consumer` and `resource`; matching names in different guilds remain distinct.
Co-occurrences include the observed interactions.

## Occurrence definitions

For Quercus, Nahuel, Gottin_HP, Gottin_PP, Garraf_HP, Garraf_PP, Garraf_PP2, Olot and
Montseny, the previous loader took species occurrence from named rows and columns
of local interaction matrices, and co-occurrence from their Cartesian product.
This does not establish independent presence sampling. Interactions were positive cells.

For Salix_Galpar, parasitoid occurrence came from recorded parasitism. Co-occurrences
paired each parasitoid with gallers recorded at its occupied sites, including gallers
without recorded parasitism. Resource occupancy came from those co-occurrence tables.
The species-retention curves inherit these definitions.

## Input decisions

- Empty species names were removed before matrix name repair, preventing the
  fictitious Gottin_HP `V1` consumer. The retained network has 108 regional links.
- Garraf_HP contains 89 links and Olot 92. Earlier comparisons with published
  degree tables differed by one link in each system. Their cause remains unresolved;
  the pipeline preserves the current site records instead of inventing links.
- The original export joined empirical interaction triples into the co-occurrence
  table and removed duplicates.
- `prepare_data.jl` verifies checksums, uniqueness, nonempty identifiers, link
  inclusion and endpoint occupancy, then writes validated copies to `processed/`.
  It does not re-import the original spreadsheets or repeat the earlier R conversion.

Raw inputs should not be edited to change an analysis. A deliberate input revision
must update its provenance and checksums. Data-specific original citations and
redistribution terms require completion before public redistribution; the repository's
MIT license covers code only.
