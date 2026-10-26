from pathlib import Path
import json,csv,sys,hashlib,shutil
import pandas as pd
out=Path(sys.argv[1]);archive=Path(sys.argv[2]);repo=Path(sys.argv[3])
e=json.loads((out/'detailed_evidence.json').read_text(encoding='utf-8'))
rows=[
 ('Quercus',22,'Yes: East/North for 22/22 trees; CRS unidentified','Integer host–enemy counts; Unpar counts','22 matched tree IDs in 2006 and 2007; within-year visits pooled','Unparasitised hosts retained; sampled host-count information; no independent adult enemy abundance','2006 protocol published; comparable 2007 effort unverified','Connectivity','Only 2006 used; weights and coordinates unused'),
 ('Nahuel',8,'Not found','Integer visitation matrices','Only pooled site matrices held; visits described in publication','No independent occupancy or abundance table found','Weekly 5-minute observation protocol published; event effort table absent','No site covariate table found','Numeric taxon IDs; weights unused'),
 ('Salix_Galpar',341,'Yes: 341/341 current sites, 374/374 raw named localities','Gall denominators and enemy counts per rearing lot','Dated lots, 1982–2010; 128 current sites on multiple dates, 114 in multiple years','Hosts recorded without parasitoids; gall counts are sampled abundance, not landscape density','Gall totals per lot; no standard search-effort or complete visit roster established','Elevation, country, region, willow identity, gall type, enemy life-history categories','Pooled across dates; only 341 of 374 raw localities used'),
 ('Gottin_HP',32,'Coordinates not found; nearest-neighbour distance available','Host cells, parasitised cells, P1/P2/P3 counts','nest_ID and spec_ID provide subsamples; no temporal date field','Host records without parasitism; sampled cell counts, not independent adult abundances','8/10/12 trap nests by size class from published protocol; records lack complete empty-unit roster','Area, connectivity, nearest-neighbour distance, land cover, flowers','P1 only in current links; other parasite slots require review'),
 ('Gottin_PP',32,'Coordinates not found; nearest-neighbour distance available','Individual visitor records; pair counts recoverable','Transect 1–12 and T_section A–F; time-code meaning unconfirmed','1,945 records lack plant identity; insect detections available but inclusion rule unresolved; flower cover site-level only','4/8/12 transects × 5 min × 6 rounds published; full site metadata held','Area, connectivity, nearest-neighbour distance, land cover, flowers','Known-plant records match current 1,987 triples; event identifiers unused'),
 ('Garraf_HP',25,'Not found; do not borrow Garraf_PP coordinates','Integer host–enemy matrix counts','One pooled matrix per site','Explicitly retained unattacked hosts; claimed nest/cell columns absent from actual matrices','Trap-nest procedure published; numerical effort table absent','No separate environmental table found','Weights unused; INFO describes columns not present'),
 ('Garraf_PP',40,'Yes: X/Y for 40/40 plots; geographic datum undocumented','Integer visitation matrices','One pooled matrix per plot','No independent species abundance table found','Plot dimensions and visitation protocol published; event-level table absent','Coordinates only in local workbook','Weights and coordinates unused'),
 ('Garraf_PP2',21,'Not found','Binary 0/1 matrices only','One pooled matrix per plot','Four all-zero plant-column presences remain uncertain; independent abundance absent','Plot dimensions/protocol published; event-level table absent','No separate environmental table found','Cannot recover counts from this workbook'),
 ('Olot',14,'Not found','Integer host–enemy matrix counts','One pooled matrix per plot','Explicitly retained unattacked hosts; claimed nest/cell columns absent from actual matrices','Trap-nest procedure published; numerical effort table absent','Taxon code lookup; no environmental table found','Weights unused; INFO describes columns not present'),
 ('Montseny',18,'Not found','Non-integer visitation weights; notes define visits / number of surveys','One pooled matrix per plot; survey counts not separately stored','Taxon lookup includes unused plants; one all-zero site plant row remains uncertain','Rates have effort adjustment, but denominator table absent','Plant/pollinator lookup, focal/unused plant notes; cc/ct codes need key','Cannot treat rates as raw counts or reconstruct visits safely'),
]
columns=['dataset','current_sites','spatial_information','interaction_quantity','within_site_repetition','occurrence_and_abundance','sampling_effort','additional_covariates','current_use_or_issue']
pd.DataFrame(rows,columns=columns).to_csv(out/'dataset_availability.csv',index=False)
issues=[
 ('HIGH','Gottin_HP','Additional parasitoid slots omitted','P1 matches 357 current triples; P2 adds 22 site–pair triples and 2 regional pairs; P3 adds no new triple','Confirm slot semantics and exclusions, then make a separately versioned corrected input','K, L, M and gap may change'),
 ('HIGH','Salix_Galpar','Analysis universe excludes sampled zero-interaction localities','374 named raw localities; all 341 with positive parasitoid records exactly match current sites','Decide whether habitat-loss universe should include the 33 additional sampled localities','K unchanged if only sites added; N and gap change'),
 ('HIGH','Salix_Galpar','Temporal pooling and unequal effort','29 years, multiple lots and dates; positive gall denominators for all 2029 lots','Separate collection amount from persistence; define spatial/temporal sampling unit','Interpretation and potentially K if unit definition changes'),
 ('MEDIUM','Gottin_PP','Missing plant identities','1945 of 7498 rows; 183 insect taxa raw versus 138 in current interaction-based inputs','Clarify whether these are valid presence detections; preserve uncertain scenario','Occupancy and C can change; known interactions unchanged'),
 ('MEDIUM','Quercus','Unused second-year files','22 matching tree IDs have 2007 workbooks; loader reads 2006 only','Confirm temporal comparability and code definitions before year comparison','Separate temporal analysis possible; do not merge years silently'),
 ('MEDIUM','Garraf_HP;Olot','Stale workbook instructions','INFO refers to first-two nest/cell columns; actual matrices begin with named enemies and have no such columns','Request earlier workbook or independent host-count tables','Host presence retained, independent abundance unavailable'),
 ('MEDIUM','Montseny','Rate denominators unavailable','Weights are visits/surveys according to tesaurus; no separate per-site denominator table','Request event totals and survey counts','Raw count and detection models not ready'),
 ('MEDIUM','Garraf_PP2;Montseny','Uncertain zero-only occupancy','Four plant columns and one plant row respectively excluded from confirmed restoration','Seek explicit presence definitions','C changes if restored; K does not unless interactions change'),
 ('MEDIUM','Garraf_HP;Olot','Historical regional-link count mismatch','Current L=89/92 versus earlier published counts 90/93','Keep current records; reconcile taxon identities and versions before publication','No invented replacement links'),
 ('MEDIUM','Quercus;Garraf_PP','Coordinate reference systems undocumented','Projected East/North versus apparently longitude/latitude X/Y','Confirm CRS/datum before mapping or inter-system distances','Distance computations need documented units'),
 ('MEDIUM','Gottin_PP','Replication codes and empty visits','T_section A–F and 4/8/12 transects fit six-round design, but codebook and empty-event roster absent','Confirm code meanings and zero-visit inclusion before detection models','Event-level sampling potential, not yet validated detection histories'),
 ('MEDIUM','Salix_Galpar','Published visit total differs from raw date grouping','783 site×date strings in held file versus 641 visits reported in source paper','Clarify visit definition, date coding and archive version','Do not equate raw distinct dates with standardized surveys'),
]
pd.DataFrame(issues,columns=['priority','dataset','issue','evidence','next_action','potential_effect']).to_csv(out/'issues_and_actions.csv',index=False)
pubs=[
 ('Galiana et al. 2024 supplement','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-023-02254-y/MediaObjects/41559_2023_2254_MOESM1_ESM.pdf','Sampling descriptions; not proof every original field is held locally'),
 ('Grass et al. 2018 original paper/reporting summary','https://www.uni-goettingen.de/de/document/download/634998daf7c35f6a8aa4c9739ef63e2a.pdf/Grass%202018%20Nature%20Ecol%20Evol.pdf','Effort scales with fragment size; six rounds and variable transect/trap-nest counts'),
 ('Kopelke et al. 2017 data paper','https://doi.org/10.1002/ecy.1832','641 site visits over 29 years; dated and georeferenced collections'),
 ('Kaartinen & Roslin 2011','https://doi.org/10.1111/j.1365-2656.2011.01811.x','2006 large-tree protocol; does not by itself establish protocol of held same-ID 2007 files')]
pd.DataFrame(pubs,columns=['source','url','use']).to_csv(out/'publication_evidence.csv',index=False)
report='''# Audit of original ecological data and current RALSA inputs

Audit date: 26 October 2026. Scope: the ten current datasets, their local original-source folders in InteractionExtinctionDebt, source processing scripts, confirmed-restoration records in RALSA, and targeted original-methods references. No original files or analysis inputs were changed. “Not found” means not found in this audited collection, not that investigators never collected it. Internet archives were consulted for documentation, not imported as substitute data.

## Main conclusion

The statement that we have neither locations nor quantitative within-site data is incorrect for the collection as a whole. Coordinates are held for three datasets; nine contain quantitative interaction information, of which eight have counts or recoverable event counts and Montseny has normalized weights. Garraf_PP2 is the only original matrix set limited to binary interactions. Temporal or subsite replication exists in several systems but is not uniformly suitable for detection modelling.

The current RALSA interface deliberately stores binary site–species and site–pair records. Those tables discard counts, dates and spatial metadata. This is a limitation of the current representation, not always of the original data. Independent abundance of BOTH guilds, complete survey histories including empty visits, and standardized effort remain largely missing.

## Quick availability table

|Dataset|Coordinates|Quantitative interactions|Repeat information actually held|
|---|---|---|---|
|Quercus|22/22; projected CRS unidentified|Counts plus unparasitised hosts|2006 and 2007 files for the same 22 tree IDs; comparability needs confirmation|
|Nahuel|Not found|Visit-count matrices|Pooled only|
|Salix_Galpar|341/341 current sites; 374 raw localities|Gall denominators and enemy counts|Dated rearing lots over 1982–2010|
|Gottin_HP|Not found; nearest-neighbour distances held|Host and attacked-cell counts|Nest/specimen IDs; no collection dates|
|Gottin_PP|Not found; nearest-neighbour distances held|Visitor event rows|Transect and A–F section codes; date interpretation unresolved|
|Garraf_HP|Not found|Integer interaction matrices|Pooled only|
|Garraf_PP|40/40; datum unconfirmed|Visit-count matrices|Pooled only|
|Garraf_PP2|Not found|Binary only|Pooled only|
|Olot|Not found|Integer interaction matrices|Pooled only|
|Montseny|Not found|Visits/survey weights, not raw counts|Pooled only; denominators absent|

## Dataset-by-dataset findings

### Quercus

**Held:** 44 workbooks, one for each of 22 tree IDs in 2006 and 2007. Sheet1 contains host rows, enemy columns and an Unpar/unpar column. The metadata has Tree, East, North and Connectivity for all 22 trees. Thousands of trailing empty rows/columns in metadata.csv are formatting/export padding, not additional sites. Only 2006 is loaded by the original processing script and current analysis.

**Usable:** within-year quantitative interaction patterns, sampled unparasitised host counts, spatial separation once the projected coordinate system is confirmed, and a potential paired-year comparison. Unpar and enemy counts should not automatically be summed into a host denominator until the count unit and possible gregarious emergence are verified.

**Missing/unresolved:** within-year dates and survey-specific tables, independent enemy abundance, CRS documentation, and the effort/protocol for the same-ID 2007 files. The published 2011 paper documents three 2006 foliage-sampling occasions for large trees and a different 2007 small-tree component; it does not establish that our held 22 same-ID 2007 files followed an identical protocol. Two years are temporal replication, not repeated detection trials under guaranteed ecological closure. Taxon codes also need a documented lookup before trait matching.

**Source:** Quercus/metadata.csv; Quercus/raw-data/Web{tree} {year} Kaartinen.xlsx; Quercus/data_processing_occurrence.R.

### Nahuel

**Held:** eight 90×14 integer matrices. Values exceed one in every matrix (maximum across matrices 673). The loader transposes them and assigns Pol1…Pol90 and Plant1…Plant14. Its site index is the file order: ag, cl, ll, mh, mnh, qh, qnh, s. One of the 90 potential insect rows has no observed interactions across the retained data; current linked consumers number 89.

**Usable:** reconstruct pooled visit counts and interaction marginals. Marginals reflect recorded visits, not independent abundance.

**Missing:** coordinates, a held species-name key, site covariates, dates, per-visit effort and independent occurrence tables. Weekly sampling is described in the literature, but the held matrices do not preserve those visits. Do not treat all-zero rows as verified local presences.

**Source:** Nahuel/raw-data/vaz_*_matr_f.txt; Nahuel/data_processing_occurrence.R.

### Salix_Galpar

**Held:** 2,029 rearing lots with unique rearing numbers, collection year/date, country, region, locality, decimal-degree latitude/longitude, elevation, willow identity, galler identity, gall type, gall totals, clean/real/parasitised/unparasitised gall fields and natural-enemy counts. Top header rows also supply enemy taxonomy, identification quality and coded life-history categories (P/I, ENDO/ECTO, KOINO/IDIO, and attacked life stage). These provide candidate compatibility-related descriptors, although category definitions and missingness require explicit handling.

All current 341 localities join exactly to raw names, have valid numeric coordinates and have a single consistent coordinate pair within each name. The raw data contain 374 named localities; exactly 341 have positive enemy counts. All 3,104 current site–pair triples match positive raw counts exactly. The other 33 localities have no recorded parasitoid interaction and are absent from the current site universe. They are sampled zero-interaction localities, not automatically missing fieldwork.

Collections span 1982–2010. Among current sites, 128 have multiple date strings and 114 have multiple years. Across all raw localities, 309 of 1,217 site–galler combinations have more than one recorded date. However, grouping raw locality×date strings yields 783 combinations, while the data paper reports 641 visits. The record-to-visit definition therefore needs clarification; do not call all these standardized repeat surveys.

**Usable:** spatial analyses, host-denominator-aware parasitism analyses, effort sensitivity based on collected gall numbers, and temporal sensitivity in subsets. Coarse life-history/gall traits are available. Counts and repeated lots should be kept separate, with gall denominators counted once per lot rather than once per enemy column.

**Missing/unresolved:** standard collection/search effort, complete visits with no collected galls, independent adult enemy abundance and an agreed survey unit. Lots are targeted host collections, not complete censuses of all pairs. Some natural enemies are inquilines or other associates rather than strict parasitoids. Pooling decades can confound spatial support with temporal coverage. The 33 additional sites would change N and removal curves even if K remained unchanged; this audit does not add them.

**Source:** Salix/raw-data/Salix_webs.csv; Salix/lib/format4R.r. RDS files and root csv/df_* tables are derived representations, not independent evidence or additional samples.

### Gottin_HP

**Held:** 3,535 rows with Site, spec_ID, nest_ID, host taxon/name, total cells, parasitised cells and three enemy identity/count slots (P1, P2, P3). One row lacks Site; spec_ID is not a globally unique row key (3,459 unique values). Host records without enemy observations preserve independent host presence and sampled brood-cell abundance.

The original loader explicitly uses P1 only. Its positive named records reproduce all 357 current site–pair triples. P2 holds 64 positive named records, adding 22 distinct site–pair triples and two regional pairs not represented currently: Melittobia.spec–Passaloecus.corniger and Melittobia.spec–Ancistrocerus.trifasciatus. P3 has one positive named record and adds no new triple. P1 also has 23 positive-count records with blank enemy names. Additional slot records are exported for review; their biological interpretation and any intended exclusions must be confirmed before correction.

**Spatial/environmental information:** metadata for the same 32 sites includes area, size class, flower cover/diversity, Hanski connectivity, nearest-neighbour distance, and natural/forest/arable cover at 250/500 scales. No coordinates or full distance matrix were found. A nearest-neighbour distance is not a pairwise distance matrix.

**Replication/effort:** nest IDs identify subsamples, not repeated temporal surveys. The original method used 8/10/12 trap nests according to fragment size. A full record of empty traps is not established from the occupied-host table, so IDs alone are not a complete sampling denominator. No independent adult enemy abundance is held.

**Source:** Gottin/raw-data/host_para_all_interactions.csv; Gottin/metadata.csv; Gottin/data_processing_occurrence.R.

### Gottin_PP

**Held:** 7,498 insect observation rows with Site, Transect, T_section, family, insect name, plant name and site covariates. Of these, 5,553 have plant identities; aggregating these reproduces all 1,987 current site–pair triples. The remaining 1,945 records lack plant identities. Raw insect names span 183 taxa, compared with 138 in the known-link analysis; these extra detections require a documented inclusion decision rather than automatic deletion or restoration.

Transect runs from 1 to 12; each site contains 4, 8 or 12 distinct transect codes. T_section takes A–F at all 32 sites. This is consistent with the published six sampling rounds, but its exact time/section meaning has not been confirmed. There are 1,393 observed site×transect×section groups. Empty sampling events cannot be inferred merely from missing rows.

**Effort correction to our earlier reminder:** the original protocol was 4/8/12 five-minute transects per round in small/medium/large fragments, repeated six times, not four transects at every site. Effort scales with area. Site size classes allow nominal effort to be reconstructed from the protocol after checking the codebook; they do not prove which empty events were completed.

**Usable:** visit counts, relative insect activity indices, habitat covariates, potentially within-site sample comparisons. No species-specific plant abundance table is held; site-level flower cover is not that denominator. The same spatial metadata and coordinate limitations apply as for Gottin_HP.

**Source:** Gottin/raw-data/plant_poll_all_interactions.csv; Gottin/metadata.csv.

### Garraf_HP

**Held:** 25 local host–enemy count matrices, INFO and a metaweb. There are 159 positive local interaction cells, 86 greater than one; values sum to 663. INFO explicitly retains hosts without parasitism, supporting the 115 restored host presences.

**Important mismatch:** INFO says the first two columns contain nest and cell totals, but actual local tables start immediately with named enemy columns. No corresponding total columns, hidden columns or hidden sheets were found. Those host abundance totals are therefore not available merely because INFO mentions them.

**Missing:** coordinates, visit/nest-level records, measured effort and independent enemy abundance. Do not assign Garraf_PP coordinates by similarly numbered plots; no crosswalk establishes shared locations.

**Source:** Garraf-Montseny-Olot/raw-data/garraf-hp.xlsx.

### Garraf_PP

**Held:** 40 local visitation-count matrices, a metaweb and a coordinates sheet. All 40 plot IDs map to coordinates. X/Y values are consistent with longitude/latitude, but the datum is not explicitly stated. Local matrices contain 747 positive cells, of which 240 exceed one; counts sum to 3,576.

**Usable:** spatial distances after documenting the coordinate interpretation, quantitative link frequencies and sampled visitation marginals. **Missing:** visit-level dates, effort denominators and independent plant/insect abundance. Species codes include abbreviated plants and broad insect/morphospecies labels; a complete external trait mapping is not ready.

**Source:** Garraf-Montseny-Olot/raw-data/garraf-pp.xlsx, coordinates and plot1…plot40.

### Garraf_PP2

**Held:** 21 local matrices; all 4,288 positive entries equal one. There are no extra metadata sheets or hidden sheets/columns. This dataset is genuinely binary in the held original file.

**Missing:** coordinates, counts, repeats, independent abundance and effort records. Four named all-zero plant columns remain ambiguous and are excluded from the confirmed occurrence scenario. Request the pre-binarisation data and plot crosswalk before quantitative or spatial extensions.

**Source:** Garraf-Montseny-Olot/raw-data/garraf-pp-2.xlsx, P1…P21.

### Olot

**Held:** 14 local host–enemy count matrices, taxon lookup/INFO and metaweb. There are 253 positive cells, 182 greater than one; counts sum to 1,695. INFO explicitly retains unattacked hosts, supporting 50 restored presences.

As in Garraf_HP, INFO describes nest/cell-total columns that are absent from the held matrices. No hidden sheets/columns recover them. **Missing:** coordinates, event-level repetition, total host abundance and measured effort. Some enemy labels include parasites/predators; some host identities are grouped taxa.

**Source:** Garraf-Montseny-Olot/raw-data/olot.xlsx.

### Montseny

**Held:** 18 local weighted matrices and a taxon/notes sheet. Notes explicitly define interaction strength as visits recorded in a plot divided by number of surveys. Of 1,451 positive cells, 1,439 are noninteger. These are normalized weights, not integer individual counts. Totals should not be interpreted as numbers of insects or visits.

**Usable:** weighted interaction structure under the stated normalization; plant/pollinator name lookup. **Missing:** raw visit counts, survey denominators, dates, coordinates and independent abundance. Do not reverse-engineer survey counts from decimal denominators. cc/ct labels and focal-plant annotations are clues, not a validated environmental treatment table. One all-zero plant row remains an unresolved occupancy addition.

**Source:** Garraf-Montseny-Olot/raw-data/montseny.xlsx, tesaurus and cc/ct sheets.

## What these data can and cannot resolve

1. **Distance:** immediately promising in Salix and Garraf_PP, and in Quercus once CRS is resolved. Coordinates do not measure environmental suitability or dispersal directly. Existing random-removal expectations depend on K and N; spatially structured removal would additionally need arrangement.
2. **Unequal sampled interaction frequency:** addressable in nine systems, but common units cannot be assumed. Counts of visits, reared enemies, attacked cells and normalized visitation rates differ. An interaction matrix's row/column sums are interaction-derived activity, not independent abundance.
3. **Opportunity versus encounter amount:** especially promising in Salix and Gottin_HP because host denominators are held. Quercus offers Unpar counts, with unit verification needed. These still lack independent abundance for both partners.
4. **Detection:** no dataset is automatically ready for a standard repeated-detection model from this audit alone. Salix temporal spacing violates simple closure; Quercus years are ecological changes as well as detection; Göttingen blocks need code meanings and completed-empty-unit records. Reared host subsamples can inform sample coverage without being repeated closed-community surveys.
5. **Compatibility:** life-history descriptors in Salix are useful, but counts alone cannot identify intrinsic affinity. Availability, phenology, sampling, environment and compatibility remain entangled. A targeted subset analysis is more defensible than forcing one uniform model across all ten.

## Existing corrections and provenance

Confirmed occupancy restoration added 629 site–species presences: Quercus 156, Gottin_HP 308, Garraf_HP 115 and Olot 50. Uncertain additions remain separate. Those corrections did not alter interaction K; this audit identifies additional candidate interaction/site-universe issues that could. The fictitious Gottin_HP V1 consumer remains removed; historical one-link differences for Garraf_HP/Olot remain unresolved. The two Göttingen networks share the same 32 sites. Catalan datasets share a provider collection, but no common plot-ID crosswalk proves exact overlap.

The inventory includes 80 initially scanned source/processing files, 49 workbooks and 13 raw text tables, plus additional current-input/derived-document hashes. Workbook inspection covers all sheets, dimensions, nonempty rows, numeric values and comments; metadata sheets and relevant processing transformations were read directly. RDS objects were inventoried, not deserialized; their raw parent CSV and formatter were inspected. No new models were fitted and no data were silently corrected.

## Priority next steps

1. Resolve Gottin_HP P2/P3 semantics and Salix's 33 zero-interaction localities before treating the current input version as final.
2. Start mechanistic work with **Salix**, which combines spatial, dated, host-count and life-history information. First reconcile visit definitions and temporal pooling.
3. Use **Quercus** for a potential paired-year check after confirming 2007 provenance; **Gottin_HP** for host-count/effort work; **Garraf_PP** for spatial checks.
4. Request precise missing pieces: coordinates and plot crosswalks; Göttingen T_section dictionary and empty-event roster; Garraf_HP/Olot host nest/cell totals; Montseny raw visits and survey denominators; Garraf_PP2 pre-binarisation data; Nahuel species/site key and event records. No messages have been sent.

## Files

- dataset_availability.csv: complete ten-dataset availability matrix.
- issues_and_actions.csv: unresolved issues, evidence, action and potential impact.
- file_inventory.csv and additional_input_versions.csv: source paths and SHA-256 hashes.
- workbook_inspection.json and tabular_inspection.json: machine-readable structural inspection.
- detailed_evidence.json: exact joins, counts, metadata excerpts and limitations.
- salix_site_metadata_audit.csv and salix_site_host_replication.csv: site and site–host repeat coverage.
- quercus_site_year_inventory.csv: the 44 tree-year workbooks.
- *_coordinates_as_recorded.csv: untransformed coordinate extractions, not CRS-certified spatial products.
- gottin_*_replication.csv: observed record grouping, not complete survey inventories.
- gottin_HP_additional_slot_records.csv: candidate omitted records for review, not approved corrections.
- matrix_quantitative_inventory.csv: local matrix numeric checks, excluding summary sheets.
- publication_evidence.csv: external method references, distinguished from locally held data.

## Literature consulted

'''
for title,url,use in pubs:report+=f'- [{title}]({url}): {use}.\n'
(out/'AUDIT_REPORT.md').write_text(report,encoding='utf-8')
for file in ['input_changes.csv','presence_record_provenance.csv']:
 p=repo/'outputs/support_prediction_zero'/file
 if p.exists():shutil.copyfile(p,out/('previous_'+file))
for file in ['metadata_relationships.csv','metadata_source_versions.csv']:
 p=repo/'outputs/statistical description of local support/opportunity_theory'/file
 shutil.copyfile(p,out/file)
print('Wrote report, availability matrix, issues and reference table')
