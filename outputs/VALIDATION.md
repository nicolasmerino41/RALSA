# Migration validation

The RALSA results were compared with the existing manuscript analysis on 2026-09-29.
Rows were matched by dataset, species identity, group and removal level as applicable.
Numerical tolerance: relative 1e-11, absolute 1e-12. All comparisons passed.

```
link_support: 4304 rows matched; maximum absolute difference 0
retention_and_states: 90 rows matched; maximum absolute difference 1.31e-14
conversion_and_pink_area: 10 rows matched; maximum absolute difference 0
pooled_retention: 9 rows matched; maximum absolute difference 1.67e-15
degree_example: 240 rows matched; maximum absolute difference 0
support_by_guild: 40 rows matched; maximum absolute difference 0
support_merged: 20 rows matched; maximum absolute difference 0
group_definitions: 40 rows matched; maximum absolute difference 0
group_curves: 180 rows matched; maximum absolute difference 0
```

The complete pipeline validates the input snapshots and regenerates all main tables
and figures. Figure numbering and layout were changed; analysis definitions were retained.
The reference was `All/outputs/43 main figures/tables` in the earlier working repository.
