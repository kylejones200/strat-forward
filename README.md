# strat-forward

A two-dimensional forward stratigraphic model of a dip section. It covers accommodation (eustasy and
tectonic subsidence), conservative diffusive sediment transport with separate grain classes, Athy
compaction and optional flexural isostasy. Every deposit, erosion and active-layer exchange is kept in an
event ledger, and the preserved stratigraphic record is built from that ledger rather than from sampled
surfaces.

The model is a numerical research instrument. It has been verified numerically (mass balance, convergence,
analytic plate solutions); it has not been validated against flume or field data.

## Build, test, run

```bash
cargo test --release
cargo run --release -- example-config > my.json
cargo run --release -- run --out out/my_run --config my.json [--deck] [--run-id UUID] [--openlineage-url URL]
```

The reference configuration is `examples/prograding_wedge.json`: a 100 km section of 200 cells, run for
2 Myr with 100 record steps of 20 kyr, sand and mud supplied from the landward edge, and a 0.5 m unit
thickness cap (`"max_unit_solid_m": null` turns the cap off). `examples/prograding_wedge_flexure.json` adds
a 20 km elastic plate.

```bash
cargo run --release -- run --out out/prograding_wedge --config examples/prograding_wedge.json --deck
cargo run --release --example sensitivity -- out/sensitivity grid
cargo run --release --example deck_variants -- out/deck_variants
```

## Model

- **Accommodation.** Sea level is a base level plus a sinusoid and a trend. Tectonic subsidence varies
  linearly from proximal to distal and is read as water-loaded subsidence.
- **Transport.** Explicit finite-volume diffusion of the surface: `K_land` above sea level and
  `K_marine·exp(−depth/λ) + K_deep` below, blended by tanh over `shoreline_width_m` (0 is a sharp switch).
  Each grain class has its own multipliers. Outgoing sediment has the composition of a well-mixed active
  layer, and each class's flux is limited to its stock, so the solid budget closes to rounding.
- **Compaction.** Athy φ(z) = φ0·exp(−z/c), solved per unit by Newton iteration.
- **Flexural isostasy** (optional). Thin elastic plate, `D w'''' + Δρ g w = q`, with an energy-based bending
  operator and a pentadiagonal LDLᵀ solve every substep. Elastic thickness 0 gives local (Airy) isostasy.
- **Event ledger and record.** Replaying the ledger reproduces every surviving unit. Each unit records its
  base contact (basement, conformable, diastem or erosional).

## Outputs

Each run writes into a staging directory and replaces `--out` as a whole only on success; the replaced
outputs are kept as `<out>.previous` until the next successful run.

- `summary.json`: configuration, mass budget, surface-reconstruction gap and its breakdown, ledger replay error
- `ledger.json`, `stratigraphy.json`, `hiatus_surfaces.json`, `wheeler.json`, `trajectory.json`,
  `surface_grid.json`, `wells.json`, `section.svg`
- `prov.jsonld`: W3C PROV-O provenance of the record (one activity per ledger entry, one entity per unit)
- `openlineage-events.json`: OpenLineage 2-0-2 run events
- `SECTION.DATA`, `deck_summary.json` (with `--deck`): an Eclipse-format corner-point grid of the sandiest
  10 km window, re-compacted to the datum depth, with Kozeny–Carman permeability

The `strat:` vocabulary in `prov.jsonld` (`urn:strat-forward:vocab#`) is local to this model: `run`, `cell`,
`xCenterM`, `recordStep`, `substeps`, `bypassM`, `deposited_<class>_m`, `removed_<class>_m`,
`surviving_<class>_m`, `baseContact`, `cutM`, `missingSteps`, `truncates`, `initialSediment`,
`recycledFraction`, `meanEntryAgeMa`, `meanTransitMyr`. Geological ages are carried as `strat:startAgeMa`
and `strat:endAgeMa`.

## Tests

- `tests/numerics.rs`: mass balance, diffusion decay against the analytic rate, surface-reconstruction gap,
  shoreline blending, ledger rules
- `tests/isostasy.rs`: plate solver against Airy and the analytic line-load profile, force and moment
  balance, shoreline continuity of the load
- `tests/provenance.rs`: origin-tracer totals and PROV document structure
- `tests/openlineage.rs`: events from real runs validated against the OpenLineage schemas
- `tests/rundir.rs`: atomic replacement of the output directory and reader behavior

## Limits

- Explicit time stepping; the substep is set by the largest diffusivity.
- Diffusion only: no channels, turbidity currents, carbonate production or third dimension.
- Facies are water-depth bins (0, 15 and 80 m), a labeling rule.
- The reservoir export's rock model (re-compaction and Kozeny–Carman) is uncalibrated.
- No comparison against published flume or field cases.

## License

MIT. The OpenLineage schemas in `tests/fixtures/openlineage` are Apache License 2.0.
