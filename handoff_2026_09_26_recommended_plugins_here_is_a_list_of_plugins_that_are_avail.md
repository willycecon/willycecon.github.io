# Handoff: STEM OPT Leaver Figures and `_staying` Extensive Margin

**Date:** 2026-09-25  
**Status:** Code edits completed; HPCC `_staying` jobs finished and GraphData synced; `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd` statically verified but **not knitted**; `Draft/Paper.tex` integration incomplete. Session stopped mid-task.

## Goal

Align the STEM OPT leaver figures with the current policy estimand—U.S. retention/staying rather than leaving—and isolate the complemented extensive-margin outcome. Specifically:

- Reframe §3C of `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd` around Kaplan–Meier U.S. retention \(S(h)\), not the inverse exit hazard.
- Remove the deprecated ACS-based over-education (`over_edu`) outcome from intensive-margin facet graphs while retaining the O*NET over-education outcome (`onet_over_edu`).
- Add a mechanical complemented extensive-margin outcome \(1 - y_{\text{abroad}}\), tagged `_staying`, preserving all legacy abroad-first results.
- Assess existing STEM×OPT DDD/event-study evidence and identify what remains undocumented in §4.3.

## Motivation

The project’s STEM×OPT estimand is now U.S. retention—“still U.S.-resident at horizon \(h\)”—not leaving. Mislabeling complements can invert the sign of a share-DD and conflate first-job U.S. placement with horizon-based retention. The figures must match the estimator before any draft language or causal interpretation is advanced.

## Dependencies

| Dependency | Path / Name | Status | Notes |
|---|---|---:|---|
| Leaver graph R Markdown | `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd` | BUILT | Edited §3C and intensive/extensive-margin readers; static checks pass; not knitted. |
| Changelog | `CHANGELOG.md` | BUILT | Entries at lines 6, 10, 13. |
| Pooled hiring-margin producer | `Codes/hpcc_ana/A_JMP4_refactor/grad/JMP4_grad_t_leavers_stemxopt_hiringmargin.R` | BUILT | `EXTMARGIN_OUTCOME=abroad|stayer` switch added around line 394. |
| Event-study producer | `Codes/hpcc_ana/A_JMP4_refactor/grad/JMP4_grad_t_leavers_stemxopt_extmargin_eventstudy.R` | BUILT | `EXTMARGIN_OUTCOME=abroad|stayer` switch added around line 160. |
| Canonical wrapper | `/mnt/gs21/scratch/willyc/soc_project/Scr_Ana/A_JMP4_refactor/JMP4_grad_leavers_stemxopt_hiringmargin.sb` | BUILT | Submitted with `EXTMARGIN_OUTCOME=stayer`; job finished. |
| §4A wrapper | `/mnt/gs21/scratch/willyc/soc_project/Scr_Ana/A_JMP4_refactor/JMP4_grad_leavers_stemxopt_hiringmargin_4a.sb` | BUILT | Submitted with `EXTMARGIN_OUTCOME=stayer`; job finished. |
| §4B wrapper | `/mnt/gs21/scratch/willyc/soc_project/Scr_Ana/A_JMP4_refactor/JMP4_grad_leavers_stemxopt_hiringmargin_4b.sb` | BUILT | Submitted with `EXTMARGIN_OUTCOME=stayer`; job finished. |
| Event-study wrapper | `/mnt/gs21/scratch/willyc/soc_project/Scr_Ana/A_JMP4_refactor/JMP4_grad_leavers_stemxopt_extmargin_eventstudy.sb` | BUILT | Submitted with `EXTMARGIN_OUTCOME=stayer`; job finished. |
| `_staying` GraphData | `_onetoe_staying` pooled outputs; `_staying` event-study outputs | BUILT / SYNCED | Eight required nonempty inputs present: canonical, §4A, §4B, event-study; reported for bachelor’s and master’s degrees. |
| Native DDD master’s table | `Output/Tables/Full_Gelbach/PDL_grad/leavers_stemxopt_ddd_native_raw_ctsXus_deg3_onetoe.csv` | BUILT | Source for master’s native-DDD results. |
| Native DDD PhD table | `Output/Tables/Full_Gelbach/PDL_grad/leavers_stemxopt_ddd_native_raw_ctsXus_deg4_onetoe.csv` | BUILT | Source for PhD native-DDD results. |
| Old abroad-first event-study master’s table | `Output/Tables/Full_Gelbach/PDL_grad/leavers_stemxopt_extmargin_eventstudy_raw_cts_deg3.csv` | BUILT | Legacy abroad-first series; do not overwrite. |
| Old abroad-first event-study PhD table | `Output/Tables/Full_Gelbach/PDL_grad/leavers_stemxopt_extmargin_eventstudy_raw_cts_deg4.csv` | BUILT | Legacy abroad-first series; do not overwrite. |
| Dynamic DD in draft | `Draft/Paper.tex` lines 664, 672, 616 | BUILT / TODO corrections | Clustering text and `US_Job` description require fixes. |
| Native DDD in draft | `Draft/Paper.tex:988` | NOT BUILT / TODO | Only commented placeholder; needs equation and prose. |
| Five-year fixed-horizon event study | Not yet coded | NOT BUILT | Suggested, not implemented. |
| True five-year employment/nonemployment outcome | Not yet coded | NOT BUILT | Requires defensible raw-panel status measure; `ue_first` is not unemployment. |
| Knitted/rendered figures | `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd` | TODO | Not knitted locally. |
| Handoff file | `/Handoffs` | TODO | This document should be stored there. |

## Proposed Changes

### Already implemented

- `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd`
  - §3C now graphs U.S. retention via Kaplan–Meier \(S(h)\), the share still U.S.-resident.
  - Removed the inverse exit-hazard panel.
  - New output names: `S3c_usRetention_stemxopt_*`.
  - Intensive-margin path reads only O*NET over-education inputs (`*_onetoe.rds`).
  - Filters ACS-based `over_edu` from legacy untagged inputs; retains `onet_over_edu`.
  - Removed ACS/O*NET source legend and duplicate-series aesthetics.
  - Extensive-margin path reads `_onetoe_staying` pooled outputs and `_staying` event-study outputs.
  - Labels the new extensive-margin figure as `P(first job in the U.S.)`, not retention.
- `Codes/hpcc_ana/A_JMP4_refactor/grad/JMP4_grad_t_leavers_stemxopt_hiringmargin.R`
  - Added `EXTMARGIN_OUTCOME=abroad|stayer`.
  - For `stayer`, estimates \(1 - y_{\text{abroad}}\) on identical rows/specification.
  - Writes only new `_staying` extensive-margin outputs; preserves abroad-first outputs.
- `Codes/hpcc_ana/A_JMP4_refactor/grad/JMP4_grad_t_leavers_stemxopt_extmargin_eventstudy.R`
  - Same `EXTMARGIN_OUTCOME=abroad|stayer` switch and `_staying` output isolation.
- `CHANGELOG.md`
  - Documented Rmd §3C retention reframing, intensive-margin ACSOE cleanup, and `_staying` extensive-margin path.

### Pending / recommended

- Knit `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd` and inspect:
  - §3C U.S.-retention Kaplan–Meier figure with `S3c_usRetention_stemxopt_*` output names.
  - Intensive-margin facet with only O*NET over-education; no ACS/O*NET legend.
  - Extensive-margin figure labeled `P(first job in the U.S.)` from `_staying` inputs.
- Correct `Draft/Paper.tex` §4.3:
  - Line 672 currently says cluster by STEM×graduation-year; replace with major-level CRV1 used in event-study code.
  - Distinguish the unpooled-pre-trend specification from the pooled-2003–06 reference used in the main dynamic graph.
  - Line 616 mistakenly says `US_Job` indicates a job outside the U.S.; correct it.
  - At line 988, replace the commented placeholder with an equation and prose for \(\mathrm{DD}_{CtS}-\mathrm{DD}_{native}\), still-U.S.-resident conditioning, native-spillover caveat, and cohort-cluster/WCR inference.
- If approved, build a fixed-horizon event study at \(h=5\):
  - Cohort-specific STEM DD, pooled 2003–06 reference, drop 2007, retain only cohorts with a fully observed five-year window.
  - Two-part figure: five-year U.S.-retention/outcome-observation rate; five-year job-quality DDs for wage, O*NET OE, mismatch, PageRank, firm value.
  - Label job-quality estimates conditional on being observed in the U.S. at \(h=5\); show retention margin alongside.
  - Follow with CtS-versus-native DDD analogue at \(h=5\).
- Do not use `ue_first` as unemployment. It means the actual first job has an SOC outside the keep-list, often an employed academic teaching/research job.

## What to Watch For

- **`stayer` vs `_staying`:** the environment switch is `EXTMARGIN_OUTCOME=stayer`, but the filename tag is `_staying`. Do not “fix” one to match the other without checking the readers and wrappers.
- **Mechanical complement vs retention:** \(1 - y_{\text{abroad}}\) is first-job U.S. placement, not “still U.S.-resident at \(h\).” The `_staying` tag does not make it the same as the Kaplan–Meier retention outcome.
- **Legacy abroad-first outputs:** keep unchanged. The `_staying` path must not overwrite `abroad_first` files.
- **`stayer10`:** not the requested mechanical complement. It is a horizon-10 outcome and mechanically falls for recent post cohorts with less than 10 years of follow-up; the producer warns against using its raw share for a DD.
- **ACSOE cleanup scope:** remove only ACS-based `over_edu`; do **not** remove the entire “Over-education” facet. Keep `onet_over_edu`.
- **Clustering:** event-study code uses major-level CRV1. The old draft language says cluster by STEM×graduation-year, which is invalid/degenerate. Native DDD requires cohort-cluster/WCR inference.
- **Pre-trends:** master’s dynamic DD pre-trend test does not reject (p=.189); PhD rejects decisively (p=.007). Do not make a causal PhD claim from the dynamic DD.
- **Timing confound:** master’s dynamic placement estimates grow from about +2 pp in 2008 to +5–7 pp in 2015–19, raising concern that the pattern partly reflects the 2016 rule change or later STEM-specific changes rather than the 2008 extension alone.
- **Selection:** master’s STEM-post graduates are +4.0 pp more likely to have an off-keep-list first job and +10.1 pp more likely never to enter the keep-list. Conditional stayer-quality results are not unconditional treatment effects.
- **`ue_first`:** not economic unemployment. A true five-year employment/nonemployment outcome needs a defensible status measure from the raw panel.
- **Degree labels:** the Rmd input check reported bachelor’s and master’s degrees; the results discussion uses master’s and PhD. Verify `deg3`/`deg4` coding and all figure labels before drafting.
- **Native DDD endpoints:** results are horizon-specific endpoints under the at-risk floor, not a joint test over all horizons. Native comparison removes common STEM shocks but not CtS-specific compositional changes or native spillovers.
- **Native DDD strength:** master’s DDD reinforces wage, firm-quality, PageRank, and mismatch results; O*NET OE is +1.5 pp, p=.086, not conventional. For PhDs, only firm value clearly survives; PageRank is borderline p=.070.
- **Employer switching:** no credible effect after rebuilding full-panel employer IDs; the earlier result is retracted.
- **Cap-exempt mechanism:** not ready to cite; destination-share evidence is inconclusive.
- **No extensive-margin DDD:** natives’ abroad-first rate is essentially zero, so there is no corresponding extensive-margin native DDD.
- **Do not execute pipeline scripts in sandbox.** Follow AGENTS: local `.py` uses `/Users/Willy/postgres_csv_env/bin/python3`; HPCC outputs go to `/mnt/gs21/scratch/willyc/soc_project/Output` and created data to `/mnt/gs21/scratch/willyc/soc_project/DataOut`.

## Execution Sequence

1. Verify every synchronized `_staying` artifact matches the filenames requested by `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd`; confirm nonempty canonical, §4A, §4B, and event-study inputs.
2. Knit/render `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd`; inspect §3C U.S.-retention, O*NET-only intensive margin, and `P(first job in the U.S.)` extensive margin.
3. Confirm no ACSOE series, ACS/O*NET legend, or duplicate-source aesthetics remain in the intensive-margin figures.
4. Correct `Draft/Paper.tex` §4.3: line 672 clustering, line 616 `US_Job` description, and the unpooled-vs-pooled reference distinction at lines 664–672.
5. Replace the native-DDD placeholder at `Draft/Paper.tex:988` with equation, prose, conditioning, spillover caveat, and inference language.
6. Decide whether to build the proposed \(h=5\) fixed-horizon event study; if approved, code it, run on HPCC, and store outputs under the established `DataOut`/`Output` paths.
7. If true unemployment is pursued, first construct a defensible five-year employment/nonemployment status variable; do not repurpose `ue_first`.
8. Update `CHANGELOG.md` for any `Draft/Paper.tex` or new-code changes, and store this handoff in `/Handoffs`.

**Next action:** Knit/render `Output/Graphs/Leavers/leaver_plots_stemxopt.Rmd` and inspect the three edited figures—§3C U.S.-retention, O*NET-only intensive margin, and `_staying` extensive-margin `P(first job in the U.S.)`—before updating `Draft/Paper.tex`.