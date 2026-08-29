# solomonR: Analyze Solomon Four-Group Designs

Implements (1) the classic Solomon analysis path (2x2 posttest ANOVA to
check pretest sensitization, then ANCOVA on pretested cells and a
posttest-only t-test), plus optional Stouffer Z "Test I" (Braver &
Braver, 1988); and (2) a unified GLM that estimates treatment, pretest
indicator, and their interaction, using ANCOVA logic for the pretested
cells and robust or permutation-based inference.

## Details

Key references: Huck & Sandler (1973) for the classic workflow; Braver &
Braver (1988) for Stouffer Z; and the Sawilowsky & Markman critiques and
replies (1988–1990) for cautions and selection bias issues. Teaching
overviews are included.
