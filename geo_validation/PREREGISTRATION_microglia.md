# Pre-registration: is the HA-receptor/DAMP signal microglia-intrinsic? (GSE96055)

Written before any expression value of GSE96055 was analysed (only sample labels and file layout seen).

Data: FACS-sorted CX3CR1-eGFP microglia, mouse SCI; non-injured (MGNI, n=3) vs injured at 72 h, 1 wk, 2 wk
(full transection MGFT + hemisection MGHS pooled as "injured"; n = 6, 5, 6). Normalized counts -> log2(x+1);
genes with mean normalized count < 5 removed.

Hypothesis context: the whole-tissue SCI/TBI data (cohorts GSE45006, GSE47681, GSE45997) show the HA
receptor/DAMP arm (Cd44, Tlr2, Tlr4, Cd14, Hmmr) going up. A critic can say this only reflects immune-cell
infiltration. If the signal is real cell-state change, the arm should rise in purified microglia themselves.

Predictions (same competitive permutation machinery as run_geo_validation.py, expression-matched null):
* M1  HA receptor/DAMP arm up in injured vs non-injured microglia at 72 h (primary), 1 wk, 2 wk. Primary = 72 h.
* M2  Full Hyaluronan axis up at 72 h.
* M3  (symmetry check) HA synthesis (Has1/2/3) and degradation (Hyal*, Cemip, Tmem2) are reported without a directional prediction.
* M4  Specificity: effect of the receptor arm is compared with the other 5 ECM axes; no claim of uniqueness is made.
Secondary: per-injury-type contrasts (full vs hemi) and per-gene effects. A negative M1 will be reported as such.
