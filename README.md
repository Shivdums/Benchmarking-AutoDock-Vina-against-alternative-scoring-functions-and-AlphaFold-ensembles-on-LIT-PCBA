# Benchmarking-AutoDock-Vina-against-alternative-scoring-functions-and-AlphaFold-ensembles-on-LIT-PCBA
A virtual-screening benchmark across five protein targets, run end to end on Google Colab. The question underneath all of it: **when AutoDock Vina fails to find the actives, what is actually to blame — the scoring function, the receptor structure, or the ligands themselves?

Everything here was produced in a single week on a Colab T4. All code is in `docking.ipynb`,
all numbers in `RESULTS/`, all figures in `FIGURES/`.
 
---
 
## Contents
 
1. [What this project is](#what-this-project-is)
2. [Assumptions and substitutions](#assumptions-and-substitutions)
3. [The five targets and why they were chosen](#the-five-targets-and-why-they-were-chosen)
4. [How to read this repository](#how-to-read-this-repository)
5. [Pipeline, step by step](#pipeline-step-by-step)
6. [Metrics](#metrics)
7. [Results](#results)
8. [Limitations](#limitations)
9. [What I would run next](#what-i-would-run-next)
10. [Biological interpretation](#biological-interpretation)
11. [Reproducing this](#reproducing-this)
---
 
## What this project is
 
The brief had three parts:
 
1. Pick ~5 proteins from a docking benchmark database and try to reproduce or beat the
   enrichment those systems report.
2. Where Vina fails significantly relative to DOCK 3.7, identify and explain the causes.
   *Bonus:* flag ligands as "DOCK does better" or "Vina does better" **a priori**, without docking.
3. Re-run everything on a predicted structural ensemble instead of the crystal.
   *Bonus:* classify the false negatives as "failed because of the ligand/software" versus
   "failed because of the protein structure", then use that to pick a smarter top-K than
   Vina's top-K, validated by leave-one-target-out.
Only the top pose per ligand is used anywhere, as specified.
 
---
 
## Assumptions and substitutions
 
These are the decisions that shape everything downstream. They are stated up front because
they are the first thing anyone reading the numbers needs to know.
 
**The LSD server was unreachable for the entire project period.**
`lsd.docking.org` timed out from three separate networks across several days. DNS resolved
correctly (169.230.75.4) and `zinc20.docking.org` returned a 302, so the outage was server-side
and partial, not a local network problem. The consequence is significant: **there are no
per-ligand DOCK 3.7 scores in this project.** The head-to-head comparison the brief asks for
in Part 2 could not be made.
 
**LIT-PCBA was used instead.**
Obtained from the HuggingFace mirror `THU-ATOM/DrugCLIP_data` (`LIT-PCBA.zip`) after the
official distribution links returned 404. LIT-PCBA is a reasonable substitute for one specific
reason: its actives and inactives are both experimentally confirmed from PubChem BioAssay
dose-response data, not computationally generated decoys, and the sets are AVE-unbiased. So
enrichment here is measured against real tested non-binders rather than property-matched
guesses. It is not a substitute for DOCK 3.7 scores, only for the benchmark systems.
 
**Part 2 is a rescoring, not a re-docking.**
AD4 and Vinardo evaluated the poses Vina had already produced and chosen. No new search was
run. This matters, and the limitations section explains why.
 
**Everything ran in the cloud.** Local machine is a Windows laptop; all compute is Colab.
 
---
 
## The five targets and why they were chosen
 
LIT-PCBA has fifteen targets. I picked five that each had a clean co-crystal ligand to define
the box, were small enough to dock ~5,000 ligands on one GPU, and — deliberately — spanned five
different *kinds* of binding site. The point was that if Vina failed, I wanted to be able to see
which kind of pocket it failed on, rather than just recording that it failed.
 
| Target | PDB | Biological role | Pocket type | What a ligand has to do there |
|---|---|---|---|---|
| **FEN1** | 5FV7 | Flap endonuclease 1. Cuts the 5′ single-stranded DNA flaps left behind during Okazaki fragment maturation and base-excision repair. | Two-metal (Mg²⁺) catalytic centre | Coordinate the metals and block the phosphodiester cut |
| **VDR** | 3A2I | Vitamin D receptor. Nuclear hormone receptor; binds calcitriol, heterodimerises with RXR, drives transcription. | Fully enclosed, hydrophobic | Mimic the secosteroid, or displace helix 12 so coactivators cannot bind |
| **MAPK1** | 4QTE | ERK2. Terminal kinase of the RAS–RAF–MEK–ERK cascade. | Classic ATP cleft | Compete with ATP, hydrogen-bond to the hinge |
| **PKM2** | 4G1N | Pyruvate kinase M2. Final step of glycolysis; the isoform tumours express and deliberately keep slow. | Allosteric, at the dimer–dimer interface | **Activate**, not inhibit — lock the enzyme in its active tetramer |
| **ALDH1A1** | 5L2O | Oxidises retinaldehyde to retinoic acid; detoxifies reactive aldehydes. | Long, solvent-exposed substrate tunnel | Plug access to the catalytic cysteine |
 
Three of these have features that stock Vina handles badly on paper: FEN1 needs metal
coordination, which Vina has no electrostatic term for; PKM2's site is an interface rather than
a cleft; ALDH1's tunnel is long and open to solvent. That was intentional. It is where the
failure analysis has something to work with.
 
PKM2 is also the one case where "better binder" means activator rather than inhibitor. Docking
scores a binding event. It knows nothing about the direction of the pharmacology.
 
---
 
## How to read this repository
 
```
.
├── README.md                       ← you are here
├── ASSUMPTIONS.md                  ← the substitutions above, standalone
├── notebook.ipynb                  ← the full pipeline, top to bottom
│
├── RESULTS/
│   ├── part1_enrichment.csv        ← headline Part 1 table, 5 targets
│   ├── part2_scoring_comparison.csv← vina / ad4 / vinardo × 5 targets
│   ├── part2_top10_hits.csv        ← actives found in each function's top 10
│   ├── part3_crystal_vs_ensemble.csv ← crystal vs AF2 ensemble, 4 targets
│   ├── part3a_false_negatives.csv  ← per-ligand rescue labels + descriptors
│   ├── part3b_loocv.csv            ← leave-one-target-out re-ranking
│   ├── key_results.json            ← every headline number in one file
│   ├── top10_all_targets.csv       ← top 10 by Vina for each target
│   │
│   ├── {TARGET}_vina_scores.csv    ← mol_id, label, vina        (crystal docking)
│   ├── {TARGET}_ad4_scoring.csv    ← lig index, ad4 score       (rescoring)
│   ├── {TARGET}_vinardo.csv        ← lig index, vinardo score   (rescoring)
│   ├── {TARGET}_all_scores.csv     ← the three above, joined
│   ├── {TARGET}_part3_scores.csv   ← ensemble docking scores
│   └── {TARGET}_descriptors.csv    ← RDKit descriptors per ligand
│
└── FIGURES/
    ├── fig1_part1_enrichment.png   ← ROC curves + logAUC bars
    ├── fig2_part2_scoring.png      ← scoring functions compared
    ├── fig3_size_bias.png          ← the VDR size-bias diagnosis
    ├── fig4_part3_ensemble.png     ← crystal vs ensemble + flexibility
    ├── fig5_part3a_classifier.png  ← rescue prediction + ablation
    └── fig6_part3b_loocv.png       ← LOOCV re-ranking
```
 
### How the files connect
 
The index columns are the thing to understand, because two different conventions are in play.
 
- `{TARGET}_vina_scores.csv` is keyed by **`mol_id`**, the PubChem CID from LIT-PCBA, and carries
  the ground-truth `label` (1 = active, 0 = inactive).
- `{TARGET}_ad4_scoring.csv` and `{TARGET}_vinardo.csv` are keyed by **`lig`**, a 1-based integer
  that is the ligand's position in the SDF handed to Uni-Dock. It is *not* the PubChem ID.
- `{TARGET}_all_scores.csv` is the join of those three, and is what every Part 2 analysis reads.
The join is positional, so it is only valid if the ligand order never changed. During the gnina
conversion step two molecules failed (one in PKM2, one in ALDH1) and silently dropped out. That
would have shifted every subsequent index by one and invalidated every label downstream. The
notebook contains a `convert_verified()` function written specifically to catch this: it
converts in small batches, identifies which individual file failed, and preserves the index
mapping. Worth knowing about if you re-run anything.
 
`mol_id` must be read as a **string**, not a float — `pd.read_csv(..., dtype={'mol_id': str})`.
Otherwise pandas turns `842469` into `842469.0` and the joins fail.
 
Flow of the whole thing:
 
```
LIT-PCBA ──► prep_target() ──► ligands.sdf + receptor.pdbqt + box.npy
                                        │
                                        ▼
                              Uni-Dock (crystal)
                                        │
                        ┌───────────────┼────────────────┐
                        ▼               ▼                ▼
              {T}_vina_scores.csv   poses ──► gnina   descriptors
                        │           --score_only          │
                        │                 │               │
                        ▼                 ▼               │
                 part1_enrichment   {T}_ad4 / {T}_vinardo │
                                          │               │
                                          ▼               │
                                  part2_comparison        │
                                                          │
ColabFold ──► 8 conformers ──► Uni-Dock (ensemble)        │
                                     │                    │
                                     ▼                    │
                           {T}_part3_scores.csv           │
                                     │                    │
                     ┌───────────────┴──────────┐         │
                     ▼                          ▼         ▼
            part3_crystal_vs_ensemble    part3a_false_negatives
                                                  │
                                                  ▼
                                           part3b_loocv
```
 
---
 
## Pipeline, step by step
 
### 1. Target preparation
 
A single reusable function, `prep_target()`, does all five. For each target it:
 
- **Defines the box** from the co-crystal ligand: takes the bound ligand's atomic coordinates,
  uses the centroid as the box centre and the bounding-box dimensions plus 8 Å padding as the
  size. This is the equivalent of eyeballing the grid in AutoDockTools, done arithmetically so
  it is reproducible.
- **Downloads the receptor** from RCSB, keeps the specified chains, strips waters and
  crystallographic additives, but **explicitly retains catalytic metals** (Mg, Zn, Mn, Fe, Co,
  Ni, Cu, Ca). This matters for FEN1, whose two Mg²⁺ ions are the catalytic machinery. The
  supplied benchmark receptor did not include them.
- **Adds hydrogens at pH 7.4** with Open Babel (`-h -p 7.4`), then converts to PDBQT with `-xr`.
- **Strips spurious metal-bonded hydrogens.** Open Babel treats Mg²⁺ as needing hydrogens and
  adds them. It should not. These are removed before conversion.
- **Builds the ligands** from SMILES: RDKit ETKDGv3 embedding with a fixed seed (42), parallel
  across cores, then Open Babel protonation at pH 7.4 and conversion to individual PDBQT files.
- **Subsamples** inactives where the full set is impractical, with a fixed seed.
PKM2 uses chains **A and B**, because its allosteric pocket sits at the interface between them.
With a single chain, the site does not exist.
 
### 2. Binding-site validation
 
Before trusting any of the boxes, each one was checked independently in **BIOVIA Discovery
Studio** (Define and Edit Binding Site → From Receptor Cavities), which finds cavities from
receptor geometry alone, with no knowledge of the co-crystal ligand.
 
| Target | Matching DS cavity | Volume (Å³) | Distance from my box centre |
|---|---|---|---|
| MAPK1 | Site 1 | 241.4 | 2.0 Å |
| PKM2 | Site 3 | 922.3 | 2.5 Å |
| VDR | Site 1 | 457.3 | 3.9 Å |
| FEN1 | *(no match)* | 46.5 | 11.0 Å |
| ALDH1 | *(no match)* | 642 | 17.6 Å |
 
Three of five agree closely. FEN1 and ALDH1 do not, and the honest reading is that geometric
cavity detection and co-crystal-ligand position disagree for those two. FEN1's catalytic centre
is a shallow, metal-defined surface rather than an enclosed pocket, so a volume-based detector
finds something small and elsewhere. ALDH1's substrate tunnel is long and open, so the detector
picks a different part of it. I kept the ligand-derived boxes, because the co-crystal ligand is
direct experimental evidence of where ligands actually bind and cavity detection is an
inference. But the disagreement is real and is recorded here rather than quietly dropped.
 
I also ran fpocket as a cross-check. It is not in apt, the conda install aborted on a libmamba
error, and a from-source build compiled but crashed with a qhull error at runtime. Abandoned.
 
### 3. Docking (Part 1)
 
Stock AutoDock Vina was benchmarked first and measured at **50.3 seconds per ligand** on
Colab's CPU. At ~5,400 ligands per target that is roughly 75 hours for one target. Not viable.
 
Switched to **Uni-Dock**, which is the Vina scoring function reimplemented in CUDA. Same
scoring function, same search algorithm, GPU-parallel. Runtime dropped to roughly 25 minutes
per target on a T4. Installed via micromamba from conda-forge, because the GitHub release
binary 404s and silently writes an HTML error page to disk — hence the file-size assertion
(`assert sz > 10_000_000`) on every binary download in the notebook.
 
```bash
unidock --receptor {rec}.pdbqt --ligand_index {ligs}.txt \
        --center_x {cx} --center_y {cy} --center_z {cz} \
        --size_x {sx} --size_y {sy} --size_z {sz} \
        --search_mode balance --num_modes 1 --dir {outdir}
```
 
`--num_modes 1` enforces the brief's top-pose-only requirement at the docking level.
 
### 4. Rescoring (Part 2)
 
The saved poses were rescored with **gnina 1.0.1** using two alternative functions:
 
```bash
gnina --score_only -r {rec}.pdbqt -l {poses}.sdf \
      --scoring ad4_scoring --cnn_scoring none
```
 
`--cnn_scoring none` disables gnina's neural network, so this is purely the classical function.
 
Three functions, one pose set:
 
| | Vina | AutoDock4 | Vinardo |
|---|---|---|---|
| Derivation | Empirical fit to measured affinities | Physics-based, AMBER-derived | Vina's form, reparameterised |
| Steric | Gaussian + repulsion | Lennard-Jones 12-6 | Reweighted Gaussian |
| H-bond | Directional distance term | Directional 12-10 | Steepened directional |
| Hydrophobic | Pairwise contact | Implicit in desolvation | Steepened contact |
| **Electrostatics** | Not modelled | **Explicit Coulomb** | Not modelled |
| **Desolvation** | Not modelled | **Explicit** | Not modelled |
 
gnina rejects multi-MODEL PDBQT input, so poses were converted to multi-molecule SDF in chunks
of 400. This is where the two molecules were lost; see the index note above.
 
### 5. Ensemble generation and docking (Part 3)
 
**ColabFold with MSA subsampling**, which is the "rMSA AF2" mode the brief points at:
 
```
--max-msa 16:32 --use-dropout --num-seeds 8
```
 
Cutting the alignment down makes AlphaFold less certain about the answer. Instead of collapsing
onto a single high-confidence structure it samples several plausible backbone states. Eight
models per target.
 
AlphaFold writes output in its own coordinate frame, so each model was superposed onto the
crystal in PyMOL and the original box reused unchanged.
 
**A note on measuring the diversity:** I first used `cmd.align`, which reported 0.5–0.9 Å and
suggested the ensemble was useless. That was wrong — `align` does outlier rejection and
iterative refinement, so it reports how well the structures *can* be made to agree, not how
different they are. Switching to `cmd.rms_cur` with no rejection showed real backbone diversity
of 2.3–5.9 Å. Worth flagging because the first number would have ended Part 3 prematurely.
 
Pocket side-chain RMSD across each ensemble:
 
| Target | Pocket RMSD (Å) |
|---|---|
| MAPK1 | 0.51 |
| FEN1 | 0.75 |
| VDR | 1.10 |
| ALDH1 | 1.15 |
| PKM2 | 3.13 |
 
Every ligand was docked against all eight conformers (`--search_mode fast`) and its **best score
across the ensemble** taken as its score. The library was reduced to ~1,300–2,000 ligands per
target to fit the compute budget.
 
**PKM2 was excluded from Part 3.** AlphaFold returned monomers — 3,964 atoms against 7,943 in
the crystal dimer. PKM2's pocket is at the dimer interface, so in the predicted models it simply
does not exist. Part 3 is four targets.
 
### 6. Machine learning (Parts 3a and 3b)
 
**3a — classifying the false negatives.** Every active that Vina ranked badly on the crystal was
labelled *rescued* (the ensemble moved it into the top 5%) or *not rescued*. A random forest was
trained to predict that label under stratified 5-fold CV, with features ablated into three
groups: crystal rank, docking scores, and RDKit ligand descriptors (heavy atom count, logP,
TPSA, rotatable bonds, HBD/HBA, Fsp3, ring count, formal charge).
 
**3b — leave-one-target-out re-ranking.** Train on four targets, predict the fifth, rotate
through all five. The model re-ranks each target's crystal-docked list and its top K is compared
against Vina's top K at the same K.
 
---
 
## Metrics
 
```python
def met(df, col):
    d = df.sort_values(col); y, s = d.label.values, -d[col].values
    fpr, tpr, _ = roc_curve(y, s); m = fpr >= 1e-3
    x  = np.concatenate(([1e-3], fpr[m]))
    yv = np.concatenate(([np.interp(1e-3, fpr, tpr)], tpr[m]))
    la = (np.trapezoid(yv, np.log10(x))/3 - 0.14463)*100
    k = max(1, int(.01*len(d)))
    return round(roc_auc_score(y,s),3), round(la,2), round((d.head(k).label.sum()/k)/y.mean(),2)
```
 
- **ROC-AUC** — standard. 0.5 is random. Weights the whole list equally, which is not what
  screening cares about.
- **Adjusted logAUC** — the Shoichet-lab metric. Log-weighted from λ = 0.001, so the top 0.1%,
  1% and 10% each contribute roughly equally. The random baseline (0.14463) is subtracted, so
  **0 means random and negative means worse than random.** This is the primary metric throughout.
- **EF1% / EF5%** — enrichment factor. Of the top 1% of the ranked list, how many times more
  actives than you would get picking at random. 1.0 is random.
EF1% and logAUC can disagree, and here they sometimes do. They are measuring different slices.
EF1% looks only at the first ~55 molecules and ignores everything below. logAUC summarises the
whole early range. A function that puts a few extra actives in the first handful but is mediocre
across the next few hundred will win on EF1% and lose on logAUC. Which one matters depends
entirely on how many compounds you can afford to buy.
 
*(`np.trapz` was removed in NumPy 2.x; `np.trapezoid` is the replacement.)*
 
---
 
## Results
 
### Part 1 — Vina on crystal structures
 
`RESULTS/part1_enrichment.csv` · `FIGURES/fig1_part1_enrichment.png`
 
| Target | N | Actives | AUC | logAUC | EF1% | EF5% |
|---|---|---|---|---|---|---|
| PKM2 | 5,543 | 546 | 0.633 | **+7.20** | 2.77 | 2.24 |
| MAPK1 | 5,298 | 308 | 0.658 | **+6.92** | 0.99 | 2.02 |
| ALDH1 | 5,992 | 999 | 0.574 | +4.13 | 2.03 | 1.62 |
| FEN1 | 5,368 | 369 | 0.553 | +2.35 | 0.82 | 1.36 |
| VDR | 5,882 | 882 | 0.427 | **−3.08** | 0.69 | 0.64 |
 
Four of five beat random. None beat it by much. VDR is below random, meaning it actively sorts
actives *below* inactives.
 
This is the realistic picture for stock Vina on an experimentally-grounded, unbiased benchmark.
Published enrichments on decoy-based sets are usually higher, and part of that gap is that
property-matched decoys are easier to reject than real tested non-binders.
 
### Part 2 — Alternative scoring functions
 
`RESULTS/part2_scoring_comparison.csv` · `FIGURES/fig2_part2_scoring.png`
 
**logAUC**
 
| Target | Vina | AD4 | Vinardo |
|---|---|---|---|
| PKM2 | **7.20** | 3.00 | 4.07 |
| MAPK1 | **6.92** | 4.13 | 5.56 |
| ALDH1 | **4.13** | 3.55 | 2.61 |
| FEN1 | **2.35** | −2.94 | 1.12 |
| VDR | −3.08 | −5.71 | **−2.48** |
 
**EF1%**
 
| Target | Vina | AD4 | Vinardo |
|---|---|---|---|
| FEN1 | 0.82 | 0.55 | **1.92** |
| VDR | **0.69** | 0.11 | 0.34 |
| MAPK1 | 0.99 | 1.32 | **1.65** |
| PKM2 | **2.77** | 0.55 | 1.11 |
| ALDH1 | **2.03** | 1.12 | 1.73 |
 
Vina wins on logAUC at every target. Vinardo wins on EF1% at three of five. AD4 — the most
physically complete of the three, the only one with electrostatics and desolvation — is the
worst performer everywhere, and on FEN1 it inverts the ranking entirely (+2.35 → −2.94).
 
Two readings, and both are probably partly true. The first is methodological: AD4 is scoring
poses that Vina optimised for Vina's function, so it is being penalised for geometry it would
never have chosen. The second is more interesting: **Vina's omissions act as regularisation.**
Having no electrostatic term means having no *wrong* electrostatic term. A charged ligand near a
charged residue gets a number Vina simply does not compute; AD4 computes it, with no explicit
water to screen the charge, and gets it wrong.
 
The clean test is a true re-dock. See limitations.
 
### Why VDR fails
 
`FIGURES/fig3_size_bias.png`
 
| | Mean heavy atoms |
|---|---|
| VDR actives | 23.0 |
| VDR inactives | 25.2 |
| **Vina's top 50** | **32.3** |
 
Correlation between Vina score and heavy atom count: **r = −0.321**.
 
Vina's score is a sum over atom pairs, and because the dominant terms are attractive, adding
atoms generally improves the score regardless of fit. This is a known property of unnormalised
docking scores.
 
Usually it is survivable, because in most libraries the actives *are* the larger molecules. On
VDR they are not — the actives are the smaller population. So the bias points away from the
right answer and the ranking inverts. That is why VDR is at −3.08 rather than merely low.
 
This also explains why rescoring did not help: AD4 and Vinardo are sums over atom pairs too, and
carry the same bias. All three are negative on VDR.
 
The standard correction is ligand efficiency (score ÷ heavy atom count). **I did not test it,
and I should have.** It is a two-minute calculation on data already in this repo.
 
### Part 3 — Predicted ensembles
 
`RESULTS/part3_crystal_vs_ensemble.csv` · `FIGURES/fig4_part3_ensemble.png`
 
| Target | N | Actives | logAUC xtal → ens | EF1% xtal → ens | Pocket RMSD |
|---|---|---|---|---|---|
| FEN1 | 1,312 | 355 | 2.45 → **9.52** | 0.85 → **2.27** | 0.75 Å |
| VDR | 1,878 | 881 | −3.04 → **0.02** | 0.59 → 0.83 | 1.10 Å |
| MAPK1 | 1,307 | 308 | 7.20 → 6.71 | 0.98 → 1.31 | 0.51 Å |
| ALDH1 | 1,992 | 998 | 5.17 → 4.76 | 1.58 → 1.58 | 1.15 Å |
 
**The ensemble helped exactly the two targets that were failing and slightly hurt the two that
were working.** FEN1 nearly quadrupled. VDR crossed from below-random to neutral. MAPK1 and
ALDH1 each lost about 0.5, which is small enough that it may be noise.
 
**Side-chain flexibility does not predict benefit.** ALDH1 moves most (1.15 Å) and gains
nothing. MAPK1 moves least (0.51 Å) and loses slightly. Four points, so this is an observation,
not a trend — but it is a useful negative result, because "measure the pocket RMSD and decide
whether to build an ensemble" is an intuitive heuristic and it does not work here.
 
**The FEN1 result is confounded.** AlphaFold predicts no metals, so the ensemble receptors have
no Mg²⁺. The gain could be genuine backbone flexibility, or simply the removal of two ions Vina
was scoring badly. These cannot be separated without transplanting the metals into the AF2
models and re-running.
 
### Part 3a — Predicting which misses get rescued
 
`RESULTS/part3a_false_negatives.csv` · `FIGURES/fig5_part3a_classifier.png`
 
Of 2,250 actives that Vina ranked badly on the crystal, **163 (7.2%) were rescued** by the
ensemble.
 
Feature ablation, ROC-AUC:
 
| Features | ROC-AUC |
|---|---|
| Crystal rank only | 0.713 |
| + docking scores | 0.772 |
| **Ligand chemistry only** | **0.848** |
| Everything | 0.886 |
 
Precision on the positive class 0.41, recall 0.37.
 
Rescue rate by target: FEN1 15%, VDR 8%, MAPK1 5%, ALDH1 5% — which tracks the Part 3
enrichment gains, so the two analyses agree with each other.
 
**The ablation is the result.** Ligand chemistry alone (0.848) beats everything the docking run
produced (0.772). Whether a different receptor conformer will rescue a ligand is mostly a
property of the molecule, not of how it scored. Small, rigid, low-polarity ligands are the ones
a slightly different pocket shape accommodates; large flexible ones fail for reasons a new
backbone does not touch.
 
Two caveats that matter. Only 7% get rescued, so the classes are badly imbalanced — AUC 0.886
looks strong, but precision 0.41 means that when the model flags something, it is right about
four times in ten. And the CV pooled all four targets, so ligands from the same target appear in
both train and test folds. The model may be partly learning target identity rather than a
general chemical rule. Leave-one-target-out would settle it, and the machinery exists in 3b.
 
### Part 3b — Leave-one-target-out re-ranking
 
`RESULTS/part3b_loocv.csv` · `FIGURES/fig6_part3b_loocv.png`
 
| Target | K (1%) | EF1% Vina → ML | EF5% Vina → ML |
|---|---|---|---|
| FEN1 | 53 | 0.82 → **1.65** | 1.36 → 1.14 |
| VDR | 58 | 0.69 → **1.26** | 0.64 → **1.50** |
| MAPK1 | 52 | 0.99 → **1.98** | 2.02 → 1.69 |
| PKM2 | 55 | 2.77 → 1.85 | 2.24 → 1.17 |
| ALDH1 | 59 | 2.03 → 1.12 | 1.62 → 1.18 |
 
Three of five improved at 1%; one of five at 5%. Mean EF1% went from 1.46 to 1.57 — barely
anything. Mean EF5% went *down*, 1.58 to 1.34.
 
But the split is clean and it is the whole result: **every target where Vina was at or below
random improved, and both targets where Vina was working got worse.** The model pulls everything
toward a middling EF of roughly 1.1–2.0.
 
That is precisely what a model trained on *other* targets should do. It has learned a
target-independent prior about what an enriched molecule looks like. A general prior reliably
beats a broken target-specific signal and reliably loses to a working one.
 
The honest framing: this is not a better screen, it is a more predictable one. After re-ranking
no target falls below random, where three of five did before. On a target where you have no
prior idea whether docking will work, capping the downside may be worth more than the upside you
give up. But it does not meet the "show that this is consistent" bar the brief set, and I would
not deploy it.
 
---
 
## Limitations
 
Ordered by how much they would change the conclusions.
 
1. **No DOCK 3.7 comparison.** The central comparison of the brief could not be made. LIT-PCBA
   substitutes for the benchmark systems, not for the scores.
2. **Part 2 is rescoring, not re-docking.** AD4 and Vinardo never searched. This biases the
   comparison against them and is the most likely reason AD4 looks as weak as it does.
3. **Single seed throughout, no replicates.** FEN1 docked with and without Mg²⁺ gave similar
   enrichment (logAUC 2.35 vs 1.53) but only **r = 0.47** on per-ligand scores. That is a lot of
   run-to-run variance, and I never quantified it with a same-receptor different-seed control.
   Several of the smaller differences in this repo may not survive one.
4. **The FEN1 ensemble gain is confounded** by the absence of metals in the AF2 models.
5. **The ensemble gets eight attempts; the crystal gets one.** Taking the best score across
   conformers mixes structural benefit with extra sampling. The control is eight crystal runs
   with different seeds. Not done.
6. **Part 3a's validation pools targets.** Possible target-identity leakage. LOOCV would fix it.
7. **PKM2 dropped from Part 3.** Part 3 is four targets, not five.
8. **Five targets.** Every cross-target pattern here rests on five points. Observations, not
   trends.
9. **Ligand efficiency never tested**, despite the size bias being the cleanest finding here.
10. **One protonation state per ligand**, pH 7.4, no tautomer enumeration, one ETKDG conformer
    as the docking input.
11. **Part 2's bonus** — an a priori classifier for "which scoring function will win on this
    ligand" — was not built. With no DOCK 3.7 labels the original version was impossible, and
    the Vina-vs-Vinardo version ran out of time.
12. **FEN1 and ALDH1 box placement disagrees with geometric cavity detection** by 11 and 18 Å.
    I trusted the co-crystal ligand; that is a defensible choice but it is a choice.
---
 
## What I would run next
 
In order of information gained per GPU-hour:
 
1. **Ligand-efficiency normalisation** on the existing scores. One afternoon, no new docking,
   directly tests the clearest finding in the project.
2. **Eight crystal runs with different seeds.** Separates real ensemble gains from sampling
   variance, and finally puts an error bar on everything.
3. **Leave-one-target-out for the 3a classifier.** Tests whether the chemistry rule generalises.
4. **True re-docking with AD4 and Vinardo.** Makes the Part 2 comparison fair.
5. **Metals transplanted into the FEN1 ensemble.** De-confounds the biggest Part 3 result.
6. **PKM2 as a dimer**, via AlphaFold-Multimer, to bring it back into Part 3.
7. **DOCK 3.7**, if LSD returns. Everything here is structured so those scores drop straight in.
---
 
## Biological interpretation
 
Running back through all five targets, the results line up with the biology more cleanly than
I expected, and in a way that says something fairly specific about what Vina is and is not
modelling.
 
**Vina does well exactly where its central assumption is true.** The scoring function is, at
bottom, a statement that affinity is dominated by shape complementarity and buried hydrophobic
surface. MAPK1 (+6.92) and PKM2 (+7.20) are the two best-performing targets, and they are the two
whose sites are enclosed and largely hydrophobic. A kinase ATP cleft is a well-defined slot; the
PKM2 allosteric pocket is a pocket between two subunits. In both, "does this molecule fit, and
does it bury a lot of greasy surface" is genuinely most of the binding story. Vina's crudeness
is not a handicap there.
 
**It does worst where the defining interaction is one it cannot see.** FEN1's active site is two
magnesium ions. A real FEN1 inhibitor binds by chelating them, that is the pharmacophore. Vina
has no electrostatic term, so it treats Mg²⁺ as a generic atom with a radius and no charge. The
single most important interaction in the site is invisible to the scoring function, and FEN1 is
the weakest of the four targets that worked (+2.35). The instructive part is what happened when I
added the missing physics: AD4 *does* have Coulomb and desolvation terms, and on FEN1 it went to
−2.94. Having an electrostatic term and applying it to a pose that was never optimised for it was
considerably worse than not having one at all. Incomplete physics, consistently applied, beat
complete physics applied to the wrong geometry.
 
**ALDH1 sits in between, and for a structural reason.** Its substrate tunnel is long, narrow and
open to solvent. Vina has no explicit water. In a deep enclosed pocket that omission is cheap,
because displaced water is roughly constant across ligands and partly absorbed into the fitted
weights. In a solvent-exposed tunnel it is not, how much ordered water a given ligand displaces
varies ligand to ligand, and that variation is exactly the part Vina cannot model. ALDH1 lands
mid-table (+4.13), which is about where you would expect a target whose binding is real but whose
energetics are partly aqueous.
 
**VDR is the most interesting failure, because it is not a structural failure at all.** The
pocket is a textbook enclosed hydrophobic cavity, precisely the geometry Vina should handle
best. And it is the one target that scores below random. The explanation turned out to be
arithmetic rather than physical: Vina's score accumulates over atom pairs, so larger molecules
score better almost regardless of fit (r = −0.321 with heavy atom count). On most libraries that
bias is harmless because actives skew large. On VDR the actives are *smaller* than the inactives
(23.0 vs 25.2 heavy atoms), because they are secosteroid analogues of a specific natural
hormone — chemically constrained to roughly calcitriol's size. So the bias points directly away
from the answer and the ranking inverts.
 
That is worth dwelling on. The failure had nothing to do with the protein, the pose, or the
quality of the structure. It came from the intersection of a known scoring artefact with an
unusual property of that target's chemical series. No better conformer fixes it — and indeed the
ensemble only lifted VDR to neutral (−3.04 → 0.02), never to good, which is exactly what you
would predict if the problem is a size bias rather than a shape mismatch.
 
**The ensemble results point the same direction.** The targets that benefited were the ones
already failing, and the amount of benefit had no relationship to how much the pocket side chains
moved. Meanwhile the classifier found that whether a missed active gets rescued is predicted
better by the molecule's own chemistry (AUC 0.848) than by anything the docking produced (0.772).
Small, rigid, low-polarity ligands are recovered by a slightly different pocket; large flexible
ones are not, because their failure was never about pocket shape in the first place.
 
**So the through-line across every part of this project is the same.** The failures tracked the
ligands more than the protein. Receptor conformation, scoring-function physics, and pose quality
each mattered less than molecular size, rigidity and polarity. For a screening campaign that is
genuinely encouraging news, because ligand descriptors are free to compute and conformational
ensembles are not. Before spending GPU time generating AlphaFold conformers, it is worth looking
at the chemical series you are screening and asking whether its failures are the kind a new
structure could plausibly fix.
 
And one caveat that applies to the whole thing: docking scores a binding event. On PKM2 the
therapeutic goal is activation, not inhibition, the ligands worth finding are the ones that lock
the enzyme in its active tetramer. Vina cannot distinguish an activator from an inhibitor from a
tight-binding inert molecule. Enrichment here measures the ability to find binders. Whether those
binders do the thing you want is a question for an assay.
 
---
 
## Reproducing this
 
Open `notebook.ipynb` in Colab, **set the runtime to GPU (T4)** — not TPU — and mount Drive.
Every cell that depends on a prior stage re-reads from Drive, so the notebook survives runtime
disconnects, which it needed to several times.
 
Install block:
 
```python
!pip -q install rdkit tqdm scikit-learn seaborn
!apt-get -qq install -y openbabel
!wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
!./bin/micromamba create -y -p /content/ud -c conda-forge unidock
```
 
Check the GPU is real before launching anything long:
 
```python
import subprocess
gpu = subprocess.run(['nvidia-smi','--query-gpu=name','--format=csv,noheader'],
                     capture_output=True, text=True).stdout.strip()
print('GPU:', gpu or '*** NONE — switch runtime to T4 ***')
```
 
Rough runtimes on a T4: preparation ~10 min per target, crystal docking ~25 min per target,
rescoring ~5 min per target, ColabFold ~20 min per target for eight models, ensemble docking
~40 min per target, ML ~2 min total.
 
Known traps, all of which cost me time:
 
- Colab defaults to TPU in some sessions. Check.
- The Uni-Dock GitHub release binary 404s and writes an HTML page. Assert file size on download.
- smina from SourceForge returns HTML too. gnina 1.0.1 from GitHub releases works.
- gnina rejects multi-MODEL PDBQT. Convert to multi-molecule SDF in chunks.
- An empty `out.tar.gz` will make a naive "already done, skip" check pass. Count archive
  contents, not file existence.
- `np.trapz` is gone in NumPy 2.x.
- `from rdkit import Descriptors` is wrong; it is `from rdkit.Chem import Descriptors`.
- A shell `!` command after a semicolon is a SyntaxError in Colab; use `os.system`.
---
 
*Shivam Nichit · November 2026*

