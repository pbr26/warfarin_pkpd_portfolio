# Learning Notes — Warfarin PK/PD, explained for a first-timer

These notes explain **what we are doing and why**, step by step, so you can
recreate the whole project and understand it — not just run it. The R scripts
themselves stay clean; all the teaching lives here.

Read this top to bottom once, then keep it open beside the scripts.

---

## 0. The big picture in plain English

**Pharmacokinetics (PK)** = what the *body does to the drug* (absorption,
distribution, clearance). Answers: "what is the drug concentration over time?"

**Pharmacodynamics (PD)** = what the *drug does to the body* (the effect).
Answers: "given that concentration, how big is the effect?"

**PK/PD modelling** links the two: dose → concentration → effect.

**Warfarin** is a blood thinner. It works by blocking the production of
clotting factors. We measure its effect as **PCA** (prothrombin complex
activity) — lower PCA means "more anticoagulated" (blood clots less easily).
In the clinic this is usually reported as **INR**; PCA and INR are related, but
our dataset records PCA, so that's what we model.

**Why warfarin is a great teaching case:** it has a narrow safety margin (too
little = clots, too much = bleeding), a well-understood mechanism, and a
clear *delay* between concentration and effect. That delay is the heart of the
modelling.

**"Population" model:** instead of fitting each person separately, we fit
everyone at once. We estimate (a) typical parameter values for the population
and (b) how much individuals vary around those typical values. This is called
a **nonlinear mixed-effects model (NLME)**. `nlmixr2` is the R package that
does this.

---

## 1. Tools you'll install

- **nlmixr2** — the engine that fits population PK/PD models. It pulls in
  `rxode2` (solves the differential equations) and `nlmixr2data` (the data).
- **rxode2** — simulates models forward in time (used in script 07).
- **ggplot2 / dplyr / tidyr** — plotting and data wrangling.
- **xpose** — extra goodness-of-fit diagnostic plots (optional).

`setup.R` installs anything missing and loads them. Run it once:
`source("setup.R")`.

> First install can take a while — nlmixr2 compiles C++ under the hood. You
> also need a working compiler (Rtools on Windows, Xcode command-line tools on
> Mac). If installation fails, that's almost always the missing compiler.

---

## 2. The data (`nlmixr2data::warfarin`)

We do **not** download anything. The dataset ships with the package.

It is in **long format** — one row per event (a dose or an observation).
Key columns:

| Column | Meaning |
|--------|---------|
| `id`   | subject number (1–32) |
| `time` | hours since dosing |
| `amt`  | dose amount (only on dosing rows) |
| `dv`   | the measured value ("dependent variable") |
| `dvid` | which value: `"cp"` = concentration, `"pca"` = effect |
| `evid` | event id: `0` = observation, `1` = a dose |
| `wt`, `age`, `sex` | covariates (subject characteristics) |

**Why long format matters:** the same subject has both concentration rows and
effect rows, distinguished by `dvid`. When we fit PK only, we keep the dosing
rows + `cp` rows. For PD only, dosing rows + `pca` rows. For the joint model,
we keep everything.

That's exactly what **`01_data_cleaning.R`** does: loads the data, snapshots a
raw CSV (so the project is self-contained), then builds three analysis sets
(`pk_data`, `pd_data`, `pkpd_data`) and saves them. The `stopifnot()` lines are
sanity checks — they stop the script if something is obviously wrong (e.g. a
dose with zero amount).

---

## 3. Exploratory data analysis (`02_EDA.R`)

Before modelling, always *look* at the data. We make four plots:

1. **Concentration profiles** — spaghetti plot of `cp` vs time. Expect a rise
   (absorption) then fall (clearance).
2. **Effect profiles** — `pca` vs time. Expect PCA to *drop* after dosing
   (drug suppresses clotting factors) then slowly recover.
3. **Covariate distributions** — histograms of weight and age, to see the
   range of people.
4. **Hysteresis loop** — the important one. Plot effect (`pca`) against
   concentration (`cp`) and connect the dots *in time order*.

**How to read the hysteresis loop:** if effect depended *instantly* on
concentration, this plot would be a single line. Instead it forms a **loop**,
because when concentration is falling the effect is still near its peak — the
effect *lags* the concentration. This delay tells us the PD model must be an
**indirect-response (turnover) model**, not a simple "effect = f(concentration
right now)" model. This is the single most important modelling decision in the
project, and we *derived it from a picture*.

---

## 4. The PK model (`03_PK_Model.R`)

We describe concentration over time with a **compartment model** — imagine the
body as one or two connected "tanks".

**One-compartment, first-order oral absorption** (our base model):

- `depot` = the gut (where the pill dissolves).
- `center` = the bloodstream.
- Drug moves gut → blood at rate `ka` (absorption).
- Drug leaves blood at rate `cl/v` (clearance ÷ volume).

In the code this is written as two differential equations:

```
d/dt(depot)  = -ka*depot                    # gut empties into blood
d/dt(center) =  ka*depot - (cl/v)*center     # blood fills, then clears
cp = center/v                                # concentration = amount/volume
```

**The `ini({ ... })` block** sets starting guesses for the parameters:
- `tka`, `tcl`, `tv` are the *typical* (population) values. We write them as
  `log(...)` because rates and volumes must be positive; fitting on the log
  scale guarantees that (`exp()` in the model turns them back).
- `eta.ka ~ 0.3` etc. are **between-subject variability** terms — how much
  individuals differ from the typical value. The `~` means "this is a random
  effect with this starting variance".
- `add.err` / `prop.err` are **residual error** (measurement noise):
  additive + proportional.

**Why compare 1- vs 2-compartment?** We don't assume the structure — we test
it. The two-compartment model adds a "peripheral tank" (tissue). We fit both
and compare **AIC** (Akaike Information Criterion). Lower AIC = better fit
*after penalising extra parameters*. The script automatically keeps whichever
wins and saves it to `models/pk_fit.rds`.

**`est = "saem"`** chooses the SAEM algorithm — a robust way to estimate
population parameters. (FOCEi is the classic alternative.)

---

## 5. The PD model (`04_PD_Model.R`)

This is the **indirect-response / turnover** model the hysteresis told us we
need.

**Mental model:** think of PCA as the level of water in a tub.
- A tap fills it at rate `kin` (the body producing clotting factors).
- A drain empties it at rate `kout` (natural loss).
- At baseline the tub sits at `kin/kout` (fill = drain).
- **Warfarin turns the tap down** — it inhibits production. So PCA slowly
  falls, and recovers slowly when the drug wears off. Slow = the delay we saw.

The inhibition uses an **Imax model**:

```
inh = 1 - (imax * cp) / (ic50 + cp)
```

- `imax` = the *maximum* fraction by which production can be inhibited (0–1).
  We fit it on a `logit` scale so it stays between 0 and 1.
- `ic50` = the concentration that gives *half* of the maximum inhibition.
- When `cp` is 0, `inh` = 1 (no inhibition, tap full open). As `cp` rises,
  `inh` drops toward `1 - imax`.

Then:

```
pca(0)    = kin/kout          # start at baseline
d/dt(pca) = kin*inh - kout*pca # tap (inhibited) minus drain
```

We fit this to the `pca` observations. This is a **sequential** PD fit — it
uses the concentration already recorded in the data. (In the next script we do
it *jointly*, letting the PK model supply the concentration.)

---

## 6. The joint PK/PD model + covariates (`05_PKPD_Model.R`)

Now we fit **everything at once**: the PK equations produce `cp`, and that same
`cp` drives the PD equations — all parameters estimated together. The model has
**two endpoints**, declared at the bottom:

```
cp  ~ add(add.err) + prop(prop.err)   # concentration observations
pca ~ add(pca.err)                     # effect observations
```

nlmixr2 knows which rows are which from the `dvid`/compartment mapping.

**Covariates.** People differ. A common, physiologically-motivated covariate is
**body weight on clearance and volume** (bigger people clear/hold more drug).
We add it *allometrically*, normalised to 70 kg:

```
cl = exp(tcl + eta.cl) * (wt/70)^wt_cl
v  = exp(tv  + eta.v)  * (wt/70)^wt_v
```

- `wt_cl` starts at 0.75 and `wt_v` at 1.0 — the classic allometric exponents.
  Estimating them (rather than fixing) lets the data speak.
- Adding a covariate should *reduce unexplained variability* (`eta.cl`) and
  improve the objective function if it truly matters.

Saved to `models/pkpd_fit.rds` — this is the model we diagnose and simulate.

---

## 7. Diagnostics — does the model actually fit? (`06_Model_Diagnostics.R`)

A fitted model is worthless until you check it. Four standard checks:

1. **Observed vs predicted (DV vs PRED / IPRED).** Points should scatter around
   the red identity line. `PRED` uses population parameters; `IPRED` uses each
   subject's individual estimates and should hug the line more tightly.

2. **CWRES vs time and vs prediction.** CWRES = conditional weighted residuals
   (standardised errors). They should be **centred on 0 with no trend** — the
   smooth line should be roughly flat at 0. A curve or funnel means the model
   is missing structure.

3. **Eta shrinkage.** If shrinkage is high (say >30%), the data can't really
   inform that individual random effect — interpret individual estimates with
   caution.

4. **Visual Predictive Check (VPC).** Simulate hundreds of virtual trials from
   the model and overlay the observed data. If the real data falls inside the
   simulated prediction intervals, the model reproduces reality. This is the
   most persuasive plot for a portfolio. (It's wrapped in `tryCatch` so the
   script won't crash if the VPC setup hiccups.)

All figures land in `figures/`.

---

## 8. Simulation — using the model to answer questions (`07_Simulation.R`)

The payoff: once the model is trusted, we can ask "what if?" without new
experiments. We rebuild the model in **rxode2** using the fitted *typical*
parameter values (`fixef()` pulls them out; we back-transform with `exp()` /
`plogis()` because they were fit on log / logit scales).

Two scenarios:
- **A: single 100 mg dose** — matches the study, a sanity check.
- **B: 5 mg daily for 10 days** — a maintenance regimen. The `et(...)` event
  table with `ii = 24, addl = 9` means "interval 24 h, 9 additional doses".

> **Honesty flag:** the source data is *single-dose only*. Simulating repeated
> dosing is an **extrapolation** — the model wasn't trained on steady-state
> data. We label it "illustrative". Saying this out loud is exactly the kind of
> scientific honesty reviewers look for.

We also compute a simple **time-to-target**: when does PCA first fall below a
chosen threshold. In real warfarin dosing this maps to "time to reach
therapeutic anticoagulation".

---

## 9. The report (`analysis.qmd`)

Quarto stitches the story together — objective, data, EDA, each model,
diagnostics, simulation, conclusions — pulling in the saved figures. Render with
`quarto::quarto_render("analysis.qmd")` (or the Render button in RStudio) to get
a single self-contained HTML in `reports/`. This is the thing you'd actually
show someone.

---

## 10. Order of operations (cheat sheet)

```r
source("setup.R")                       # once
source("scripts/01_data_cleaning.R")    # -> data/processed/
source("scripts/02_EDA.R")              # -> figures/02_*
source("scripts/03_PK_Model.R")         # -> models/pk_fit.rds
source("scripts/04_PD_Model.R")         # -> models/pd_fit.rds
source("scripts/05_PKPD_Model.R")       # -> models/pkpd_fit.rds
source("scripts/06_Model_Diagnostics.R")# -> figures/06_*
source("scripts/07_Simulation.R")       # -> figures/07_*
quarto::quarto_render("analysis.qmd")   # -> reports/analysis.html
```

Each script reads what the previous one saved, so run them in order.

---

## 11. Vocabulary quick-reference

- **Compartment** — a modelled "tank" the drug moves through.
- **ka / CL / V** — absorption rate / clearance / volume of distribution.
- **kin / kout** — production and loss rates in a turnover model.
- **Imax / IC50** — max inhibition and the concentration giving half of it.
- **eta (η)** — a subject's deviation from the typical value (random effect).
- **Residual error** — leftover noise between model and observation.
- **OFV / AIC** — model fit scores; lower is better (AIC penalises complexity).
- **CWRES** — conditional weighted residuals (standardised errors).
- **Shrinkage** — how little the data informs an individual estimate.
- **VPC** — visual predictive check (simulate-and-overlay validation).
- **SAEM / FOCEi** — algorithms that estimate the population parameters.

---

## 12. Common gotchas

- **Compiler missing** → nlmixr2 install fails. Install Rtools (Windows) /
  Xcode CLT (Mac) first.
- **Model won't converge / crazy estimates** → your `ini()` starting guesses
  are far off. Nudge them toward the EDA (e.g. eyeball the half-life for CL/V).
- **Fitting on the natural scale** → parameters can go negative and break the
  ODE. Always fit rates/volumes on the `log` scale (we do).
- **Confusing PCA with INR** → this dataset is PCA. Don't relabel it INR.
- **Simulating multi-dose from single-dose data** → fine to show, but label it
  as extrapolation.
- **`object 'CWRES' not found` in diagnostics** → SAEM fits don't compute CWRES
  automatically. Call `fit <- addCwres(fit)` first (script 06 now does this).
- **Sequential PD needs a `cp` column** → the long data has cp and pca on
  separate rows, so script 04 interpolates each subject's concentration onto
  the PD times before fitting.
