---
title: "Warfarin Population PK/PD Modelling"
subtitle: "A Complete, First-Principles Tutorial — From Zero to a Published Analysis"
author: "Pramod BR"
date: "2026-07-24"
lang: en
toc-title: "Contents"
---

# Preface

This book is the companion to a real, working pharmacometrics project: a
population pharmacokinetic/pharmacodynamic (PK/PD) analysis of the anticoagulant
drug **warfarin**, built in R with the `nlmixr2` package and published as a live
report on Posit Connect Cloud.

It is written for one specific reader: someone who is **new to pharmacometrics**
but serious about learning it properly — well enough to recreate this project
from scratch, explain every decision, and go on to build their own. Nothing is
assumed beyond curiosity and a willingness to think. Every term is defined the
first time it appears. Every modelling choice is justified with the scientific
reason behind it, not just the code that implements it.

The book has three intertwined threads:

1. **The science** — what pharmacokinetics and pharmacodynamics actually are,
   the biology of warfarin, and the mathematics of compartment and
   turnover models, derived from first principles.
2. **The craft** — the exact R code, script by script, that turns raw data into
   a fitted model, diagnostic plots, and simulations, with each figure from the
   project reproduced and interpreted.
3. **The profession** — how this same workflow is used inside the
   pharmaceutical industry to design doses and support drug approvals, and a
   concrete, self-guided path to learn the field independently.

Read it start to finish and you will understand not only *what* the project
does, but *why* every line exists and *where* the discipline sits in the real
world.

\newpage

# How to use this book

The chapters build on each other. Part I develops the concepts (pharmacology,
the warfarin story, the mathematics of PK, PD, and population modelling). Part
II walks through the project itself, one script at a time, with the generated
figures embedded and explained. Part III steps back to the industry context and
a learning roadmap, and closes with an extensive question-and-answer chapter, a
glossary, and the complete code listings.

If you have the project open in RStudio alongside this book, you can run each
script as you reach its chapter and watch the same outputs appear. If you only
have the book, every figure and result is reproduced here.

A note on notation. Mathematical symbols appear in maths formatting, for example
$C_p$ for plasma concentration or $\eta$ (the Greek letter *eta*) for an
individual's random deviation. R code appears in fixed-width blocks. Key terms
are **bold** on first use and collected in the glossary.

\newpage

# Part I — Foundations

## Chapter 1. Pharmacology from zero

### 1.1 The central question

When a person takes a drug, two things happen at once. The body acts on the
drug — absorbing it, moving it around, breaking it down, and excreting it. And
the drug acts on the body — producing the effect we actually care about, whether
that is lowering blood pressure, killing bacteria, or, in our case, thinning the
blood. Pharmacology is the study of both directions of that relationship.

The first direction — *what the body does to the drug* — is
**pharmacokinetics**, abbreviated **PK**. It is usually summarised by the four
letters **ADME**:

- **Absorption** — how the drug gets from the site of administration (a
  swallowed tablet in the gut) into the bloodstream.
- **Distribution** — how it spreads from blood into tissues and organs.
- **Metabolism** — how the body chemically transforms it, mostly in the liver.
- **Excretion** — how it leaves the body, mostly via the kidneys.

The practical output of PK is a single, measurable quantity that changes over
time: the **concentration** of the drug in the blood plasma, written $C_p$. Plot
that concentration against time after a dose and you get the **concentration–time
profile**, the central object of all of PK.

The second direction — *what the drug does to the body* — is
**pharmacodynamics**, abbreviated **PD**. Its output is the **effect**: some
biological signal that responds to the drug. For warfarin the effect is a
measure of how readily the blood can clot.

### 1.2 Why link them?

You could study concentration alone, or effect alone. The power of PK/PD comes
from *linking* them: building a mathematical chain

$$\text{Dose} \;\rightarrow\; \text{Concentration}(t) \;\rightarrow\; \text{Effect}(t).$$

Once you have that chain as equations with estimated numbers in them, you can
answer questions you never directly measured. *What dose gives the right effect?
How fast? What happens in a heavier patient? What if a dose is missed?* This is
the difference between describing data and having a **model** — a compact set of
equations that reproduces the data and, crucially, predicts situations you have
not yet observed. Prediction is the entire point.

### 1.3 A vocabulary you will keep using

A few PK quantities recur so often they deserve early introduction, because the
whole project is ultimately about estimating them:

- **Clearance ($CL$)** — the volume of blood cleared of drug per unit time
  (units like L/h). It is the single most important PK parameter because it
  governs how much drug you need to maintain a target concentration.
- **Volume of distribution ($V$)** — a proportionality constant linking the
  amount of drug in the body to the measured plasma concentration
  ($C_p = A/V$). It is not a literal anatomical volume; it is the volume the
  drug *appears* to occupy.
- **Absorption rate constant ($k_a$)** — how quickly an oral dose moves from
  gut to blood.
- **Half-life ($t_{1/2}$)** — the time for concentration to fall by half, a
  convenient summary derived from $CL$ and $V$.

Hold these lightly for now; Chapter 4 derives them properly.

\newpage

## Chapter 2. The warfarin story

### 2.1 Why this drug is the perfect teacher

Every field has its canonical teaching example. In population PK/PD, one of the
most enduring is warfarin, and for good reasons that are worth understanding
because they *are* the modelling challenge.

Warfarin is an oral **anticoagulant** — a blood thinner — prescribed to prevent
dangerous clots in conditions such as atrial fibrillation, deep-vein thrombosis,
and mechanical heart valves. It has been in clinical use since the 1950s, so it
is exceptionally well characterised.

Three features make it an ideal case study:

1. **A narrow therapeutic window.** Too little warfarin and the patient can
   throw a clot (stroke, embolism); too much and they can bleed dangerously.
   The margin between "not enough" and "too much" is small, which makes precise,
   individualised dosing genuinely important — exactly the problem PK/PD
   modelling exists to solve.
2. **A clear, well-understood mechanism.** Warfarin blocks an enzyme
   (vitamin-K epoxide reductase) that the liver needs to manufacture several
   clotting factors. Less enzyme activity means fewer functional clotting
   factors, which means blood that clots less readily. This is an *indirect*
   mechanism — the drug does not neutralise clotting factors directly; it slows
   their *production*. That single fact dictates the entire shape of the PD
   model, as we will see.
3. **A pronounced delay between concentration and effect.** Because warfarin
   works by throttling the *production* of clotting factors, and because the
   existing factors take time to degrade, the anticoagulant effect lags the
   blood concentration by many hours to days. This delay — visible as a *loop*
   when you plot effect against concentration — is the signature phenomenon the
   model must reproduce.

### 2.2 Measuring the effect: PCA and INR

To model the effect we need to measure it. Two related quantities appear
throughout:

- **PCA — prothrombin complex activity.** A laboratory measure of how active the
  vitamin-K-dependent clotting factors are, expressed as a percentage of
  normal. High PCA means normal clotting; low PCA means the blood is
  anticoagulated. **In the dataset we use, PCA is the measured effect**, so PCA
  is what our PD model predicts.
- **INR — international normalised ratio.** The number used in the clinic to
  monitor warfarin. It is derived from prothrombin time and is, roughly,
  inversely related to PCA. INR is what a physician reads off a chart; PCA is
  the more mechanistic quantity.

A point of scientific honesty that recurs in the project: our data records PCA,
not INR. We therefore model PCA and *discuss* the relationship to INR rather
than inventing INR values we did not measure. Fabricating a variable to look
more clinical would be exactly the kind of shortcut a careful analyst refuses.

### 2.3 The dataset

The project uses the `warfarin` dataset that ships with the `nlmixr2data` R
package — the classic O'Reilly study of 32 healthy subjects who each received a
single oral dose. For every subject it records, over time, the plasma
concentration (`cp`), the effect (`pca`), and covariates: body weight, age, and
sex. Because it is bundled with the software, the analysis is perfectly
reproducible: anyone who installs the package has the identical data, with no
download and no licensing friction.

One limitation shapes a later chapter: the study is **single-dose**. That is
fine for estimating the model, but it means any simulation of *repeated* dosing
is an extrapolation beyond the observed data — a caveat we label honestly rather
than hide.

\newpage

## Chapter 3. The mathematics of pharmacokinetics

This chapter builds the PK model from first principles. The only prerequisite is
comfort with the idea of a rate of change; the calculus is kept gentle and every
step is explained.

### 3.1 The compartment idea

We cannot track every molecule in the body, so we simplify. Imagine the body as
one or more well-stirred **compartments** — think of them as tanks that the drug
flows between. This is a deliberate abstraction: a "central compartment" lumps
together blood and the tissues that equilibrate quickly with it. Astonishingly,
this crude picture predicts real concentration–time data remarkably well, which
is why it has been the backbone of PK for decades.

For an oral drug we need at least two tanks:

- a **depot** compartment representing the gut, where the swallowed dose waits
  to be absorbed, and
- a **central** compartment representing the bloodstream, where we measure
  concentration and from which the drug is cleared.

### 3.2 First-order kinetics

The key modelling assumption is that transfer rates are **first-order**: the
rate at which drug moves out of a compartment is proportional to the amount
currently in it. If $A_{gut}$ is the amount in the gut, absorption proceeds at
rate $k_a A_{gut}$. This is the same mathematics as radioactive decay or
compound interest, and it produces the familiar exponential rise-and-fall of a
drug profile.

Writing $A_{gut}$ and $A_{cen}$ for the amounts in gut and central compartments,
the model is a pair of **ordinary differential equations** (ODEs):

$$\frac{dA_{gut}}{dt} = -k_a\,A_{gut},$$

$$\frac{dA_{cen}}{dt} = k_a\,A_{gut} - \frac{CL}{V}\,A_{cen}.$$

Read them in words. The gut only loses drug, at a rate proportional to what it
holds (first equation). The central compartment gains what the gut gives up and
loses drug to clearance at rate $\frac{CL}{V}$ per unit amount (second
equation). The measured concentration is simply the central amount divided by
the volume:

$$C_p = \frac{A_{cen}}{V}.$$

### 3.3 What the parameters mean, and why we fit on the log scale

The unknowns are $k_a$, $CL$, and $V$. Estimating them from data is the heart of
"fitting a PK model." Two practical points that show up directly in the code:

- All three parameters are **strictly positive** — a negative clearance or
  volume is physically meaningless. To guarantee positivity during fitting, we
  estimate their **logarithms** and exponentiate inside the model
  ($CL = e^{\theta_{CL}}$). The optimiser is then free to search over all real
  numbers without ever producing an impossible value.
- Half-life falls out of the estimates as $t_{1/2} = \ln(2)\,V/CL$. You never
  fit it directly; it is a *derived* quantity.

### 3.4 One compartment or two?

Sometimes a drug leaves the blood in two distinct phases: a rapid initial drop
as it distributes into tissues, then a slower decline as it is cleared. Capturing
that requires a second, "peripheral" tank exchanging drug with the central one —
a **two-compartment** model, with extra parameters for the inter-compartmental
clearance ($Q$) and peripheral volume ($V_2$).

Which structure is right is an empirical question, not an assumption. The
project fits **both** a one- and a two-compartment model and compares them
objectively (Chapter 10). For this warfarin data the one-compartment model wins,
and we let the data say so rather than deciding in advance.

\newpage

## Chapter 4. The mathematics of pharmacodynamics

### 4.1 The naive model, and why it fails for warfarin

The simplest way to link effect to concentration is to assume the effect at any
instant depends only on the concentration at that same instant — a **direct**
model, often an **Emax** curve:

$$E = E_0 - \frac{E_{max}\,C_p}{EC_{50} + C_p}.$$

Here $E_{max}$ is the largest achievable effect and $EC_{50}$ the concentration
producing half of it. Direct models work for many drugs. For warfarin they fail,
and the failure is diagnostic.

If effect depended only on the current concentration, then plotting effect
against concentration would trace a single curve — the same effect value every
time the concentration passed through a given level. Instead, warfarin data
traces a **loop**: at a given concentration the effect is different depending on
whether the concentration is rising or falling. This phenomenon is called
**hysteresis**, and it is unmistakable evidence that the effect *lags* the
concentration. A direct model cannot produce a loop; something with memory —
a delay — is required.

### 4.2 The turnover (indirect-response) model

The delay is not a mathematical trick; it reflects the biology. Warfarin does
not destroy clotting factors, it slows their *synthesis*. The pool of clotting
activity (our PCA) behaves like the water level in a tub with a tap and a drain:

- a tap fills the tub at a **production rate** $k_{in}$ (the body making
  clotting factors),
- a drain empties it at a **first-order loss rate** $k_{out}$ (natural
  degradation),
- at rest the level sits where inflow equals outflow, i.e. baseline
  $PCA_0 = k_{in}/k_{out}$.

Warfarin **turns the tap down**. Because the tub drains only gradually, the
level falls slowly after dosing and recovers slowly as the drug clears — exactly
the observed delay. Mathematically:

$$\frac{d\,PCA}{dt} = k_{in}\cdot\big(1 - I(C_p)\big) - k_{out}\cdot PCA,$$

with the inhibition of production given by an $I_{max}$ relationship

$$I(C_p) = \frac{I_{max}\,C_p}{IC_{50} + C_p}.$$

When concentration is zero there is no inhibition and the tub sits at baseline.
As concentration rises, production is throttled by up to a fraction $I_{max}$,
and PCA falls. This **indirect-response** or **turnover** model is the
scientifically correct structure for warfarin, and choosing it is *derived from
the hysteresis plot*, not assumed. That is the single most important modelling
decision in the entire project, and we reached it by looking at a picture.

### 4.3 Keeping parameters in bounds

As with PK, the turnover parameters must stay physical. Rates ($k_{in}$,
$k_{out}$, $IC_{50}$) are positive, so we fit them on the log scale. $I_{max}$ is
a fraction between 0 and 1, so we fit it on the **logit** scale and squash it
back with the inverse-logit (expit) function. These transformations are not
cosmetic; they stop the optimiser from proposing impossible values and are the
reason the code is littered with `exp()` and `expit()`.

\newpage

## Chapter 5. Population modelling and mixed effects

### 5.1 Why not just fit each person separately?

You could fit the PK model to each of the 32 subjects individually and average
the results. Population modelling does something smarter and more powerful: it
fits **everyone at once**, in a single model that simultaneously estimates the
*typical* parameter values and the *variability* around them. This approach —
**nonlinear mixed-effects (NLME) modelling** — is the standard of the field, and
`nlmixr2` is an R package built to do it.

The payoff is threefold. Sparse individuals (few samples) borrow strength from
the group. You get an explicit, quantified estimate of how much patients differ
from one another. And you can bring in **covariates** — patient characteristics
like body weight — to explain *why* they differ.

### 5.2 Fixed effects, random effects, and residual error

"Mixed effects" means the model mixes two kinds of parameters:

- **Fixed effects** (also called $\theta$, theta) — the population-typical
  values: the typical clearance, the typical volume, and so on. One number each
  for the whole population.
- **Random effects** (also called $\eta$, eta) — each individual's personal
  deviation from the typical value. If the typical log-clearance is
  $\theta_{CL}$, then subject $i$ has $CL_i = e^{\theta_{CL} + \eta_{i}}$. The
  $\eta$'s are assumed drawn from a normal distribution with mean zero and a
  variance the model estimates. That variance *is* the between-subject
  variability.

On top of the individual predictions sits **residual error** — the leftover
mismatch between the model's prediction for a sample and the actual measured
value, capturing assay noise, sampling-time errors, and model imperfection. It
is usually split into an **additive** part (constant size) and a
**proportional** part (grows with the concentration).

So the full hierarchy is: population typical value → individual value (typical
plus that person's $\eta$) → predicted measurement → observed measurement
(prediction plus residual error). Every diagnostic plot later in the book is a
way of interrogating one layer of this hierarchy.

### 5.3 How the fitting actually works (intuition)

Estimating all of this means finding the fixed effects, the random-effect
variances, and the residual-error parameters that make the observed data most
plausible. The engine does it by maximising a likelihood, but the likelihood
involves integrating over every individual's unknown $\eta$'s, which has no neat
closed form. Two families of algorithms handle it:

- **FOCEi** (First-Order Conditional Estimation with interaction) — the classic
  approach, linearising the model around each individual's estimated $\eta$.
- **SAEM** (Stochastic Approximation Expectation–Maximisation) — a
  simulation-based approach that is robust for tricky, nonlinear models. The
  project uses SAEM.

You do not need to implement either; you need to understand that they are doing
maximum-likelihood estimation over a hierarchical model, and that their output
is a set of parameter estimates with uncertainty.

### 5.4 Shrinkage — a concept you must know

When the data on an individual are sparse, the model cannot really tell how that
person differs from the typical value, so their estimated $\eta$ is pulled
("shrunk") toward zero. High **shrinkage** (say above 30%) is a warning that
individual-level estimates and the diagnostic plots based on them are
unreliable. We report shrinkage and interpret it in Chapter 12; for now, file it
as "how much the data failed to inform each individual."

\newpage

# Part II — The project, step by step

## Chapter 6. Project architecture and reproducibility

Before any modelling, the project is organised so that a stranger — or you, six
months later — can reproduce every result by running scripts in order. The
layout separates concerns:

```
warfarin_pkpd_portfolio/
├── data/raw/          # untouched snapshot of the source data
├── data/processed/    # cleaned, analysis-ready data
├── scripts/           # numbered pipeline 00–08
├── models/            # saved fitted model objects (.rds)
├── figures/           # generated plots
├── reports/           # rendered report + result tables
├── setup.R            # installs/loads packages
├── analysis.qmd       # the published report
├── LEARNING_NOTES.md  # beginner walkthrough
└── PUBLISHING.md      # how to publish online
```

Three principles are baked in. **Numbered scripts** encode the order of
operations, each reading what the previous one saved. **Separation of raw and
processed data** means the original is never overwritten, so mistakes are always
recoverable. And **saved model objects** mean expensive fits run once and are
reused by the diagnostics, report, and simulation steps rather than being
recomputed. These are not fussy details; reproducibility is a core professional
expectation in pharmacometrics, where analyses support regulatory submissions
and must be auditable years later.

The remaining chapters in this part walk the pipeline script by script.

\newpage

## Chapter 7. Data cleaning (`01_data_cleaning.R`)

The first script loads the `warfarin` data, takes a raw snapshot, and produces
three analysis-ready datasets: one for PK, one for PD, and one combined for the
joint fit.

The data arrive in **long format**: one row per event, where each row is either
a dose or a single observation, and a column `dvid` says whether an observation
is a concentration (`cp`) or an effect (`pca`). This long, event-based layout is
the universal convention in pharmacometrics (it is how NONMEM and Monolix expect
data too), because it flexibly represents any schedule of doses and mixed
measurement types.

Cleaning here is deliberately light — the dataset is already tidy — but the
script still does the professionally important things: it snapshots the raw data
to `data/raw/`, splits the data into the three analysis sets by filtering on
`evid` (event id: dose vs observation) and `dvid`, runs sanity checks with
`stopifnot()` that will halt the pipeline if an assumption breaks (no missing
times, all doses positive), and saves the results for downstream scripts. The
guiding idea is *fail loudly and early*: a cheap assertion now prevents a
baffling error three scripts later.

\newpage

## Chapter 8. Exploratory analysis (`02_EDA.R`)

You never fit a model to data you have not looked at. Exploratory data analysis
(EDA) builds intuition and, in this project, actually *determines the model
structure*.

### 8.1 Concentration profiles

![Individual warfarin plasma-concentration profiles with the population mean and standard-deviation band. Concentration rises as the oral dose is absorbed, peaks, then declines as the drug is cleared.](figures/02_pk_profiles.png)

Each faint line is one subject; the bold line and shaded band are the population
mean plus or minus one standard deviation. The classic absorb-then-eliminate
shape is exactly what the one-compartment oral model of Chapter 3 produces,
which is our first hint that the structure will fit.

### 8.2 Effect profiles

![Individual PCA (effect) profiles with population mean and SD band. The effect falls after dosing and recovers slowly, hinting at the delayed, indirect mechanism.](figures/02_pd_profiles.png)

The effect does not track the concentration in lockstep. PCA drifts down and
recovers gradually over days — the visual signature of a slow turnover process,
foreshadowing the indirect-response model.

### 8.3 The hysteresis loop — the decisive plot

![Concentration–effect relationship for one subject, coloured by time. The path forms a loop rather than a single curve, proving the effect lags the concentration and mandating a delayed (turnover) PD model.](figures/02_hysteresis.png){width=70%}

This is the plot that settles the modelling strategy. Because the trajectory
forms a loop — the effect at a given concentration differs on the way up versus
the way down — a direct concentration-to-effect model is ruled out and the
turnover model of Chapter 4 is required. We *derived* the model structure from
data, which is exactly how it should be done.

### 8.4 Covariates

![Distributions of body weight and age across the 32 subjects.](figures/02_covariates.png)

Knowing the spread of weight and age tells us which covariate relationships we
can realistically estimate later. A covariate with almost no variability cannot
explain between-subject differences, so this plot is a feasibility check as much
as a description.

Aesthetically, all figures share a single theme defined in
`00_plotting_theme.R` — consistent colours (blue for PK, red for PD), fonts, and
high-resolution export. Consistency is not vanity: a reviewer who can read your
plots quickly trusts your analysis more.

\newpage

## Chapter 9. The PK model (`03_PK_Model.R`)

This script estimates the pharmacokinetic parameters by fitting the
one-compartment oral model and, for rigour, a two-compartment alternative, then
choosing between them objectively.

### 9.1 Writing a model in nlmixr2

An `nlmixr2` model is an R function with two blocks. The `ini({...})` block sets
initial estimates: the typical values (on the log scale), the between-subject
variability terms ($\eta$), and the residual-error parameters. The `model({...})`
block writes the ODEs and the concentration equation almost exactly as they
appear in Chapter 3, then declares the residual-error structure with
`cp ~ add(add.err) + prop(prop.err)`. The closeness of the code to the
mathematics is deliberate and is one of `nlmixr2`'s strengths.

### 9.2 Choosing between structures with AIC

Fitting both models, we compare them with the **Akaike Information Criterion
(AIC)** — a score that rewards good fit but penalises extra parameters, so a more
complex model must earn its complexity with a real improvement. Lower AIC is
better. In this run the one-compartment model scored **1132** against the
two-compartment model's **1287**, so the simpler model wins decisively and the
script saves it as the chosen PK model. Preferring the simplest model that
adequately describes the data is the principle of **parsimony**, and it is
everywhere in pharmacometrics.

### 9.3 Reading the estimates

The fitted typical values were a clearance of about **0.135 L/h**, a volume of
about **7.6 L**, and an absorption rate around **0.59 /h**. Clearance was
estimated very precisely (relative standard error near 3%), while absorption was
fuzzier — unsurprising, since absorption is captured only in the first few hours
after dosing where samples are sparser. The full parameter table appears in the
published report and in Appendix A.

\newpage

## Chapter 10. The PD model (`04_PD_Model.R`)

With PK in hand, this script fits the turnover PD model of Chapter 4 to the
effect data — a **sequential** fit, meaning the concentration is treated as a
known input driving the effect.

### 10.1 A subtlety the error messages taught us

The turnover model needs the drug concentration at *every effect observation
time*. But in long-format data, concentration and effect sit on separate rows —
there is no concentration value sitting next to each PCA measurement. A first
attempt failed with a missing-column error precisely because of this. The fix is
instructive: for each subject we take their observed concentrations and linearly
**interpolate** them onto the effect times, creating the `cp` column the model
needs. This is a genuine feature of sequential PD fitting, not a quirk of our
data, and understanding *why* the error appeared teaches more than the fix
itself.

### 10.2 What the fit shows, honestly

The turnover structure captured the delayed effect, but two parameters — $IC_{50}$
and $I_{max}$ — were estimated very imprecisely, with enormous confidence
intervals. This is expected and worth stating plainly: with single-dose data and
crude interpolated concentrations, the data cannot cleanly separate "how potent"
($IC_{50}$) from "how much maximum inhibition" ($I_{max}$). The honest conclusion
is that the *joint* fit of the next chapter, which uses the model's own
concentration predictions instead of interpolated ones, behaves far better — and
saying so is a strength, not an admission of failure. Reviewers trust an analyst
who flags the weaknesses in their own results.

\newpage

## Chapter 11. The joint PK/PD model (`05_PKPD_Model.R`)

The centrepiece. Here PK and PD are fit **simultaneously**: the PK equations
generate the concentration, that same concentration drives the PD turnover
equations, and every parameter is estimated together against both the
concentration and the effect data at once. The model declares two endpoints —
`cp ~ add + prop` and `pca ~ add` — and `nlmixr2` knows from the data which rows
belong to which.

### 11.1 Adding a covariate: body weight

This is where the population approach earns its keep. Heavier people tend to
clear and hold more drug, so the model lets clearance and volume scale with body
weight using the classic **allometric** form,

$$CL_i = CL_{typ}\left(\frac{WT_i}{70}\right)^{\!\theta_{CL}},\qquad
V_i = V_{typ}\left(\frac{WT_i}{70}\right)^{\!\theta_{V}},$$

with weight normalised to a reference of 70 kg. Rather than fix the exponents at
their textbook values (0.75 for clearance, 1 for volume), we *estimate* them and
let the data speak. They came out at about **0.71** and **1.05** — reassuringly
close to theory, which is a nice internal validation.

### 11.2 Why the joint model is better

Two signs told us the joint model was superior. The exponents landed near their
physiological values, and — more tellingly — the unexplained between-subject
variability on clearance collapsed from around 25% in isolation to about **6%**
once weight was included, because weight now *explains* much of that variability.
A covariate that genuinely matters should reduce the random-effect variance it
speaks to, and here it did. This fitted joint model, saved to `models/`, is the
one we diagnose and simulate.

\newpage

## Chapter 12. Model diagnostics (`06_Model_Diagnostics.R`)

A fitted model means nothing until you have checked it. Diagnostics ask: does the
model actually describe the data, and are its assumptions satisfied?

### 12.1 Goodness of fit

![Observed versus population predictions for both endpoints. Points should scatter evenly around the dashed identity line.](figures/06_dv_vs_pred.png)

If the model is unbiased, observations scatter symmetrically around the line of
identity. Systematic curvature or a lean off the line would signal a structural
problem. Population predictions (using only typical values) scatter more widely
than individual predictions (using each subject's $\eta$), which is expected —
the individual fit should always hug the line more tightly.

### 12.2 Residuals

![Conditional weighted residuals against time. They should sit around zero with no trend; the smooth line is a visual check.](figures/06_cwres_time.png)

**Conditional weighted residuals (CWRES)** are standardised errors: if the model
is right they behave like standard-normal noise, centred on zero with roughly 95%
falling within ±2, and — critically — showing *no trend* against time or
prediction. A drift or a fan shape would reveal something the model is missing,
such as an unmodelled time dependence.

A technical note that, again, an error taught us: SAEM fits do not compute CWRES
by default, so the script must add them explicitly with `addCwres()` before
plotting. The failure message `object 'CWRES' not found` was not a bug in the
data but a reminder that the diagnostic must be requested.

### 12.3 Visual predictive check

![Visual predictive check: many datasets simulated from the model, summarised as prediction intervals, overlaid on the observed data.](figures/06_vpc.png){width=80%}

The **visual predictive check (VPC)** is the most persuasive diagnostic. You
simulate hundreds of virtual trials from the fitted model and check that the real
data fall within the simulated intervals. If they do, the model reproduces both
the central trend and the variability of reality — not just the average, but the
spread. For a portfolio or a regulatory submission, a good VPC is the plot that
convinces.

### 12.4 Shrinkage

The script also prints eta shrinkage. Low shrinkage on clearance in the joint
model means the data genuinely informed each subject's clearance, so the
individual predictions and the plots built on them are trustworthy.

\newpage

## Chapter 13. Simulation (`07_Simulation.R`)

Diagnostics confirm the model is trustworthy; simulation is where it pays off. We
rebuild the fitted structure in `rxode2` (the simulation engine underneath
`nlmixr2`) using the estimated typical parameters, and ask "what if" questions
without running a new experiment.

![Simulated anticoagulant effect over time for two regimens. The shaded band marks an illustrative therapeutic zone; the dotted line marks when the effect first reaches it.](figures/07_dosing_scenarios.png)

Two scenarios are shown: the single 100 mg dose that matches the study (a sanity
check — the simulation should resemble the observed data), and an illustrative 5
mg daily maintenance regimen. The second is explicitly **an extrapolation**: the
data were single-dose, so multi-dose predictions go beyond what we observed, and
we label them as illustrative rather than pretending they are validated. That
label is a mark of integrity, not weakness.

![Concentration and effect on a common timeline for the maintenance regimen, showing how repeated dosing builds concentration while the effect follows with a delay.](figures/07_conc_effect.png){width=70%}

Placing concentration above effect on one timeline makes the central theme of the
whole book visible in a single picture: the effect is a delayed, smoothed
response to the concentration that drives it. The script also computes a
**time-to-target** — when the effect first crosses an illustrative therapeutic
threshold — which is the kind of quantity that, with proper data, directly
informs clinical dosing.

\newpage

## Chapter 14. From analysis to a published report (`08_report_tables.R`, `analysis.qmd`)

The final pipeline step turns results into a shareable artefact. `08_report_tables.R`
exports the fitted parameter tables and the time-to-target as small markdown
files, and `analysis.qmd` is a **Quarto** document that assembles the figures and
tables into a clean, navigable HTML report.

A deliberate design choice makes online publishing painless: the report contains
*no model fitting*. It only embeds pre-computed figures and tables, so it renders
in seconds and needs none of the heavy modelling packages installed on the
server. The project publishes to **Posit Connect Cloud**, which renders the
Quarto document straight from a GitHub repository and serves it at a public URL.
Separating the expensive computation (done once, locally) from the lightweight
presentation (rendered anywhere) is a pattern worth internalising well beyond
this project.

\newpage

# Part III — Context and mastery

## Chapter 15. How this works in the pharmaceutical industry

Everything in Part II is a small, self-contained version of a discipline that
sits at the centre of modern drug development: **pharmacometrics**, the science
of quantitative models of drug behaviour, practised by specialists often called
pharmacometricians or clinical pharmacology modellers.

### 15.1 Model-informed drug development

Regulators and companies now expect quantitative models to inform decisions
across a drug's life, an approach the FDA and EMA call **model-informed drug
development (MIDD)** (the industry also writes MID3). Rather than picking doses
by trial and error, teams build PK and PK/PD models from early data and use them
to choose doses, design the next trial, and justify labelling. The workflow you
just walked — structure selection, covariate modelling, diagnostics, simulation
— is precisely the workflow used on real programmes, only with more data, more
endpoints, and far more scrutiny.

### 15.2 Where models are used along the pipeline

In **preclinical and first-in-human** work, PK models establish how a molecule
behaves and support the very first dose selection. In **early clinical
development**, PK/PD and exposure–response models link concentration to efficacy
and safety, driving dose selection for the pivotal trials — often the single
highest-value decision in a programme. In **late development and submission**,
population PK analyses quantify how age, weight, organ function, and drug
interactions shift exposure, feeding directly into the drug label. After
approval, models support dosing in special populations such as children or
patients with kidney impairment, where running full trials may be impractical.

### 15.3 Deliverables, software, and standards

The concrete outputs are population PK and PK/PD analysis reports, exposure–
response analyses, simulation-based dose justifications, and contributions to
regulatory submission documents. The dominant software has historically been
**NONMEM**, with **Monolix** as a major alternative and open-source tools like
**nlmixr2** (used here) and **Pumas** rising quickly. The modelling logic is the
same across all of them; only the syntax differs, which is why learning the
concepts matters more than learning one tool.

Because these analyses support regulatory decisions, they are held to strict
standards of **reproducibility and traceability** — version-controlled code,
documented data provenance, and auditable results, sometimes under formal
good-practice (GxP) expectations. That is exactly why the project you built
emphasises a raw-data snapshot, numbered scripts, saved model objects, and a
reproducible published report. Those habits are not academic; they are the daily
reality of the job.

### 15.4 The role and the people

A pharmacometrician typically blends pharmacology, statistics, and programming,
works alongside clinical pharmacologists, statisticians, and clinicians, and
translates messy trial data into decisions a drug-development team can act on. It
is an unusually leveraged role: a single well-built model can redirect a
multi-year, multi-hundred-million-dollar programme. Demand consistently outstrips
supply, which is part of why a well-executed portfolio project like this one is a
genuine door-opener.

\newpage

## Chapter 16. Learning pharmacometrics independently

You can learn this field on your own. Here is an honest, sequenced path.

### 16.1 Prerequisites, and how much you actually need

You need less mathematics than you might fear and more *comfort with concepts*
than raw technique. Aim for a working grasp of exponential/first-order processes
and the idea of a differential equation (you do not need to solve them by hand —
the software does), basic probability and regression (means, variances,
distributions, what "fitting" means), and enough R to manipulate data and make
plots. If any of these is shaky, shore it up in parallel rather than waiting
until you feel "ready"; you learn fastest with a concrete project in front of
you — which is the whole point of the one in this book.

### 16.2 A suggested twelve-week plan

A realistic self-study arc, assuming a few focused hours a week:

- **Weeks 1–2 — Pharmacology and PK basics.** Learn ADME, clearance, volume,
  half-life, and one-compartment kinetics until you can sketch a
  concentration–time curve and explain each parameter. Re-read Chapters 1–3
  here.
- **Weeks 3–4 — R and data handling.** Get fluent in reading data, `dplyr`
  wrangling, and `ggplot2`. Reproduce the EDA chapter on the warfarin data
  yourself.
- **Weeks 5–6 — Fitting your first PK model.** Install `nlmixr2`, refit the
  one-compartment model from Chapter 9, and change things deliberately (initial
  values, error model) to see what happens.
- **Weeks 7–8 — PD and PK/PD.** Work through turnover models and the hysteresis
  idea; refit Chapters 10–11 and try altering the covariate model.
- **Weeks 9–10 — Diagnostics and simulation.** Learn to read GOF plots, CWRES,
  and VPCs, and run your own dosing simulations.
- **Weeks 11–12 — Communicate.** Rebuild the report and publish it. Being able
  to *explain* a model clearly is as valued as building it.

### 16.3 Resources worth your time

For textbooks, *Clinical Pharmacokinetics and Pharmacodynamics* by Rowland and
Tozer is the gentle, authoritative starting point; *Pharmacokinetic–
Pharmacodynamic Data Analysis* by Gabrielsson and Weiss is the practical
modelling bench-book; and Bonate's *Pharmacokinetic–Pharmacodynamic Modeling and
Simulation* goes deeper on the statistics. For the tools, the `nlmixr2` and
`rxode2` project websites carry worked examples, and the open pharmacometrics
community publishes tutorials and model libraries. Above all, learn by
**reproducing published analyses** on open datasets — the fastest route from
reading to doing.

### 16.4 Building a portfolio

The project in this book is itself the template: pick a well-characterised drug
with open data, take it end-to-end from EDA to a published report, and — this is
the part most people skip — *write up the reasoning*, not just the code. A
reviewer can see that you can run software; what distinguishes you is showing
that you understood *why* each step was necessary. That is precisely what this
book demonstrates for warfarin, and what you should aim to reproduce for a drug
of your own.

\newpage

## Chapter 17. Questions and answers

A wide-ranging FAQ, from beginner conceptual questions to the specific technical
snags this project hit.

**What is the difference between PK and PD?**
PK is what the body does to the drug (concentration over time); PD is what the
drug does to the body (effect over time). PK/PD links dose to concentration to
effect.

**Why model a population instead of individuals?**
To estimate both the typical response and the variability around it, to let
sparse individuals borrow strength from the group, and to explain differences
between people using covariates. Individual fits give you none of that.

**What exactly is an "eta"?**
$\eta$ is a single individual's deviation from the population-typical value of a
parameter, assumed drawn from a normal distribution with mean zero. Its estimated
variance is the between-subject variability.

**Why fit parameters on the log scale?**
Because clearance, volume, and rate constants must be positive. Estimating their
logarithms lets the optimiser range over all real numbers while the
exponentiation guarantees a positive result.

**Why did the PD model need an indirect-response structure?**
Because the concentration–effect data form a hysteresis loop, proving the effect
lags the concentration. A direct model cannot produce a loop; a turnover model,
which reflects warfarin's inhibition of clotting-factor *production*, can.

**What is hysteresis, in one sentence?**
A different effect at the same concentration depending on whether concentration
is rising or falling — the fingerprint of a delay between concentration and
effect.

**How do you choose between a one- and two-compartment model?**
Fit both and compare with AIC (or the objective function), preferring the
simpler model unless the complex one improves fit enough to justify its extra
parameters. Here the one-compartment model won.

**What is AIC and why penalise complexity?**
The Akaike Information Criterion scores fit while subtracting a penalty for the
number of parameters, guarding against overfitting; lower is better.

**Why were IC50 and Imax estimated so poorly in the sequential PD fit?**
Single-dose data over a limited concentration range, plus interpolated rather
than model-predicted concentrations, cannot separate potency ($IC_{50}$) from
maximum inhibition ($I_{max}$). The joint fit, using model-predicted
concentrations, does much better.

**What does adding body weight as a covariate accomplish?**
It explains part of why people differ, scaling clearance and volume
allometrically. A good covariate reduces the unexplained between-subject
variability — here, clearance variability fell from about 25% to about 6%.

**Why did I get `object 'CWRES' not found`?**
SAEM fits do not compute conditional weighted residuals by default. Call
`addCwres(fit)` before plotting them; the project's diagnostics script now does
this automatically.

**Why did the PD script first fail with a missing `cp` column?**
In long-format data, concentration and effect are on separate rows, so there is
no concentration value beside each effect observation. A sequential PD fit needs
one, so the script interpolates each subject's concentrations onto the effect
times.

**Why did the hysteresis plot throw a `list` / `is.finite` error at first?**
Reshaping the data produced duplicate time points, so the pivot created
list-columns instead of numbers. Collapsing duplicates (averaging) and keeping
only times with both a concentration and an effect fixes it.

**What is a VPC and why is it so valued?**
A visual predictive check overlays observed data on prediction intervals from
many simulated trials. It tests whether the model reproduces both the trend and
the variability of the data — the most convincing single diagnostic.

**What is shrinkage and when should I worry?**
Shrinkage measures how far individual random-effect estimates are pulled toward
zero because the data could not inform them. Above roughly 30%, treat
individual-level estimates and the plots based on them with caution.

**Is PCA the same as INR?**
No. PCA (prothrombin complex activity) is the mechanistic effect measure in our
data; INR is the related clinical index. We model PCA and discuss INR rather than
fabricating INR values.

**Can I trust the multi-dose simulation?**
Treat it as illustrative. The data were single-dose, so repeated-dose
predictions extrapolate beyond what was observed and are labelled accordingly.

**Do I need to know NONMEM to get a job?**
It helps, because much of industry runs on it, but the concepts transfer across
tools. Demonstrating that you understand the modelling — as this project does —
matters more than any single package, and NONMEM syntax is quick to pick up once
the ideas are solid.

**Is `nlmixr2` "good enough" for real work?**
Yes for many purposes, and increasingly used in industry; it implements the same
estimation methods as the commercial tools. Some regulatory settings still expect
NONMEM, but the modelling skill is what transfers.

**How long does it take to become employable in pharmacometrics?**
With focused study and a couple of solid portfolio projects, many people reach
an entry-level, assistant-modeller standard in several months to a year; depth
and judgement then grow with real datasets over years.

**What mathematics do I really need?**
An intuitive grasp of first-order/exponential processes and differential
equations, and basic probability and regression. You do not need to solve ODEs
by hand — the software integrates them.

**What is the single most important idea in this whole project?**
That model *structure* should be derived from the data — the hysteresis loop
telling us to use a turnover model — rather than assumed. Everything else is
mechanics.

\newpage

## Chapter 18. Glossary

**ADME** — Absorption, Distribution, Metabolism, Excretion; the four processes of
pharmacokinetics.

**AIC** — Akaike Information Criterion; a model-comparison score balancing fit
against the number of parameters (lower is better).

**Allometric scaling** — expressing a parameter as a power function of body
size, e.g. $CL \propto (WT/70)^{0.75}$.

**Clearance ($CL$)** — volume of blood cleared of drug per unit time; the primary
determinant of maintenance dose.

**Compartment** — a modelled "tank" representing a kinetically homogeneous part
of the body.

**Covariate** — a subject characteristic (weight, age, sex) used to explain
between-subject variability.

**CWRES** — Conditional Weighted Residuals; standardised errors used to diagnose
model fit.

**Emax model** — a direct concentration–effect model with a maximal effect
$E_{max}$ and a half-maximal concentration $EC_{50}$.

**Eta ($\eta$)** — an individual's random deviation from a population-typical
parameter value.

**Fixed effect ($\theta$)** — a population-typical parameter value.

**FOCEi / SAEM** — algorithms for estimating nonlinear mixed-effects models.

**Half-life ($t_{1/2}$)** — time for concentration to halve; $\ln(2)\,V/CL$.

**Hysteresis** — a loop in the concentration–effect plot indicating the effect
lags concentration.

**IC50 / Imax** — the concentration giving half-maximal inhibition, and the
maximum fractional inhibition, in an inhibitory turnover model.

**INR** — International Normalised Ratio; the clinical index for monitoring
warfarin.

**Indirect-response / turnover model** — a PD model where the drug modulates the
production or loss of a response that turns over with its own kinetics.

**NLME** — Nonlinear Mixed-Effects modelling; the statistical framework of
population PK/PD.

**PCA** — Prothrombin Complex Activity; the measured anticoagulant effect in this
dataset.

**Pharmacometrics** — the science of quantitative models of drug behaviour.

**Residual error** — the unexplained difference between predicted and observed
measurements.

**Shrinkage** — the degree to which individual random-effect estimates are pulled
toward zero when data are uninformative.

**Volume of distribution ($V$)** — apparent volume linking amount of drug to
plasma concentration.

**VPC** — Visual Predictive Check; a simulation-based diagnostic comparing model
predictions to observed data.

\newpage

# Appendix A. Complete code listings

The full source of every pipeline script follows, exactly as used to produce the
results and figures in this book. Read them alongside the corresponding chapters
in Part II.

\newpage

## `setup.R`

```r
# setup.R — install and load everything the project needs.
# Author: Pramod BR
# Date:   2026-07-24
# Run this once before the pipeline: source("setup.R")

pkgs <- c(
  "nlmixr2",      # model fitting (pulls in rxode2, nlmixr2data, etc.)
  "nlmixr2data",  # the warfarin dataset
  "rxode2",       # ODE simulation
  "ggplot2",      # plots
  "dplyr",        # data wrangling
  "tidyr",        # data reshaping
  "xpose",        # goodness-of-fit diagnostics
  "xpose.nlmixr2" # xpose bridge for nlmixr2 fits
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(pkgs, library, character.only = TRUE))

# A place to keep helper settings used across scripts.
set.seed(1234)
theme_set(theme_bw())

message("Setup complete. Packages loaded.")

```

\newpage

## `scripts/00_plotting_theme.R`

```r
# 00_plotting_theme.R
# Author: Pramod BR
# Date:   2026-07-24
# Shared plotting style so every figure in the project looks consistent.
# Sourced by 02, 06 and 07.

library(ggplot2)

# --- Palette ----------------------------------------------------------------
pk_col   <- "#2C6E9B"   # concentration (blue)
pd_col   <- "#B5453B"   # effect / PCA (red)
accent   <- "#E4A700"   # highlights
ink      <- "#1A1A1A"   # text / lines

# --- Theme ------------------------------------------------------------------
theme_pk <- function(base_size = 13) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 3,
                                     hjust = 0, margin = margin(b = 4)),
      plot.subtitle   = element_text(colour = "grey35", size = base_size - 1,
                                     hjust = 0, margin = margin(b = 10)),
      plot.caption    = element_text(colour = "grey55", size = base_size - 3,
                                     hjust = 1, margin = margin(t = 8)),
      axis.title      = element_text(colour = "grey20", face = "bold"),
      axis.text       = element_text(colour = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.4),
      strip.text      = element_text(face = "bold", size = base_size,
                                     margin = margin(4, 4, 4, 4)),
      strip.background = element_rect(fill = "grey95", colour = NA),
      plot.margin     = margin(14, 16, 12, 14),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "top",
      legend.title    = element_text(face = "bold", size = base_size - 2),
      legend.text     = element_text(size = base_size - 2)
    )
}

theme_set(theme_pk())

# --- Helper: consistent, high-resolution export -----------------------------
save_fig <- function(plot, file, width = 7.5, height = 5) {
  ggsave(file, plot, width = width, height = height, dpi = 300,
         bg = "white")
}

```

\newpage

## `scripts/01_data_cleaning.R`

```r
# 01_data_cleaning.R
# Author: Pramod BR
# Date:   2026-07-24
# Load the warfarin data, inspect it, label PK vs PD, save analysis-ready sets.

library(nlmixr2data)
library(dplyr)

data(warfarin)

# Snapshot the raw data so the project is self-contained.
write.csv(warfarin, "data/raw/warfarin_raw.csv", row.names = FALSE)

# Quick look
str(warfarin)
summary(warfarin)

# The dataset is long format. Key columns:
#   id, time, amt, dv, dvid ("cp" = concentration, "pca" = effect),
#   evid (0 = observation, 1 = dose), wt, age, sex.

clean <- warfarin %>%
  mutate(
    dvid = as.character(dvid),
    sex  = as.factor(sex)
  ) %>%
  arrange(id, time)

# Dosing records (shared by both PK and PD models).
dose_records <- clean %>% filter(evid != 0)

# PK analysis set: dosing rows + concentration observations.
pk_data <- clean %>%
  filter(evid != 0 | dvid == "cp")

# PD analysis set: dosing rows + pca observations.
pd_data <- clean %>%
  filter(evid != 0 | dvid == "pca")

# Combined set for the joint PK/PD fit (both endpoints kept).
pkpd_data <- clean

# Basic data-quality checks
stopifnot(!any(is.na(pk_data$time)))
stopifnot(all(dose_records$amt > 0))

# Save processed sets
write.csv(pk_data,   "data/processed/pk_data.csv",   row.names = FALSE)
write.csv(pd_data,   "data/processed/pd_data.csv",   row.names = FALSE)
write.csv(pkpd_data, "data/processed/pkpd_data.csv", row.names = FALSE)
saveRDS(list(pk = pk_data, pd = pd_data, pkpd = pkpd_data),
        "data/processed/analysis_sets.rds")

message("01 done: ", nrow(pk_data), " PK rows, ", nrow(pd_data), " PD rows.")

```

\newpage

## `scripts/02_EDA.R`

```r
# 02_EDA.R
# Author: Pramod BR
# Date:   2026-07-24
# Exploratory plots to understand the data before modelling.

library(dplyr)
library(tidyr)
library(ggplot2)
source("scripts/00_plotting_theme.R")

sets <- readRDS("data/processed/analysis_sets.rds")
obs  <- sets$pkpd %>% filter(evid == 0)

# Mean +/- SD summary at each nominal time, for an overlay on the profiles.
summ <- function(df) {
  df %>% group_by(time) %>%
    summarise(m = mean(dv, na.rm = TRUE), sd = sd(dv, na.rm = TRUE),
              .groups = "drop")
}

# 1. Concentration profiles: faint per-subject lines + mean trend + SD ribbon.
cp <- obs %>% filter(dvid == "cp")
p_pk <- ggplot(cp, aes(time, dv)) +
  geom_line(aes(group = id), colour = pk_col, alpha = 0.18) +
  geom_point(colour = pk_col, alpha = 0.20, size = 1) +
  geom_ribbon(data = summ(cp), aes(time, ymin = pmax(m - sd, 0), ymax = m + sd),
              inherit.aes = FALSE, fill = pk_col, alpha = 0.15) +
  geom_line(data = summ(cp), aes(time, m), inherit.aes = FALSE,
            colour = pk_col, linewidth = 1.1) +
  labs(title = "Warfarin plasma concentration",
       subtitle = "Individual profiles with population mean ± SD",
       x = "Time (h)", y = "Concentration (mg/L)",
       caption = "nlmixr2data::warfarin — 32 subjects, single oral dose")
save_fig(p_pk, "figures/02_pk_profiles.png")

# 2. Effect (PCA) profiles.
pca <- obs %>% filter(dvid == "pca")
p_pd <- ggplot(pca, aes(time, dv)) +
  geom_line(aes(group = id), colour = pd_col, alpha = 0.18) +
  geom_point(colour = pd_col, alpha = 0.20, size = 1) +
  geom_ribbon(data = summ(pca), aes(time, ymin = m - sd, ymax = m + sd),
              inherit.aes = FALSE, fill = pd_col, alpha = 0.15) +
  geom_line(data = summ(pca), aes(time, m), inherit.aes = FALSE,
            colour = pd_col, linewidth = 1.1) +
  labs(title = "Prothrombin complex activity (effect)",
       subtitle = "Effect falls then slowly recovers — note the delay",
       x = "Time (h)", y = "PCA (%)",
       caption = "Lower PCA = greater anticoagulant effect")
save_fig(p_pd, "figures/02_pd_profiles.png")

# 3. Covariate distributions: histogram + density, clean facets.
covs <- sets$pkpd %>% distinct(id, wt, age, sex)
cov_long <- covs %>%
  pivot_longer(c(wt, age), names_to = "cov", values_to = "value") %>%
  mutate(cov = recode(cov, wt = "Weight (kg)", age = "Age (years)"))
p_cov <- ggplot(cov_long, aes(value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12,
                 fill = pk_col, alpha = 0.35, colour = "white") +
  geom_density(colour = ink, linewidth = 0.9) +
  facet_wrap(~ cov, scales = "free") +
  labs(title = "Covariate distributions",
       subtitle = "Demographics across the 32 subjects",
       x = NULL, y = "Density")
save_fig(p_cov, "figures/02_covariates.png", height = 4)

# 4. Hysteresis: effect vs concentration, coloured by time.
#    A loop (not a single line) => effect lags concentration => use an
#    indirect-response (turnover) PD model.
wide <- obs %>%
  filter(!is.na(dv), dvid %in% c("cp", "pca")) %>%
  group_by(id, time, dvid) %>%
  summarise(dv = mean(dv), .groups = "drop") %>%
  pivot_wider(names_from = dvid, values_from = dv) %>%
  filter(!is.na(cp), !is.na(pca)) %>%
  arrange(id, time)

p_hyst <- wide %>%
  filter(id == unique(id)[1]) %>%
  ggplot(aes(cp, pca, colour = time)) +
  geom_path(linewidth = 1, arrow = grid::arrow(length = unit(0.18, "cm"),
                                                type = "closed")) +
  geom_point(size = 2.6) +
  scale_colour_viridis_c(option = "magma", end = 0.9) +
  labs(title = "Concentration–effect hysteresis",
       subtitle = "Subject 1 — the loop is why we need a delayed PD model",
       x = "Concentration (mg/L)", y = "PCA (%)", colour = "Time (h)")
save_fig(p_hyst, "figures/02_hysteresis.png", width = 6.5, height = 5.2)

message("02 done: EDA figures written to figures/")

```

\newpage

## `scripts/03_PK_Model.R`

```r
# 03_PK_Model.R
# Author: Pramod BR
# Date:   2026-07-24
# Fit the population PK model: how the body absorbs and clears warfarin.

library(nlmixr2)
library(dplyr)

sets    <- readRDS("data/processed/analysis_sets.rds")
pk_data <- sets$pk %>% filter(dvid == "cp" | evid != 0)

# One-compartment, first-order oral absorption.
pk_1cmt <- function() {
  ini({
    tka  <- log(0.5)    # absorption rate (ka)
    tcl  <- log(0.13)   # clearance (CL)
    tv   <- log(8)      # volume (V)
    eta.ka ~ 0.3
    eta.cl ~ 0.1
    eta.v  ~ 0.1
    add.err  <- 0.1
    prop.err <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / v) * center
    cp <- center / v
    cp ~ add(add.err) + prop(prop.err)
  })
}

# Two-compartment alternative, to justify the structural choice.
pk_2cmt <- function() {
  ini({
    tka <- log(0.5); tcl <- log(0.13); tv <- log(8)
    tq  <- log(0.5); tv2 <- log(4)
    eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1
    add.err <- 0.1; prop.err <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    q  <- exp(tq)
    v2 <- exp(tv2)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl/v)*center - (q/v)*center + (q/v2)*periph
    d/dt(periph) <-  (q/v)*center - (q/v2)*periph
    cp <- center / v
    cp ~ add(add.err) + prop(prop.err)
  })
}

fit_1 <- nlmixr2(pk_1cmt, pk_data, est = "saem",
                 control = saemControl(print = 0))
fit_2 <- nlmixr2(pk_2cmt, pk_data, est = "saem",
                 control = saemControl(print = 0))

# Compare by objective function / AIC (lower = better, penalising complexity).
cat("1-cmt AIC:", AIC(fit_1), "  2-cmt AIC:", AIC(fit_2), "\n")
pk_fit <- if (AIC(fit_1) <= AIC(fit_2)) fit_1 else fit_2

print(pk_fit)
saveRDS(pk_fit, "models/pk_fit.rds")
message("03 done: PK model saved to models/pk_fit.rds")

```

\newpage

## `scripts/04_PD_Model.R`

```r
# 04_PD_Model.R
# Author: Pramod BR
# Date:   2026-07-24
# Fit the PD model: how warfarin concentration drives the effect (PCA).
# Warfarin inhibits synthesis of clotting factors -> an indirect-response
# (turnover) model with inhibition of the production rate (kin).

library(nlmixr2)
library(dplyr)

sets <- readRDS("data/processed/analysis_sets.rds")

# In a SEQUENTIAL PD fit the model needs the drug concentration (cp) at every
# record. The raw data is long format (cp and pca on separate rows), so we
# build a cp column by linearly interpolating each subject's observed
# concentrations onto all of their PD record times.
conc <- sets$pkpd %>%
  filter(dvid == "cp", evid == 0) %>%
  select(id, time, cp = dv)

pd_data <- sets$pd %>%
  group_by(id) %>%
  group_modify(function(df, key) {
    ci <- conc[conc$id == key$id, ]
    df$cp <- if (nrow(ci) >= 2) {
      approx(ci$time, ci$cp, xout = df$time, rule = 2)$y
    } else {
      0
    }
    df
  }) %>%
  ungroup()
pd_data$cp[is.na(pd_data$cp)] <- 0

# Turnover model with Imax inhibition on production.
#   baseline PCA = kin/kout
#   drug reduces kin via  (1 - Imax*cp/(IC50+cp))
pd_turnover <- function() {
  ini({
    tkin  <- log(1)      # zero-order production rate
    tkout <- log(0.05)   # first-order loss rate
    tic50 <- log(1)      # concentration giving half-maximal inhibition
    timax <- logit(0.9)  # maximum inhibition (0-1)
    eta.kout ~ 0.1
    pca.err <- 5
  })
  model({
    kin  <- exp(tkin)
    kout <- exp(tkout + eta.kout)
    ic50 <- exp(tic50)
    imax <- expit(timax)
    # cp is supplied per record from the observed/PK-predicted concentration.
    inh  <- 1 - (imax * cp) / (ic50 + cp)
    pca(0) <- kin / kout
    d/dt(pca) <- kin * inh - kout * pca
    pca ~ add(pca.err)
  })
}

# Sequential PD fit: uses the interpolated cp column built above.
# (In 05 we fit PK and PD jointly; here we isolate the PD structure.)
pd_fit <- nlmixr2(pd_turnover, pd_data, est = "saem",
                  control = saemControl(print = 0))

print(pd_fit)
saveRDS(pd_fit, "models/pd_fit.rds")
message("04 done: PD model saved to models/pd_fit.rds")

```

\newpage

## `scripts/05_PKPD_Model.R`

```r
# 05_PKPD_Model.R
# Author: Pramod BR
# Date:   2026-07-24
# Joint PK/PD model: fit concentration (cp) and effect (pca) simultaneously,
# and add a covariate (weight on clearance and volume).

library(nlmixr2)
library(dplyr)

sets      <- readRDS("data/processed/analysis_sets.rds")
pkpd_data <- sets$pkpd

pkpd_model <- function() {
  ini({
    # PK
    tka <- log(0.5)
    tcl <- log(0.13)
    tv  <- log(8)
    wt_cl <- 0.75      # allometric-style weight effect on CL
    wt_v  <- 1.0       # weight effect on V
    eta.ka ~ 0.3
    eta.cl ~ 0.1
    eta.v  ~ 0.1
    # PD
    tkin  <- log(1)
    tkout <- log(0.05)
    tic50 <- log(1)
    timax <- logit(0.9)
    eta.kout ~ 0.1
    # residual error
    add.err  <- 0.1
    prop.err <- 0.1
    pca.err  <- 5
  })
  model({
    # PK with weight normalised to 70 kg
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl) * (wt / 70) ^ wt_cl
    v  <- exp(tv  + eta.v)  * (wt / 70) ^ wt_v
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / v) * center
    cp <- center / v

    # PD driven by model-predicted cp
    kin  <- exp(tkin)
    kout <- exp(tkout + eta.kout)
    ic50 <- exp(tic50)
    imax <- expit(timax)
    inh  <- 1 - (imax * cp) / (ic50 + cp)
    pca(0) <- kin / kout
    d/dt(pca) <- kin * inh - kout * pca

    # two endpoints
    cp  ~ add(add.err) + prop(prop.err)
    pca ~ add(pca.err)
  })
}

pkpd_fit <- nlmixr2(pkpd_model, pkpd_data, est = "saem",
                    control = saemControl(print = 0))

print(pkpd_fit)
saveRDS(pkpd_fit, "models/pkpd_fit.rds")
message("05 done: joint PK/PD model saved to models/pkpd_fit.rds")

```

\newpage

## `scripts/06_Model_Diagnostics.R`

```r
# 06_Model_Diagnostics.R
# Author: Pramod BR
# Date:   2026-07-24
# Check whether the joint model fits well: GOF plots, residuals, VPC, shrinkage.

library(nlmixr2)
library(ggplot2)
library(dplyr)
source("scripts/00_plotting_theme.R")

pkpd_fit <- readRDS("models/pkpd_fit.rds")

# SAEM fits don't carry CWRES by default; add them (computed via FOCEi).
if (!"CWRES" %in% names(as.data.frame(pkpd_fit))) {
  pkpd_fit <- addCwres(pkpd_fit)
}

d <- as.data.frame(pkpd_fit) %>%
  mutate(CMT = recode(as.character(CMT),
                      cp = "PK: concentration", pca = "PD: effect (PCA)"))

# --- 1. Observed vs predicted (population and individual) --------------------
gof <- function(xvar, title) {
  ggplot(d, aes(.data[[xvar]], DV)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55",
                linetype = "dashed") +
    geom_point(aes(colour = CMT), alpha = 0.45, size = 1.6) +
    geom_smooth(method = "loess", se = FALSE, colour = ink,
                linewidth = 0.8, formula = y ~ x) +
    facet_wrap(~ CMT, scales = "free") +
    scale_colour_manual(values = c("PK: concentration" = pk_col,
                                   "PD: effect (PCA)" = pd_col),
                        guide = "none") +
    labs(title = title, x = xvar, y = "Observed (DV)")
}

p_pred  <- gof("PRED",  "Observed vs population prediction") +
  labs(subtitle = "Points should scatter evenly around the dashed identity line")
save_fig(p_pred, "figures/06_dv_vs_pred.png", width = 9, height = 4.4)

p_ipred <- gof("IPRED", "Observed vs individual prediction") +
  labs(subtitle = "Individual predictions should hug the identity line closely")
save_fig(p_ipred, "figures/06_dv_vs_ipred.png", width = 9, height = 4.4)

# --- 2. Conditional weighted residuals --------------------------------------
cwres_plot <- function(xvar, title) {
  ggplot(d, aes(.data[[xvar]], CWRES)) +
    geom_hline(yintercept = 0, colour = "grey55") +
    geom_hline(yintercept = c(-2, 2), colour = "grey80", linetype = "dotted") +
    geom_point(aes(colour = CMT), alpha = 0.45, size = 1.6) +
    geom_smooth(method = "loess", se = FALSE, colour = accent,
                linewidth = 1, formula = y ~ x) +
    scale_colour_manual(values = c("PK: concentration" = pk_col,
                                   "PD: effect (PCA)" = pd_col),
                        name = NULL) +
    labs(title = title, y = "CWRES")
}

p_cwres_t <- cwres_plot("TIME", "Conditional weighted residuals vs time") +
  labs(subtitle = "Should sit around 0 with no trend (yellow = loess)",
       x = "Time (h)")
save_fig(p_cwres_t, "figures/06_cwres_time.png", width = 8, height = 4.4)

p_cwres_p <- cwres_plot("PRED", "Conditional weighted residuals vs prediction") +
  labs(subtitle = "Should sit around 0 with no trend", x = "Population prediction")
save_fig(p_cwres_p, "figures/06_cwres_pred.png", width = 8, height = 4.4)

# --- 3. Eta shrinkage --------------------------------------------------------
print(pkpd_fit$shrink)

# --- 4. Visual predictive check ---------------------------------------------
vpc_res <- tryCatch({
  v <- vpcPlot(pkpd_fit, n = 300, show = list(obs_dv = TRUE)) +
    theme_pk() +
    labs(title = "Visual predictive check",
         subtitle = "Observed data should fall within the simulated intervals")
  save_fig(v, "figures/06_vpc.png", width = 8, height = 5)
  "VPC written"
}, error = function(e) paste("VPC skipped:", conditionMessage(e)))
message(vpc_res)

message("06 done: diagnostics written to figures/")

```

\newpage

## `scripts/07_Simulation.R`

```r
# 07_Simulation.R
# Author: Pramod BR
# Date:   2026-07-24
# Use the fitted model to simulate dosing scenarios with rxode2.
# NOTE: source data is single-dose, so multi-dose results are illustrative.

library(nlmixr2)
library(rxode2)
library(ggplot2)
library(tidyr)
library(dplyr)
source("scripts/00_plotting_theme.R")

pkpd_fit <- readRDS("models/pkpd_fit.rds")

# Build an rxode2 model matching the fitted structure, using the fitted
# typical-value (fixed-effect) parameter estimates.
sim_model <- rxode2({
  cl <- tcl * (WT / 70) ^ 0.75
  v  <- tv
  ka <- tka
  d/dt(depot)  <- -ka * depot
  d/dt(center) <-  ka * depot - (cl / v) * center
  cp <- center / v
  inh <- 1 - (imax * cp) / (ic50 + cp)
  pca(0) <- kin / kout
  d/dt(pca) <- kin * inh - kout * pca
})

# Pull typical values from the fit (back-transformed).
fe  <- fixef(pkpd_fit)
tka <- exp(fe[["tka"]]); tcl <- exp(fe[["tcl"]]); tv <- exp(fe[["tv"]])
kin <- exp(fe[["tkin"]]); kout <- exp(fe[["tkout"]]); ic50 <- exp(fe[["tic50"]])
imax <- plogis(fe[["timax"]])

params <- c(tka = tka, tcl = tcl, tv = tv,
            kin = kin, kout = kout, ic50 = ic50, imax = imax, WT = 70)

# Scenario A: single 100 mg dose (matches the study).
evA <- et(amt = 100, cmt = "depot") %>% et(seq(0, 240, by = 1))
simA <- rxSolve(sim_model, params, evA)

# Scenario B: 5 mg daily maintenance for 10 days (illustrative extrapolation).
evB <- et(amt = 5, cmt = "depot", ii = 24, addl = 9) %>% et(seq(0, 336, by = 1))
simB <- rxSolve(sim_model, params, evB)

target <- 25  # illustrative "therapeutic" PCA level

plot_df <- bind_rows(
  transform(as.data.frame(simA), scenario = "A: single 100 mg"),
  transform(as.data.frame(simB), scenario = "B: 5 mg daily x10 (illustrative)")
)

# First time PCA drops to/below target, per scenario (for annotation).
ttp <- plot_df %>%
  group_by(scenario) %>%
  filter(pca <= target) %>%
  summarise(t = ifelse(n() > 0, min(time), NA_real_), .groups = "drop")

# --- Plot 1: effect (PCA) over time, with therapeutic threshold -------------
p_sim <- ggplot(plot_df, aes(time, pca)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = target,
           fill = pd_col, alpha = 0.06) +
  geom_hline(yintercept = target, colour = pd_col, linetype = "dashed") +
  geom_line(colour = pd_col, linewidth = 1) +
  geom_vline(data = ttp, aes(xintercept = t), colour = accent,
             linetype = "dotted", na.rm = TRUE) +
  facet_wrap(~ scenario, scales = "free_x") +
  annotate("text", x = Inf, y = target, label = paste0("  target ", target, "%"),
           hjust = 1, vjust = -0.5, size = 3.2, colour = pd_col) +
  labs(title = "Simulated anticoagulant effect over time",
       subtitle = "Shaded band = therapeutic zone; dotted line = time to reach it",
       x = "Time (h)", y = "PCA (%)",
       caption = "Typical-value simulation. Scenario B extrapolates beyond single-dose data.")
save_fig(p_sim, "figures/07_dosing_scenarios.png", width = 9.5, height = 4.6)

# --- Plot 2: concentration and effect together (scenario B) -----------------
long_B <- as.data.frame(simB) %>%
  select(time, Concentration = cp, PCA = pca) %>%
  pivot_longer(-time, names_to = "series", values_to = "value")
p_pkpd <- ggplot(long_B, aes(time, value, colour = series)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(Concentration = pk_col, PCA = pd_col),
                      guide = "none") +
  labs(title = "Concentration drives effect (Scenario B)",
       subtitle = "Repeated dosing builds concentration; effect follows with a delay",
       x = "Time (h)", y = NULL)
save_fig(p_pkpd, "figures/07_conc_effect.png", width = 8, height = 5.4)

message("Time to reach PCA<=", target, "%:")
print(ttp)

saveRDS(plot_df, "data/processed/simulation_results.rds")
message("07 done: simulation figures written to figures/")

```

\newpage

## `scripts/08_report_tables.R`

```r
# 08_report_tables.R
# Author: Pramod BR
# Date:   2026-07-24
# Export fitted-model results as markdown tables so the published report can
# embed them WITHOUT re-fitting (fast, reliable render on Posit Connect Cloud).
# Run this locally after 03-07, then commit reports/results/ and figures/.

library(nlmixr2)

dir.create("reports/results", recursive = TRUE, showWarnings = FALSE)

# Turn a fit's parameter table into a tidy data frame.
tidy_params <- function(fit) {
  pf <- as.data.frame(fit$parFixedDf)
  pf <- cbind(Parameter = rownames(pf), pf)
  rownames(pf) <- NULL
  num <- vapply(pf, is.numeric, logical(1))
  pf[num] <- lapply(pf[num], function(x) signif(x, 3))
  pf
}

# Write a data frame as a GitHub-style markdown table with a caption.
write_md <- function(df, file, caption) {
  md <- knitr::kable(df, format = "pipe", caption = caption, row.names = FALSE)
  writeLines(md, file)
  message("wrote ", file)
}

pk_fit   <- readRDS("models/pk_fit.rds")
pd_fit   <- readRDS("models/pd_fit.rds")
pkpd_fit <- readRDS("models/pkpd_fit.rds")

write_md(tidy_params(pk_fit),   "reports/results/pk_params.md",
         "PK model — population parameter estimates")
write_md(tidy_params(pd_fit),   "reports/results/pd_params.md",
         "PD model — population parameter estimates")
write_md(tidy_params(pkpd_fit), "reports/results/pkpd_params.md",
         "Joint PK/PD model — population parameter estimates")

# Time-to-target from the simulation step.
sim <- readRDS("data/processed/simulation_results.rds")
target <- 25
tt <- do.call(rbind, lapply(split(sim, sim$scenario), function(s) {
  hit <- s$time[s$pca <= target]
  data.frame(Scenario = unique(s$scenario),
             `Time to PCA <= 25% (h)` =
               if (length(hit)) min(hit) else NA_real_,
             check.names = FALSE)
}))
write_md(tt, "reports/results/time_to_target.md",
         "Time to reach the illustrative therapeutic target (PCA <= 25%)")

message("08 done: report tables written to reports/results/")

```

\newpage

## `analysis.qmd` (published report source)

```yaml
---
title: "Warfarin Population PK/PD Analysis"
author: "Pramod BR"
date: "2026-07-24"
format:
  html:
    theme: cosmo
    toc: true
    toc-depth: 2
    toc-location: left
    number-sections: true
    embed-resources: true
    fig-align: center
    code-tools: false
---

::: {.callout-note appearance="simple"}
A reproducible population PK/PD analysis of warfarin built in R with
**nlmixr2**. Figures and parameter tables are pre-computed by the pipeline
(`scripts/01`–`08`) and embedded here, so this report renders quickly without
re-fitting the models.
:::

# Objective

Describe the pharmacokinetics of warfarin (how the body absorbs and clears it),
link concentration to its anticoagulant effect — prothrombin complex activity
(**PCA**) — and use the fitted model to simulate dosing scenarios.

Warfarin is an ideal case study: a narrow therapeutic window, a well-understood
mechanism (inhibition of clotting-factor synthesis), and a clear *delay*
between concentration and effect that the model must capture.

# Data

Source: the `warfarin` dataset from **`nlmixr2data`** (O'Reilly study — 32
subjects, single oral dose, plasma concentration `cp` and effect `pca`, with
weight, age and sex). It ships with the package, so nothing is downloaded.

# Exploratory analysis

Individual profiles with the population mean ± SD. Concentration rises then
falls; effect falls then slowly recovers.

::: {layout-ncol=2}
![Plasma concentration](figures/02_pk_profiles.png)

![Effect (PCA)](figures/02_pd_profiles.png)
:::

The concentration–effect relationship forms a **loop**, not a single line —
effect lags concentration. This hysteresis is why the PD model must be a
delayed, indirect-response (turnover) model rather than a direct one.

![Concentration–effect hysteresis](figures/02_hysteresis.png){width=70%}

![Covariate distributions](figures/02_covariates.png)

# PK model

A one-compartment model with first-order oral absorption was selected over a
two-compartment alternative (lower AIC). Parameters: absorption rate (`ka`),
clearance (`CL`) and volume (`V`).

{{< include reports/results/pk_params.md >}}

# PD model

An indirect-response **turnover** model: warfarin inhibits production (`kin`)
of the clotting-factor pool via an Imax relationship, so PCA declines with a
delay and recovers as the drug clears.

{{< include reports/results/pd_params.md >}}

# Joint PK/PD model

PK and PD fit simultaneously, with body weight added on clearance and volume.
The weight exponents came out close to the classic allometric values, and
between-subject variability on clearance dropped sharply versus the separate
fits.

{{< include reports/results/pkpd_params.md >}}

# Diagnostics

Goodness-of-fit and residual plots confirm the joint model describes both
endpoints well.

::: {layout-ncol=2}
![Observed vs population prediction](figures/06_dv_vs_pred.png)

![CWRES vs time](figures/06_cwres_time.png)
:::

![Visual predictive check](figures/06_vpc.png){width=80%}

# Simulation

Using the fitted typical-value parameters, we simulate a single 100 mg dose
(matching the study) and an illustrative 5 mg daily maintenance regimen. The
shaded band marks an illustrative therapeutic zone.

![Dosing scenarios](figures/07_dosing_scenarios.png)

![Concentration drives effect, with a delay](figures/07_conc_effect.png){width=70%}

{{< include reports/results/time_to_target.md >}}

# Conclusions

This workflow demonstrates an end-to-end population PK/PD analysis: data
handling, exploratory analysis, structural and covariate model building,
diagnostics, and simulation. The measured readout is PCA; INR relates to PCA
clinically but is not measured in this dataset, so it is not reported here. The
maintenance-dose simulation extrapolates beyond the single-dose data and is
labelled as illustrative.

::: {.callout-tip appearance="simple"}
Full method and beginner-friendly explanation: see `LEARNING_NOTES.md`.
Reproduce everything by running `scripts/01`–`08` in order.
:::

```
