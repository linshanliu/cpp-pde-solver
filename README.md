# Finite-Difference PDE Option Pricer

A C++ library that prices options by solving the Black–Scholes partial
differential equation numerically on a grid, using finite differences.

It is built to be extended: the same solver handles vanilla calls and puts
today, and is structured so that local-volatility dynamics or barrier features
can be added later without rewriting the core. Every result is checked against
an independent closed-form formula.

---

## Table of contents

1. [The problem in one paragraph](#1-the-problem-in-one-paragraph)
2. [The idea behind the method](#2-the-idea-behind-the-method)
3. [How the code is organised](#3-how-the-code-is-organised)
4. [Walking through one pricing run](#4-walking-through-one-pricing-run)
5. [The three design decisions that matter](#5-the-three-design-decisions-that-matter)
6. [The parts a reviewer should look at](#6-the-parts-a-reviewer-should-look-at)
7. [Known limitations](#7-known-limitations)
8. [Building and running](#8-building-and-running)

---

## 1. The problem in one paragraph

An **option** is a contract whose payoff depends on the future price of some
asset (a stock, say). A European call option pays `max(S − K, 0)` at a fixed
future date `T`: if the stock `S` ends above the strike `K`, you collect the
difference; otherwise you get nothing. The question this library answers is:
**what is that contract worth today**, before we know where the stock will end
up?

There is a famous closed-form answer for the simplest case (the Black–Scholes
formula). But the moment the contract gets slightly more complicated — a
barrier that knocks the option out if the stock touches some level, a
volatility that changes with price — the formula stops existing. So instead of
a formula, we solve the underlying *equation* numerically. That is what this
library does. Black–Scholes is used here as the test case precisely *because*
the formula exists, so we can check the numerical answer against it.

---

## 2. The idea behind the method

### 2.1 The price obeys a PDE

The theory (a result called **Feynman–Kac**) says the option price `V`, viewed
as a function of the stock price `S` and time `t`, satisfies a partial
differential equation:

```
∂V/∂t  +  ½ σ²S² ∂²V/∂S²  +  (r−q)S ∂V/∂S  −  rV  =  0
```

- `σ` (sigma) is volatility — how much the stock jumps around,
- `r` is the interest rate, `q` the dividend yield,
- the three terms with derivatives are, in order: **diffusion** (spreading),
  **convection** (drift), and **reaction** (discounting).

That "diffusion + convection + reaction" structure is a standard shape in
physics (it describes heat spreading in a moving fluid). The code leans on that
analogy — the base class is literally called `ConvectionDiffusionPDE` — because
all the numerical tricks for that family of equations are well known.

We know one thing for certain: the value **at the expiry date** `t = T` is just
the payoff, `max(S − K, 0)`. So we know the answer at the end, and the PDE tells
us how it propagates backwards. The whole method is: **start from the payoff at
expiry, and march backwards in time to today.**

### 2.2 From a PDE to a grid

We can't solve for `V` at every possible `(S, t)` — there are infinitely many.
So we lay down a **grid**: a finite set of stock prices `S₀, S₁, …, S_M` and a
finite set of time points `t₀, …, t_N`. We only ever compute `V` at grid
points.

To step backwards from one time layer to the previous one, we need the
derivatives `∂V/∂S` and `∂²V/∂S²`. On a grid, a derivative is approximated by
**differences between neighbouring points** — this is the "finite difference"
in the name. For example, on an evenly spaced grid:

```
∂²V/∂S²  at point i   ≈   (V[i−1] − 2·V[i] + V[i+1]) / h²
```

where `h` is the spacing. Plug those approximations into the PDE and it turns
into a system of linear equations linking the grid values at one time layer to
the next. Solving that system, layer by layer, from expiry back to today, gives
the price.

### 2.3 The one subtlety worth knowing: implicit vs explicit

When we step from one time layer to the previous one, we have a choice.

- **Explicit:** compute the new layer directly from the old one. Simple, but
  *unstable* unless the time step is tiny — if you take steps too big, rounding
  errors blow up and the answer becomes garbage.
- **Implicit:** the new layer is defined by an equation that mixes its own
  values together, so you have to *solve* a linear system at each step. More
  work per step, but stable at any step size.

This library implements both, plus the in-between (**Crank–Nicolson**), through
a single dial called `θ` (theta):

| `θ`   | Scheme          | Trade-off                                   |
|-------|-----------------|---------------------------------------------|
| `0`   | Explicit        | Cheapest per step, needs tiny steps         |
| `0.5` | Crank–Nicolson  | Most accurate, but can "ring" (see §5.2)    |
| `1`   | Fully implicit  | Rock-solid stable, less accurate            |

The linear system to be solved at each implicit step is **tridiagonal** (each
grid point only talks to its two neighbours), which can be solved very cheaply
by the **Thomas algorithm** — that's `ThomasSolver`.

### 2.4 What the three schemes actually do, step by step

To make the trade-offs concrete, here is what one backward time step looks like
in each case. Write `L` for the operation "apply the finite-difference
derivatives at every grid point" — it takes the whole layer of values and
returns the PDE's right-hand side. We're stepping from a known layer `Vᵒˡᵈ`
(closer to expiry) to an unknown one `Vⁿᵉʷ` (one step closer to today).

The general rule, with the dial `θ`, is:

```
Vⁿᵉʷ  =  Vᵒˡᵈ  +  dt · [ θ · L(Vⁿᵉʷ)  +  (1−θ) · L(Vᵒˡᵈ) ]
```

Read it as: the change over one step is a blend of the derivatives evaluated at
the *old* layer and at the *new* layer, mixed in proportion `(1−θ) : θ`.

**Explicit (`θ = 0`).**

```
Vⁿᵉʷ  =  Vᵒˡᵈ  +  dt · L(Vᵒˡᵈ)
```

Everything on the right is known, so `Vⁿᵉʷ` is just read off directly — one
multiply-add per grid point, no system to solve. Fast and simple.

The catch is **stability**. Each step effectively takes a step of size
proportional to `dt / h²` along the diffusion. If that ratio is too big, the
grid points start *overcorrecting* each other: a small error at one point
overshoots at its neighbour on the next step, which overshoots back even
harder, and the whole solution oscillates and explodes. Concretely you need

```
dt  ≤  h² / (σ²·Sₘₐₓ²)          (roughly)
```

so halving the grid spacing `h` forces *quartering* the time step. Refining
space to get accuracy makes the explicit scheme brutally expensive. This is why
the old convergence test used `N ∝ M²` time steps — that ratio is not a choice,
it's the stability limit.

**Fully implicit (`θ = 1`).**

```
Vⁿᵉʷ  −  dt · L(Vⁿᵉʷ)  =  Vᵒˡᵈ
```

Now `Vⁿᵉʷ` appears on both sides, tangled together through `L` (each point's new
value depends on its neighbours' new values). You can't read it off; you have to
**solve** for the whole layer at once. Because each point only couples to its
two neighbours, that "solve" is a tridiagonal system — cheap, handled by
`ThomasSolver` in one pass.

The payoff for the extra work is **unconditional stability**: any `dt` you like,
however large, and it never blows up. The overcorrection problem simply can't
happen, because every point is solved in mutual consistency with its neighbours
rather than reacting to their stale values. The price is accuracy: it's only
**first-order** in time (halving `dt` roughly halves the error), because
evaluating the derivatives entirely at the new layer is a slightly lopsided
approximation of the step.

**Crank–Nicolson (`θ = 0.5`).**

```
Vⁿᵉʷ  −  ½ dt · L(Vⁿᵉʷ)  =  Vᵒˡᵈ  +  ½ dt · L(Vᵒˡᵈ)
```

The average of the two. Still an implicit solve (a tridiagonal system, same cost
as fully implicit), but by centring the derivatives halfway between the old and
new layers it becomes **second-order** accurate — halving `dt` quarters the
error. For smooth problems it's the sweet spot: implicit stability *and* high
accuracy.

The "but" — and it's the whole of §5.2 — is that Crank–Nicolson is stable in the
sense that it won't *explode*, but it doesn't actively *damp* the roughest
wobbles either; it leaves them bouncing at constant size, flipping sign each
step. On a smooth solution there are no such wobbles and this never bites. But a
payoff has a sharp corner, which is rough, so those wobbles get excited and show
up as "ringing" — hence the fixes in §5.2 (put a node on the corner; or take a
few fully-implicit steps first to damp the corner, then switch to
Crank–Nicolson).

**Summary of the mechanism:**

| Scheme          | One step is…                    | Cost/step      | Error   | Stability                         |
|-----------------|---------------------------------|----------------|---------|-----------------------------------|
| Explicit        | a direct formula                | cheapest       | `O(dt)` | only if `dt ≲ h²` — very limiting |
| Crank–Nicolson  | a tridiagonal solve (centred)   | one solve      | `O(dt²)`| stable, but doesn't damp roughness|
| Fully implicit  | a tridiagonal solve (new layer) | one solve      | `O(dt)` | stable at any `dt`, damps roughness|

All three are the *same code path* in `ThetaScheme` with a different `θ` — which
is exactly why they live in one class instead of three.

---

## 3. How the code is organised

The design principle is **separation of concerns**: each class knows one thing
and nothing else. A physicist's equation, a contract's payoff, the shape of the
grid, and the linear-algebra solver are four unrelated ideas, so they are four
unrelated classes. You can swap any one without touching the others.

```
ConvectionDiffusionPDE   ← abstract: "a PDE of the diffusion/convection shape"
      └── BlackScholesPDE ← the specific coefficients for Black–Scholes

Payoff                   ← abstract: "what does this contract pay at expiry?"
      ├── VanillaPayoff   ← max(S−K,0) etc.
      └── DigitalPayoff   ← pays 1 if in the money, else 0

BoundaryCondition        ← abstract: "what is V at the edges of the grid?"
      ├── ZeroBC
      ├── DiscountedConstantBC
      └── LinearityBC

Grid                     ← abstract: "where are the grid points?"
      ├── UniformGrid     ← evenly spaced
      └── SinhGrid        ← clustered around a chosen price

FDMBase                  ← the engine: owns the solution, marches through time
      └── ThetaScheme     ← the actual explicit/CN/implicit stepping rule

ThomasSolver             ← pure linear algebra, knows nothing about finance
AnalyticBlackScholes     ← the closed-form formula, used ONLY to check answers
```

The key relationships:

- **`FDMBase` owns nothing but the answer.** The PDE, the payoff, the boundary
  conditions and the grid are all handed to it from outside (dependency
  injection). It just orchestrates: apply the payoff at expiry, then step
  backwards.
- **`ThomasSolver` is deliberately ignorant.** Its header says: *"if the word
  'option' ever appears in this file, something has gone wrong."* It solves a
  tridiagonal system; it has no idea what the numbers mean. That isolation is
  what lets you trust it.
- **`AnalyticBlackScholes` never touches the solver.** It exists purely as an
  independent yardstick. If the exact formula leaked into the solver, checking
  the solver against it would be circular.

---

## 4. Walking through one pricing run

Here is what `main.cpp` does, end to end, to price a one-year call:

```cpp
BlackScholesPDE    pde(r, sigma, q);          // the equation
VanillaPayoff      payoff(K, OptionType::Call);// the contract
Domain             dom{ 0.0, 4.0 * K, T };    // solve for S in [0, 400], t in [0,1]
LinearityBC        lo, hi;                     // edge behaviour
BoundaryConditions bc{ lo, hi };
UniformGrid        grid(dom, 400, 2000);       // 400 price steps, 2000 time steps

ThetaScheme fd(pde, payoff, bc, grid, 0.5);    // Crank–Nicolson
fd.Solve();                                    // march from expiry back to today
double price = fd.Price(100.0);                // read off the value at S=100
```

Under the hood, `Solve()`:

1. **`ApplyTerminal()`** — fills the grid at `t = T` with the payoff.
2. Loops backwards over time layers. For each step, **`ThetaScheme::StepBack`**:
   - asks the PDE for its coefficients at each grid point (`AssembleOperator`),
   - builds the tridiagonal system for this step,
   - patches the two edge rows using the boundary conditions,
   - hands it to `ThomasSolver` to solve,
   - overwrites the solution with the result.
3. Returns the solution at `t = 0`. `Price(S)` interpolates it at the `S` you
   asked for.

`main.cpp` then prints the finite-difference price, delta and gamma next to the
closed-form values, so you can see them agree.

---

## 5. The three design decisions that matter

These are the parts that separate this from a copy-paste textbook solver. Each
is a place where the "obvious" approach is subtly wrong.

### 5.1 The boundary at `S = 0` needs no boundary condition

Every grid has edges, and normally you must specify what the solution does
there. But the **lower** edge, `S = 0`, is special. Look at the Black–Scholes
coefficients:

```
diffusion  = ½σ²S²  →  0   at S = 0
convection = (r−q)S →  0   at S = 0
```

Both vanish. The PDE **degenerates** into a simple equation that involves no
`S`-derivatives at all — it just says the value discounts at rate `r`. In other
words, **the equation determines its own value at that edge**; you are not
allowed to impose anything.

This is a known result (the *Fichera* condition for degenerate PDEs), and the
code exploits it: the `LinearityBC` placed at `S = 0` automatically produces `0`
for a call and `K·e^(−r(T−t))` for a put, with no special-casing. That is why
`main.cpp` uses the *same* boundary class at both ends and it just works.

> This only holds because the lower edge is at `S = 0`. For a barrier option
> where the lower edge is at some `B > 0`, the coefficients don't vanish and you
> *do* need a real boundary condition — which is what `ZeroBC` is for.

### 5.2 Crank–Nicolson can "ring", and the grid can fix it

The payoff `max(S − K, 0)` has a sharp corner at `S = K`. Sharp corners are hard
for numerical methods: Crank–Nicolson (the most accurate scheme) can produce
spurious oscillations near that corner — small wobbles that shouldn't be there.
They show up most visibly in **gamma** (the second derivative), which magnifies
any error.

There are two defences, both present in the code:

1. **Put a grid point exactly on the corner.** `SinhGrid` clusters points around
   a chosen `focus` and pins one exactly there. For a vanilla option you focus
   on `K`; the corner is resolved instead of straddled.
2. **Rannacher start-up** (documented in `ThetaScheme`): take a couple of fully
   implicit steps first — which smooth out the corner — then switch to
   Crank–Nicolson. Because `θ` is just a parameter, this is a two-line change,
   not a new class. *(This is described in the code comments as the intended
   extension; the current `main` runs plain Crank–Nicolson.)*

The `SinhGrid` clustering matters even more for **barriers**: if the barrier
level isn't exactly on a grid point, the error there pollutes the whole
solution and accuracy halves. So `SinhGrid` guarantees the focus lands on a
node.

### 5.3 Greeks come free from the grid

**Delta** (sensitivity of price to the stock) and **gamma** (sensitivity of
delta) are what traders actually hedge with. The naive way to get them is
*bump-and-revalue*: reprice at `S+ε` and `S−ε` and take the difference. That
means solving the whole PDE three times, and you have to pick `ε` carefully.

This library instead reads delta and gamma straight off the grid, reusing the
same finite-difference weights that built the PDE (`FDMBase::Delta`,
`FDMBase::Gamma`). Zero extra cost — the solution is already sitting on the grid
— and consistent with how the price was computed.

---

## 6. The parts a reviewer should look at

If you want to judge the quality of the code, these are the files that carry the
real thinking, in order:

- **`Grid.cpp` / `Grid.hpp`** — the finite-difference weights are derived for a
  *non-uniform* grid (the standard formulas only work on evenly spaced grids;
  getting this right is where most homemade solvers quietly go wrong). The
  `SinhGrid` construction places the focus exactly on a node.
- **`ThetaScheme.cpp`** — the time-stepping, including how the boundary rows are
  folded into the tridiagonal system, and the degenerate `S = 0` handling.
- **`ThomasSolver.hpp`** — the header comment documents *when* the solver is
  guaranteed stable (a condition called the cell Péclet number), which is a
  statement about the method's limits, not just its behaviour.
- **`BlackScholesPDE.hpp`** — the class comments explain the design boundaries:
  why a barrier is *not* a subclass of the PDE, and why the interface carries a
  time argument the current code doesn't use (so local volatility drops in
  later).
- **`AnalyticBlackScholes.hpp`** — the independent benchmark, and the discipline
  of keeping it out of the solver path.

---

## 7. Known limitations

Stated plainly, because knowing a method's limits is part of using it:

- **Convection dominance near `S = 0`.** The scheme is only guaranteed
  oscillation-free when the cell Péclet number `|b|·h / a < 2`. For
  Black–Scholes this can fail at the very smallest stock prices, where the
  option is worth ≈ 0 anyway, so it causes no visible harm. The proper fix
  (upwind differencing) is noted but not implemented.
- **Digitals need special care.** `DigitalPayoff` has a *jump*, not a corner.
  Putting the jump on a grid node actually *hurts* accuracy here; the correct
  remedy (payoff averaging) is not yet in the code.
- **Calendar / American features.** The solver prices European options.
  American early-exercise and multi-asset PDEs are natural extensions the
  architecture is built for but does not yet include.
- **`ConvergenceTest.cpp` predates the current interface** and does not compile
  against it; it is being rewritten to run a proper convergence study (halving
  the grid spacing and confirming the error drops as expected).

---

## 8. Building and running

The repository currently ships a Visual Studio solution (`PDEProject.sln`).
To build with Visual Studio, open the solution and build the `PDEProject`
target. `main.cpp` prints a table comparing the finite-difference price, delta
and gamma against the closed-form Black–Scholes values.

A cross-platform CMake build is planned so the project builds on Linux and macOS
without Visual Studio.
