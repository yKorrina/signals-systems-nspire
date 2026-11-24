# User Guide

## Getting started

1. Send `SignalsSystems.tns` to the calculator and open it.
2. The home screen lists everything. Arrow down to a tool and press enter, or just start typing to search.
3. Fill in the boxes and press enter. The answer shows up with the steps underneath.
4. Press esc to go back and change your input.

If you're not sure where to start, open **Example problems**, pick one, and look at how the input is written. You can edit it and run it again.

## Moving around

| Where | Keys |
|---|---|
| Home | up/down, enter to open, type letters to search |
| Input screen | tab or up/down to switch boxes, left/right to move the cursor, enter to solve, esc to go home |
| Results | up/down to scroll, left/right to jump a page, G for the graph, esc to edit |
| Graph | left/right to trace, up/down for the next plot, + and - to zoom, 4 6 8 2 to pan, A to reset, P for phase |

The menu key also has every tool, sorted into groups.

## Writing signals

Write things close to how you would on paper. Use parentheses for continuous time and brackets for discrete time.

```
e^(-2t)u(t)
3cos(5t+pi/2)
(1/2)^n u[n]
d(t-1) + 2d(t-2)      d is the impulse δ
rect(t/2)  tri(t)
(abs(t)<1)            1 when |t|<1, 0 otherwise
(0<=t)(t<2)           1 for 0 ≤ t < 2
t(u(t)-u(t-1))
```

Other things you can use:

- `j` or `i` for the imaginary unit, `pi`, `inf`
- `int(f, tau, a, b)` for an integral, `sum(f, k, a, b)` for a sum
- `der(x(t))` or `x'(t)` for a derivative
- `conv(f, g)` to convolve right inside an expression

For transforms, write `jw` for jω, `W` for Ω (so `e^(-jW)`), and `s` or `z` as usual.

## Definitions

The **Definitions** screen is where you name signals and constants once and reuse them everywhere. Put one per line:

```
x(t)=(t+2)(-2<=t)(t<0)+2(0<=t)(t<1)
h[n]=(1/2)^n u[n]
a=3
```

After that you can type `x(2t+1)`, `x(-t)u(t)`, `h[n-1]` or `e^(-a t)u(t)` in any tool. Loading an example may replace your definitions with the ones that example needs.

## The tools

**Signal analyzer:** Type any signal, like `x(2t+1)`. You get the simplified form, a piecewise listing, how the shift and scale move the graph, whether it's even or odd (with both parts), whether it's periodic, and its energy or power.

**Running integral / sum:** Finds ∫ from -∞ to t (or Σ up to n), interval by interval.

**Integral / sum evaluator:** Definite integral or sum between two limits. Impulses are handled by sifting.

**Convolution:** Enter x and h. It flips and shifts h, finds each overlap interval, and does each integral or sum term by term.

**System properties:** Enter something like `y(t)=x(t-1)+t x(t)` or `y[n]=x[-n]`. It checks memoryless, time-invariant, linear, causal, stable and invertible. Each one gets a reason, and when a property fails you get a counterexample.

**LTI system from h:** Give h (or the step response). It tells you if the system is causal, stable and memoryless, and finds the step response. If you also give an input it computes the output.

**Interconnected systems:** Cascade or parallel h1 and h2.

**Fourier transform, DTFT, Laplace, Z-transform:** Work out the transform piece by piece. Laplace and Z also give the ROC.

**Inverse transforms:** Enter X. It does partial fractions, shows the coefficient for each pole, and turns each term into a time signal. For inverse Laplace and Z you can give a ROC like `>-1`, `<2` or `-2<x<-1`. If you leave it blank, it assumes the signal is causal.

**Frequency response H:** Enter a differential or difference equation, like `y''+3y'+2y=x` or `y[n]-1/2y[n-1]=x[n]`, or a system like `y(t)=x(t)+x(t-1)`. You get H, its magnitude and phase, the impulse response, and the output for an input if you give one.

**Fourier series:** Give one period and the period length. Works for continuous and discrete signals.

**Sampling & aliasing:** Give x(t) and Ts (or `fs=10`). Shows the Nyquist check, x[n], which frequencies alias, and what an ideal reconstruction gives back.

**Convolution property check:** Solves the same convolution two ways (directly, and by multiplying transforms) and checks that they match.

**Complex numbers:** Converts between rectangular and polar form and shows the steps.

**Graph signals:** Plots up to three signals at once.

## Reference tables

Transform pairs and properties, convolution rules, impulse rules, system property tests, trig identities, series and sampling, all in one place on the calculator. Search finds them too.

## When something goes wrong

- **"Unknown symbol":** a constant like `a` doesn't have a value yet. Add `a=2` (or whatever value) in Definitions.
- **"Not closed-form":** signals like sin(t)/t or Gaussians can't be solved exactly. The Fourier tool still gives numbers and a plot.
- **"is not defined":** you used `x(...)` or `h(...)` without defining it first.
- **Nothing seems to happen:** some tools, mainly system properties, can take a few seconds on the calculator.
