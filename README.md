# Signals & Systems for TI-Nspire CX

A signals and systems toolkit for the TI-Nspire CX / CX II (CAS or non-CAS). You type signals the way you'd write them on paper, like `e^(-2t)u(t)`, `(1/2)^n u[n]` or `x(2t+1)`, and it works out the answer with the steps shown.

It doesn't use the calculator's CAS, so it runs the same on non-CAS models.

## What it does

- **Signal analyzer:** exact form, piecewise breakdown, shifting/scaling of a defined signal, even and odd parts, periodicity and fundamental period, energy and power
- **Running integral / running sum**
- **Definite integrals and sums**, including impulses (sifting)
- **Convolution**, continuous and discrete, worked interval by interval
- **System properties:** memoryless, time-invariant, linear, causal, BIBO stable, invertible, each with a reason or counterexample
- **LTI systems from h(t) or h[n]:** causality, stability, step response, output for a given input. You can also start from the step response.
- **Fourier transform, DTFT, Laplace and Z-transform**, with regions of convergence for Laplace and Z
- **Inverse transforms** using partial fractions; delays, impulses and band-limited (rectangular) spectra are supported
- **Frequency response** from a differential or difference equation, or from a system like `y(t)=x(t)+x(t-1)`, with magnitude and phase
- **Fourier series:** continuous and discrete
- **Sampling and aliasing**
- **Cascade and parallel systems**
- **Convolution property check:** solves the same convolution directly and with transforms and compares the results
- **Complex numbers:** rectangular and polar forms
- **Graphing:** impulses are drawn as arrows and discrete signals as stem plots. You can trace, zoom and pan.

There are 140+ worked example problems you can open and run, a set of reference tables (transform pairs, properties, identities), and a search that covers all of it.

## Installing

Copy `dist/SignalsSystems.tns` to your calculator with TI-Nspire Computer Link or the Student Software, then open it from My Documents. OS 3.2 or newer is needed.

See [GUIDE.md](GUIDE.md) for a full walkthrough.

## Typing things in

| What | Type |
|---|---|
| step, impulse | `u(t)`, `d(t)` (or `delta(t)`) |
| discrete time | `u[n]`, `d[n-2]`, `x[n-1]` |
| exponentials | `e^(-2t)`, `(1/2)^n` |
| imaginary unit | `j` or `i` |
| pulses and windows | `rect(t)`, `tri(t)`, `(abs(t)<1)`, `(0<=t)(t<2)` |
| integral, sum | `int(f, tau, a, b)`, `sum(f, k, a, b)`, use `inf` for ∞ |
| derivative | `der(x(t))` or `x'(t)` |
| convolution inline | `conv(f, g)` |
| frequency variables | `jw` for jω, `e^(-jW)` for e^(−jΩ), `s`, `z` |

Use the **Definitions** screen to set things like `x(t)=...`, `h[n]=...` or `a=2`. Every tool can then use them, for example `x(2t+1)` or `h(t-1)`.

## Keys

- **Home:** type to search, or move up and down and press enter
- **In a tool:** tab or up/down to change field, enter to solve, esc to go back
- **Results:** up/down scrolls, left/right pages, G opens the graph
- **Graph:** left/right to trace, up/down to switch plots, `+`/`-` to zoom, `4 6 8 2` to pan, `A` to auto-scale, `P` for phase

## Limitations

- **Closed-form answers:** the solver works with sums of t^k·e^(st) pieces and impulses. That covers most textbook signals.
- **Numeric-only results:** sinc functions and Gaussians get a numeric answer instead of an exact one.
- **Symbolic constants:** a constant like `a` has to be given a number in Definitions.
- **System properties:** these are worked out from the expression and also checked numerically with test signals.

## Building from source

The code lives in `src/` as plain Lua 5.1. To build:

```
./build.sh
```

This writes `dist/SignalsSystems.lua` (the whole app in one file) and, if [Luna](https://github.com/ndless-nspire/Luna) is installed, `dist/SignalsSystems.tns`. Set `LUNA=/path/to/luna` if it isn't on your PATH. You can also paste `dist/SignalsSystems.lua` into the Script Editor in the TI-Nspire Student Software.

To run every example, and the set of hand-checked answers, on a computer:

```
cat tests/stubs.lua dist/SignalsSystems.lua tests/run_examples.lua > /tmp/t.lua
lua5.1 /tmp/t.lua

cat tests/stubs.lua dist/SignalsSystems.lua tests/known_answers.lua > /tmp/t.lua
lua5.1 /tmp/t.lua
```

Every transform and inverse transform also checks its own answer: it computes the value numerically at a test point and compares it with the formula. The result shows "(OK)" at the bottom, or a "Check FAILED" line if something is off.

### Source layout

| File | Contents |
|---|---|
| `01_core.lua` | complex numbers, number formatting |
| `02_parse.lua` | expression parser and printer |
| `03_eval.lua` | numeric evaluation |
| `04_signal.lua`, `05_sigio.lua` | exact signal representation, convolution, integration, printing |
| `06_poly.lua` – `08_inverse.lua` | polynomials, forward and inverse transforms |
| `09b_check.lua` | numeric self-checks for transform answers |
| `09_ctx.lua`, `10_tools_a.lua` – `12b_tools_d.lua` | the tools |
| `13_ui.lua` | screens, keys, graphing |
| `14_examples.lua`, `15_ref.lua` | example problems and reference tables |

## License

MIT, see [LICENSE](LICENSE).
