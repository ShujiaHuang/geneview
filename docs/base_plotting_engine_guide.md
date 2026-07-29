# Geneview Base Plotting Engine (`geneview._core`) — Developer Guide

> Audience: contributors who want to add a new plotting function to geneview,
> or refactor an existing one, on top of the shared base plotting engine.

## 1. Why an engine?

Historically every top-level plotting function in geneview (`manhattanplot`,
`qqplot`, `admixtureplot`, `venn`, `karyoplot`, ...) re-implemented the same
boilerplate:

```python
def someplot(data, ..., ax=None, style=None, **kwargs):
    # 1. validate inputs
    with use_style(style):                 # enter the style context
        if ax is None:                     # create a canvas
            _, ax = subplots(figsize=(...))
        ...                                # draw
        ax.spines["top"].set_visible(False)  # spine bookkeeping (sometimes)
        ax.set_title(...); ax.set_xlabel(...)
    return ax
```

That "enter style → get axes → enforce spines" lifecycle was copied (slightly
differently) into each module, drifting over time: some functions hid the
top/right spines, some did not; some used `plt.gca()`, some used
`subplots(...)`; `karyoplot` did not even support the `style=` keyword.

`geneview._core` collects that lifecycle into a **single, tested engine** so
that:

* every plot function behaves consistently,
* new plot functions are a few lines of drawing code, and
* engine-wide improvements (e.g. a new default) land everywhere at once.

The engine is **purely additive** — the public API of every migrated function
is unchanged (`karyoplot` merely *gained* an optional `style=` keyword).

## 2. Where it sits in the architecture

```
┌─────────────────────────────────────────────────────────────┐
│  Plot modules:  gwas / popgene / baseplot / karyotype / ...   │  ← draw
├─────────────────────────────────────────────────────────────┤
│  geneview._core   (styled_plot, get_or_create_axes, colors)   │  ← engine
├─────────────────────────────────────────────────────────────┤
│  Infrastructure:  plotstyle / palette / utils                 │  ← shared
└─────────────────────────────────────────────────────────────┘
```

`_core` depends **only** on `plotstyle` and `palette`; it never imports a plot
module, so there are no import cycles. Plot modules import from `.._core`.

> **`_core` vs. `genometracks`:** this engine serves "one function → one Axes"
> statistical plots; `geneview.genometracks` (the multi-track genome browser)
> is a separate, sibling entry point that **does not reuse** `_core` but
> **shares the same `plotstyle` foundation**. For why they divide up this way,
> plus their applicable scenarios and call chains, see
> [Plotting architecture: `_core` and `genometracks`](./plotting_architecture_design.md).

## 3. Public building blocks

Import everything from the package root:

```python
from geneview._core import styled_plot, get_or_create_axes, color_cycle, resolve_colors
```

### 3.1 `@styled_plot(...)` — the decorator (start here)

Wraps a plotting function with the shared style + canvas lifecycle. Inside a
decorated function you may assume **`ax` is always a real `Axes`** and the
**style context is already active** — you only write the drawing code.

```python
@styled_plot(figsize=(9, 3), subplot_kws={"facecolor": "w", "edgecolor": "k"})
def manhattanplot(data, ..., ax=None, style=None, **kwargs):
    ax.scatter(...)          # ax is guaranteed non-None
    ax.set_xlabel(...)
    return ax                # always return the Axes
```

What the decorator does, in order, on every call:

1. reads the `style` keyword and enters `use_style(style)` (a no-op when
   `style is None`, so the currently active global style still applies);
2. reads the `ax` keyword; when it is `None`, creates a figure/axes via
   `get_or_create_axes` using `figsize` / `subplot_kws`;
3. when `apply_spines=True` (the default), enforces the active style's
   spine-visibility flags on the axes (`PlotStyle.apply_to_axes`);
4. calls your function body;
5. if the body raises **and** the decorator created the figure, closes that
   figure so invalid input never leaks half-drawn figures.

The public signature and docstring are preserved via `functools.wraps`, so
callers and Sphinx see no difference.

**Decorator parameters**

| Parameter      | Default | Meaning |
| -------------- | ------- | ------- |
| `figsize`      | `None`  | Default `(w, h)` when the function creates its own figure. `None` falls back to the active style's `figure_figsize`. |
| `apply_spines` | `True`  | Enforce the style's spine rules. Set `False` for plots that manage their own frame (bar plots, axis-off diagrams). |
| `use_gca`      | `False` | When `True`, an omitted `ax` reuses `plt.gca()` instead of a new figure (legacy `karyoplot` behaviour). |
| `subplot_kws`  | `None`  | Extra kwargs for `plt.subplots` (`facecolor`, `constrained_layout`, ...). |

**Requirement:** the decorated function *must* expose an `ax` keyword argument
(the decorator raises `TypeError` otherwise). A `style` keyword is honoured
when present; functions without one simply always use the active style.

### 3.2 `get_or_create_axes(...)` — the canvas helper

The reusable replacement for `if ax is None: _, ax = subplots(...)`. Call it
directly when a function needs finer control than the decorator gives (for
example an *internal* draw helper that must also work when called standalone —
see `_draw_admixtureplot`).

```python
def get_or_create_axes(ax=None, *, figsize=None, style=None,
                       apply_spines=False, use_gca=False, **subplot_kws) -> Axes
```

* `ax` provided → returned as-is (spines optionally enforced);
* `ax is None`, `use_gca=False` → `plt.subplots(figsize=..., **subplot_kws)`;
* `ax is None`, `use_gca=True` → `plt.gca()`;
* `figsize is None` → falls back to the resolved style's `figure_figsize`;
* `style is None` → uses the currently active style (so it is correct when
  called inside a `use_style(...)` context, e.g. from within `@styled_plot`).

Note the defaults differ from the decorator: here `apply_spines` defaults to
`False` (a low-level helper stays neutral), while the decorator defaults to
`True` (the common case for a finished plot).

### 3.3 `color_cycle(color)` — colour cycling

Returns an infinite iterator, reproducing the historical manhattan/qq rules
exactly:

| Input                     | Result |
| ------------------------- | ------ |
| `"#3B5488,#53BBD5"`       | cycles `["#3B5488", "#53BBD5"]` |
| `"rb"` (no comma)         | cycles the **characters** `"r", "b", ...` |
| `["r", "g", "b"]`         | cycles the list elements |

```python
colors = color_cycle(color)
for group in groups:
    c = next(colors)
```

### 3.4 `resolve_colors(palette, n_colors, alpha=1.0)` — palette → colours

Thin, shared wrapper over `geneview.palette.generate_colors_palette`. Accepts a
colormap name, an explicit list, or a `Colormap`, and returns a list of colours
(possibly shorter than `n_colors` if the palette cannot supply enough — check
the length and warn, as `admixtureplot` does).

## 4. Writing a new plotting function (recipe)

```python
# geneview/<subpackage>/_myplot.py
import numpy as np
from .._core import styled_plot, color_cycle   # + get_or_create_axes / resolve_colors as needed


@styled_plot(figsize=(6, 4), subplot_kws={"facecolor": "w"})
def myplot(data, ax=None, color="#3B5488,#53BBD5",
           title=None, xlabel=None, ylabel=None, style=None, **kwargs):
    """One-line summary.

    Parameters
    ----------
    ...
    style : str, PlotStyle, or None, optional
        Plot style to apply. A registered style name ("nature", "science",
        "cell"), a PlotStyle object, or None to use the active style.
    ax : matplotlib axis, optional
        Axis to plot on; a new one is created when omitted.

    Returns
    -------
    ax : matplotlib Axes
    """
    # 1. Validate inputs (raise early with a clear message).
    #    Tip: cheap validation can happen before drawing; the decorator will
    #    close the auto-created figure if you raise here.

    # 2. Draw. `ax` is guaranteed non-None; the style context is active.
    colors = color_cycle(color)
    ax.plot(...)

    # 3. Titles / labels.
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax
```

Then export it from the sub-package `__init__.py` and the top-level
`geneview/__init__.py`, exactly as the existing functions do.

### Conventions

* **Always** accept `ax=None` and `style=None`, and **always** return the `Axes`.
* Keep `**kwargs` forwarding to the underlying matplotlib call (e.g. `scatter`).
* Put domain logic (chromosome offsets, λ inflation, clustering, Venn geometry)
  in the module — the engine is deliberately generic and knows nothing about
  genomics.
* Add tests under `geneview/tests/test_<name>.py`.

## 5. `apply_spines` decision guide

The one real decision when decorating is whether the style should own the
spines. The migrated functions set the precedent:

| Function        | `apply_spines` | Why |
| --------------- | -------------- | --- |
| `manhattanplot` | `True`         | Finished framed plot; the geneview default hides top/right (matches the old hand-written spine code). |
| `qqplot` / `qqnorm` | `True`     | Same framed-plot contract; now consistent with the active style. |
| `admixtureplot` | `False`        | Draws its own box: sets a custom line width on **all four** spines. |
| `venn` (`vennx`) | `False`       | Calls `ax.set_axis_off()`; spine visibility is irrelevant. |
| `karyoplot`     | `False` + `use_gca=True` | Preserves the legacy `plt.gca()` fallback and its custom axis; does not want top/right hidden. |

Rule of thumb: **`True`** for a standard framed x/y plot; **`False`** when the
function draws its own frame or turns the axis off.

## 6. Internal draw helpers that must also run standalone

Some modules split a public entry point from an internal `_draw_*` worker that
is exercised directly by unit tests (e.g. `_draw_admixtureplot`). The decorator
only wraps the *public* function, so the worker still needs to acquire an axes
when called on its own. Route that through `get_or_create_axes` with the same
defaults, so it is a no-op when the decorator already supplied an axes:

```python
def _draw_admixtureplot(..., ax=None):
    ax = get_or_create_axes(ax, figsize=(14, 2), apply_spines=False,
                            facecolor="w", constrained_layout=True)
    ...
```

## 7. Behaviour notes / compatibility

* **API unchanged.** All migrated functions keep their signatures and return
  values. `karyoplot` only *gained* an optional `style=` keyword.
* **`qqplot`/`qqnorm` spines:** previously these left all four spines visible
  (matplotlib default); they now follow the active style and hide top/right,
  matching `manhattanplot` and the journal styles. This is the intended,
  consistent behaviour, not an API change.
* **Figure sizes preserved.** Each decorated function passes its historical
  default `figsize`, so standalone (no-`ax`) output is visually unchanged
  (manhattan `9x3`, qq `5x5`, admixture `14x2`, venn `7x7`).
* **No leaked figures.** If validation raises after the decorator created a
  figure, the decorator closes it.

## 8. Testing

Run the migrated modules' tests plus the full suite:

```bash
python -m pytest geneview/tests/test_manhattan.py geneview/tests/test_qq.py \
    geneview/tests/test_venn.py geneview/tests/test_admixture.py \
    geneview/tests/test_karyotype.py -q
python -m pytest geneview/tests/ -q
```

## 9. File map

```
geneview/_core/
├── __init__.py       # public exports
├── canvas.py         # get_or_create_axes
├── colors.py         # color_cycle, resolve_colors
└── decorators.py     # styled_plot
```
