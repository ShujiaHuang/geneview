# Geneview Plotting Architecture: why `_core` and `genometracks` share `plotstyle` but not each other

> Audience: contributors who want to understand how geneview's two plotting
> entry points (single-axis statistical plots vs. the multi-track genome
> browser) divide responsibilities and why they are not merged. Read alongside
> [`base_plotting_engine_guide.md`](./base_plotting_engine_guide.md).

## 1. One-sentence conclusion

`geneview._core` (single-axis decorator) and `geneview.genometracks`
(multi-track orchestrator) are **two sibling consumers of `plotstyle`**, each
serving a different plot topology. They reuse the same **style foundation**
(`plotstyle`) but **do not reuse each other's drawing / axes-creation logic** —
because `_core`'s "one function → one Axes" abstraction does not map onto
genometracks' "many tracks → many Axes" orchestration.

This is by design, not an oversight: both already reuse what is genuinely
shared (`plotstyle`), while each stays focused on its own rendering paradigm.

## 2. Layered view

```
                     plotstyle (shared style foundation)
               use_style / get_style / PlotStyle
               ├─ apply_to_axes(ax)     ← single-axis: spine / tick visibility
               └─ to_track_params()     ← multi-track: colours/fonts/widths per track
                    /                              \
      geneview._core (single-axis)        geneview.genometracks (multi-track)
      @styled_plot + get_or_create_axes    plot_tracks + Track class hierarchy
      manhattan / qq / venn /              GenomeAxisTrack / AnnotationTrack /
      admixture / karyoplot                DataTrack / GeneRegionTrack / ...
```

## 3. Why both reuse `plotstyle`

`plotstyle` is the topology-agnostic **style authority**: journal themes
(nature/science/cell), fonts, line widths, export settings, spine-visibility
policy, and so on. Both entry points need those conventions, so both enter
style via the `use_style(...)` context manager, ensuring the rcParams apply
during figure creation and drawing.

The only difference is **how the style is landed onto axes** — which is exactly
where the topology diverges (see below).

## 4. Why they don't reuse each other

| Aspect          | `_core`                                | `genometracks`                                  |
| --------------- | -------------------------------------- | ----------------------------------------------- |
| Topology        | one function → **one** Axes            | one Figure → **GridSpec multi-panel** (per track) |
| Organisation    | functional (`@styled_plot` decorator)  | object-oriented (`Track` hierarchy, each `plot()`) |
| Axes creation   | `get_or_create_axes()` makes **one** ax | GridSpec builds an N×2 panel matrix (data + title) |
| Style landing   | `PlotStyle.apply_to_axes(ax)`          | `PlotStyle.to_track_params()` pushed to each track |
| Spine policy    | hides top/right by default (framed plot) | each track owns its spines (often all hidden / per-type) |

Key points:

* **Axes creation does not generalise.** `get_or_create_axes()` only produces a
  single simple Axes, whereas genometracks needs a GridSpec panel matrix; the
  decorator's figure-creation branch is meaningless to it.
* **The spine policy is opposite.** `_core`'s default of hiding top/right is
  *wrong* for genome tracks — track panels frequently hide all spines or
  customise them per track type, and forcing the policy would break rendering.
* **Style lands differently.** Single-axis uses `apply_to_axes(ax)`;
  multi-track uses `to_track_params()` to push colours/fonts/widths to each
  track as a *floor* (values the user set explicitly are not overridden).

## 5. When to use which

**Use `_core` (`@styled_plot`):** a **statistical plot** where one function
draws one kind of data onto one Axes — e.g. manhattan, qq, venn, admixture,
karyoplot. To add such a function, follow the
[base plotting engine guide](./base_plotting_engine_guide.md).

**Use `genometracks` (`plot_tracks`):** a genome-browser-style view that
**vertically stacks multiple tracks on shared genomic coordinates** (axis,
annotation, gene models, coverage, alignments, lollipop, ...). To add a new
track type, subclass `Track` and implement its drawing protocol rather than
applying `@styled_plot`.

## 6. Call chains

### `_core` single-axis path

```
manhattanplot(data, ..., style="nature", ax=None)
  └─ @styled_plot wrapper
       ├─ use_style("nature")                       # enter the plotstyle context
       ├─ get_or_create_axes(ax=None, figsize=...)  # make a single Axes
       ├─ PlotStyle.apply_to_axes(ax)               # enforce spines when apply_spines=True
       └─ call the body to draw on ax → return ax
```

### `genometracks` multi-track path

```
plot_tracks([GenomeAxisTrack(), AnnotationTrack(...), DataTrack(...)], style="nature")
  ├─ get_style("nature") → resolved_style
  ├─ _apply_style_to_tracks(tracks, style)          # to_track_params() pushes per-track params
  └─ with use_style(resolved_style):                # enter the same plotstyle context
       ├─ _plot_full_layout(...)  → GridSpec builds an N×2 panel matrix
       └─ call each track's plot() to draw on its own Axes → return the axes list
```

Both chains meet `plotstyle` at the **`use_style(...)` entry point**, then
diverge into two non-shared paths: "single-axis `apply_to_axes`" and
"multi-track `to_track_params` + GridSpec".

## 7. Infrastructure layer: why `plotstyle` and `palette` stay separate

The infrastructure layer holds `plotstyle` and, alongside it, a parallel
`palette` module. They are **not merged**; each carries a distinct concern:

* **`palette` = colour data + generation utility (a value / primitive layer):**
  `xkcd_rgb` (named-colour catalogue), `circos`, `CYTOBAND_COLORS` /
  `get_cytoband_color` (genomics-domain constants), and
  `generate_colors_palette` (colormap → colour list). All topology-agnostic and
  style-agnostic.
* **`plotstyle` = journal style policy (a policy layer):** fonts, line widths,
  spine visibility, rcParams, the `PlotStyle` registry. A palette is just one of
  its many attributes, and a *curated* short ordered list (Wong / Okabe-Ito /
  Cell) — different in purpose from `palette`'s general catalogue.

Why `palette` is *not* folded into `plotstyle`:

* **The dependency runs the other way.** `plotstyle` **never imports**
  `palette` — each journal style hardcodes its palette inline as a local
  constant (e.g. `_WONG_PALETTE` in `_nature.py`). `palette` is a zero-dependency
  leaf primitive, so merging would remove no duplication.
* **Most consumers don't go through `plotstyle`.** `palette` is consumed
  directly by `karyotype`, `genometracks` (ideogram cytobands), `baseplot`
  (venn), and `_core`; if it lived inside `plotstyle`, those modules would have
  to depend on the whole style registry just to fetch a colour constant.
* **Different abstraction levels.** One is topology/journal-agnostic colour
  *data*; the other is holistic figure *policy*.

```
Infrastructure:  plotstyle (style policy)        palette (colour primitives)
                      ↑                                ↑
        ┌─────────────┴──────────┬─────────────────────┴─────────────┐
     _core (single-axis)   genometracks (multi-track)   karyotype / baseplot ...
```

> Optional reverse tweak: the curated journal palettes *could* move into
> `palette` as named palettes and be referenced by `plotstyle`, but that is a
> minor "where the constant lives" refactor with limited payoff — not required.

## 8. See also

* [Base plotting engine guide (`geneview._core`)](./base_plotting_engine_guide.md)
* Source: `geneview/_core/` (decorators / canvas / colors),
  `geneview/genometracks/_track_plot.py` (`plot_tracks`),
  `geneview/plotstyle/` (`use_style` / `PlotStyle`)
