# Geneview 绘图架构设计：`_core` 与 `genometracks` 为何共享 `plotstyle` 却不互相复用

> 读者对象：想理解 geneview 两套绘图入口（单轴统计图 vs 多轨道基因组浏览器）
> 如何分工、为何不合并的贡献者。配合
> [`base_plotting_engine_guide.zh.md`](./base_plotting_engine_guide.zh.md) 阅读。

## 1. 一句话结论

`geneview._core`（单轴装饰器）与 `geneview.genometracks`（多轨道编排器）是
**`plotstyle` 的两个并列消费者**，服务于两种不同的图形拓扑。它们复用同一个
**样式底座**（`plotstyle`），但**不互相复用绘图/建轴逻辑**——因为
`_core` 的"单函数→单 Axes"抽象无法映射到 genometracks 的"多轨道→多 Axes"编排上。

这不是遗漏，而是符合设计：两者已经复用了真正共享的东西（`plotstyle`），
而各自专注于不同的绘图范式。

## 2. 分层视图

```
                    plotstyle（共享样式底座）
              use_style / get_style / PlotStyle
              ├─ apply_to_axes(ax)      ← 单轴：应用脊柱/刻度可见性
              └─ to_track_params()      ← 多轨：下发颜色/字号/线宽到每个 track
                   /                              \
      geneview._core（单轴装饰器）        geneview.genometracks（多轨编排器）
      @styled_plot + get_or_create_axes   plot_tracks + Track 类层次
      manhattan / qq / venn /             GenomeAxisTrack / AnnotationTrack /
      admixture / karyoplot               DataTrack / GeneRegionTrack / ...
```

## 3. 为什么都复用 `plotstyle`

`plotstyle` 是与图形拓扑无关的**样式规范中心**：期刊主题（nature/science/cell）、
字体、线宽、导出参数、脊柱可见性策略等。两套入口都需要这些规范，因此都以
`use_style(...)` 上下文管理器进入样式，保证 rcParams 在建图与绘制期间生效。

区别只在于**如何把样式落到坐标轴上**——这正是拓扑差异所在（见下）。

## 4. 为什么不互相复用

| 维度         | `_core`                              | `genometracks`                               |
| ------------ | ------------------------------------ | -------------------------------------------- |
| 图形拓扑     | 单函数 → **单** Axes                 | 一个 Figure → **GridSpec 多面板**（每轨一行）|
| 组织方式     | 函数式（`@styled_plot` 装饰器）      | 面向对象（`Track` 类层次，各自 `plot()`）    |
| 建轴         | `get_or_create_axes()` 造**一块** ax | GridSpec 造 N×2 面板矩阵（数据面板+标题面板）|
| 样式落地     | `PlotStyle.apply_to_axes(ax)`        | `PlotStyle.to_track_params()` 下发到每个 track |
| 脊柱策略     | 默认隐藏 top/right（成品带框图）     | 每个轨道自管脊柱（常需全隐藏或按类型定制）   |

关键点：

* **建轴不通用。** `get_or_create_axes()` 只会产出一块简单 Axes，而 genometracks
  需要 GridSpec 面板矩阵，装饰器的自建图分支对它毫无意义。
* **脊柱策略相反。** `_core` 默认隐藏 top/right 对基因组轨道是**错误的**——
  轨道面板往往要隐藏全部脊柱或按类型定制，硬套反而破坏渲染。
* **样式落地方式不同。** 单轴用 `apply_to_axes(ax)`；多轨用 `to_track_params()`
  把颜色/字号/线宽作为"下限"下发到每个 track（用户已显式设置的不覆盖）。

## 5. 适用场景

**用 `_core`（`@styled_plot`）：** 一个函数把一类数据画在一块坐标轴上的**统计图**。
例如 manhattan、qq、venn、admixture、karyoplot。新增此类绘图函数请参照
[基础绘图引擎指南](./base_plotting_engine_guide.zh.md)。

**用 `genometracks`（`plot_tracks`）：** 需要在**共享基因组坐标**上**垂直堆叠多个
轨道**的基因组浏览器式视图（坐标轴、注释、基因模型、覆盖度、比对、lollipop 等）。
新增轨道类型请继承 `Track` 并实现其绘制协议，而非套用 `@styled_plot`。

## 6. 调用链路

### `_core` 单轴路径

```
manhattanplot(data, ..., style="nature", ax=None)
  └─ @styled_plot 包装器
       ├─ use_style("nature")                       # 进入 plotstyle 上下文
       ├─ get_or_create_axes(ax=None, figsize=...)  # 造单块 Axes
       ├─ PlotStyle.apply_to_axes(ax)               # apply_spines=True 时应用脊柱
       └─ 调用函数体在 ax 上作画 → 返回 ax
```

### `genometracks` 多轨路径

```
plot_tracks([GenomeAxisTrack(), AnnotationTrack(...), DataTrack(...)], style="nature")
  ├─ get_style("nature") → resolved_style
  ├─ _apply_style_to_tracks(tracks, style)          # to_track_params() 下发每轨参数
  └─ with use_style(resolved_style):                # 进入同一 plotstyle 上下文
       ├─ _plot_full_layout(...)  → GridSpec 造 N×2 面板矩阵
       └─ 对每个 track 调用其 plot()，各自在自己的 Axes 上作画 → 返回 axes 列表
```

两条链路在**入口 `use_style(...)`** 处交汇于 `plotstyle`，此后分别走
"单轴 `apply_to_axes`" 与 "多轨 `to_track_params` + GridSpec" 两条互不复用的路径。

## 7. 基础设施层：`plotstyle` 与 `palette` 为何各自独立

基础设施层里除了 `plotstyle`，还有一个平行的 `palette` 模块。二者**不合并**，
分别承担不同职责：

* **`palette` = 颜色数据 + 生成工具（值层 / 原语）**：`xkcd_rgb`（命名色库）、
  `circos`、`CYTOBAND_COLORS` / `get_cytoband_color`（基因组学领域常量）、
  `generate_colors_palette`（colormap → 颜色列表工具）。与图形拓扑、与"哪本期刊"
  都无关。
* **`plotstyle` = 出版样式策略（策略层）**：字体、线宽、脊柱可见性、rcParams、
  `PlotStyle` 注册表。调色板只是它众多属性之一，且是**特定策展**的有序短列表
  （Wong / Okabe-Ito / Cell），与 `palette` 的通用色库用途不同。

为何不把 `palette` 并入 `plotstyle`：

* **依赖方向相反。** `plotstyle` **从不导入** `palette`——各 journal 样式把调色板
  内联硬编码为本地常量（如 `_nature.py` 的 `_WONG_PALETTE`）。`palette` 是零依赖
  的叶子原语，并入不会消除任何重复。
* **消费者大多不经过 `plotstyle`。** `palette` 被 `karyotype`、`genometracks`
  （ideogram 染色带）、`baseplot`（venn）、`_core` 直接消费；若 `palette` 变成
  `plotstyle` 的子模块，这些模块仅为取一个颜色常量就得反向依赖整套样式注册机制。
* **抽象层级不同。** 一个是与拓扑/期刊无关的颜色**值**，一个是整体图形**策略**。

```
基础设施层：  plotstyle（样式策略）        palette（颜色原语）
                     ↑                          ↑
       ┌─────────────┴──────────┬───────────────┴─────────────┐
    _core（单轴）          genometracks（多轨）      karyotype / baseplot ...
```

> 可选的反向微调：journal 样式内联的策展色板*可以*收纳进 `palette` 成为命名调色板
> 再由 `plotstyle` 引用，但这只是"常量存放位置"的小重构、收益有限，非必须。

## 8. 参见

* [基础绘图引擎指南（`geneview._core`）](./base_plotting_engine_guide.zh.md)
* 源码：`geneview/_core/`（decorators / canvas / colors）、
  `geneview/genometracks/_track_plot.py`（`plot_tracks`）、
  `geneview/plotstyle/`（`use_style` / `PlotStyle`）
