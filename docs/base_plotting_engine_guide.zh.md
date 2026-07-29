# Geneview 基础绘图引擎（`geneview._core`）开发者指南

> 读者对象：希望在共享的基础绘图引擎之上，为 geneview 新增绘图函数，
> 或重构现有绘图函数的贡献者。

## 1. 为什么需要引擎？

在历史代码中，geneview 的每一个顶层绘图函数（`manhattanplot`、`qqplot`、
`admixtureplot`、`venn`、`karyoplot` 等）都重复实现了同一套样板逻辑：

```python
def someplot(data, ..., ax=None, style=None, **kwargs):
    # 1. 校验输入
    with use_style(style):                 # 进入样式上下文
        if ax is None:                     # 创建画布
            _, ax = subplots(figsize=(...))
        ...                                # 绘图
        ax.spines["top"].set_visible(False)  # 脊柱处理（有的函数有，有的没有）
        ax.set_title(...); ax.set_xlabel(...)
    return ax
```

这套“进入样式 → 获取坐标轴 → 处理脊柱”的生命周期被（略有差异地）复制到
每个模块中，久而久之出现了漂移：有的函数隐藏 top/right 脊柱，有的没有；
有的用 `plt.gca()`，有的用 `subplots(...)`；`karyoplot` 甚至根本不支持
`style=` 参数。

`geneview._core` 把这套生命周期收敛为**一个经过测试的引擎**，从而实现：

* 所有绘图函数行为一致；
* 新增绘图函数只需写几行绘图代码；
* 引擎级的改进（例如新的默认行为）一次落地、处处生效。

该引擎是**纯增量式**的——所有被迁移函数的对外 API 保持不变
（`karyoplot` 只是**新增**了一个可选的 `style=` 关键字参数）。

## 2. 在整体架构中的位置

```
┌─────────────────────────────────────────────────────────────┐
│  绘图模块层：gwas / popgene / baseplot / karyotype / ...       │  ← 绘图
├─────────────────────────────────────────────────────────────┤
│  geneview._core（styled_plot、get_or_create_axes、colors）     │  ← 引擎
├─────────────────────────────────────────────────────────────┤
│  基础设施层：plotstyle / palette / utils                       │  ← 共享
└─────────────────────────────────────────────────────────────┘
```

`_core` **只**依赖 `plotstyle` 与 `palette`，从不导入任何绘图模块，因此不存在
循环导入。绘图模块通过 `from .._core import ...` 使用引擎。

> **`_core` 与 `genometracks` 的关系：** 本引擎服务于"单函数→单 Axes"的统计图；
> `geneview.genometracks`（多轨道基因组浏览器）是另一套并列入口，它
> **不复用** `_core`，但与 `_core` **共享同一个 `plotstyle` 样式底座**。
> 两者为何这样分工、各自的适用场景与调用链路，见
> [绘图架构设计：`_core` 与 `genometracks`](./plotting_architecture_design.zh.md)。

## 3. 对外可用的构件

从包根统一导入：

```python
from geneview._core import styled_plot, get_or_create_axes, color_cycle, resolve_colors
```

### 3.1 `@styled_plot(...)` —— 装饰器（从这里开始）

它为绘图函数包裹上共享的“样式 + 画布”生命周期。在被装饰的函数体内，你可以
认为 **`ax` 一定是一个真实的 `Axes`**，且**样式上下文已经激活**——你只需要
编写绘图代码。

```python
@styled_plot(figsize=(9, 3), subplot_kws={"facecolor": "w", "edgecolor": "k"})
def manhattanplot(data, ..., ax=None, style=None, **kwargs):
    ax.scatter(...)          # ax 一定非 None
    ax.set_xlabel(...)
    return ax                # 始终返回 Axes
```

装饰器在每次调用时依次完成：

1. 读取 `style` 关键字并进入 `use_style(style)`（当 `style is None` 时为空操作，
   此时仍沿用当前全局激活的样式）；
2. 读取 `ax` 关键字；当其为 `None` 时，用 `figsize` / `subplot_kws` 通过
   `get_or_create_axes` 创建图形/坐标轴；
3. 当 `apply_spines=True`（默认）时，把激活样式的脊柱可见性规则应用到坐标轴上
   （`PlotStyle.apply_to_axes`）；
4. 调用你的函数体；
5. 如果函数体抛出异常**且**图形是由装饰器创建的，则关闭该图形，避免非法输入
   遗留半成品图形。

对外签名与 docstring 通过 `functools.wraps` 得到保留，因此调用方与 Sphinx
文档都感知不到差异。

**装饰器参数**

| 参数           | 默认值  | 含义 |
| -------------- | ------- | ---- |
| `figsize`      | `None`  | 函数自建图形时的默认 `(宽, 高)`；`None` 时回退到激活样式的 `figure_figsize`。 |
| `apply_spines` | `True`  | 是否强制应用样式的脊柱规则。对自行管理边框的图（条形图、关闭坐标轴的图）设为 `False`。 |
| `use_gca`      | `False` | 为 `True` 时，省略 `ax` 将复用 `plt.gca()` 而非新建图形（`karyoplot` 的历史行为）。 |
| `subplot_kws`  | `None`  | 传给 `plt.subplots` 的额外参数（`facecolor`、`constrained_layout` 等）。 |

**要求：** 被装饰函数**必须**暴露 `ax` 关键字参数（否则装饰器抛出 `TypeError`）。
若存在 `style` 关键字则会被识别使用；没有该关键字的函数则始终沿用激活样式。

### 3.2 `get_or_create_axes(...)` —— 画布助手

用于替代 `if ax is None: _, ax = subplots(...)` 的可复用函数。当函数需要比
装饰器更精细的控制时可直接调用（例如一个既被装饰器调用、又需要独立运行的
**内部**绘图助手——参见 `_draw_admixtureplot`）。

```python
def get_or_create_axes(ax=None, *, figsize=None, style=None,
                       apply_spines=False, use_gca=False, **subplot_kws) -> Axes
```

* 传入了 `ax` → 原样返回（可选地应用脊柱规则）；
* `ax is None` 且 `use_gca=False` → `plt.subplots(figsize=..., **subplot_kws)`；
* `ax is None` 且 `use_gca=True` → `plt.gca()`；
* `figsize is None` → 回退到解析出的样式的 `figure_figsize`；
* `style is None` → 使用当前激活样式（因此在 `use_style(...)` 上下文中调用时——
  例如在 `@styled_plot` 内部——结果是正确的）。

注意其默认值与装饰器不同：这里 `apply_spines` 默认 `False`（底层助手保持中性），
而装饰器默认 `True`（成品图的常见需求）。

### 3.3 `color_cycle(color)` —— 颜色循环

返回一个无限迭代器，完全复现历史上 manhattan/qq 的规则：

| 输入                       | 结果 |
| -------------------------- | ---- |
| `"#3B5488,#53BBD5"`        | 循环 `["#3B5488", "#53BBD5"]` |
| `"rb"`（无逗号）           | 循环**字符** `"r", "b", ...` |
| `["r", "g", "b"]`          | 循环列表元素 |

```python
colors = color_cycle(color)
for group in groups:
    c = next(colors)
```

### 3.4 `resolve_colors(palette, n_colors, alpha=1.0)` —— 调色板 → 颜色

对 `geneview.palette.generate_colors_palette` 的轻量共享封装。接受色图名称、
显式颜色列表或 `Colormap`，返回颜色列表（若调色板无法提供足够颜色，长度可能
小于 `n_colors`——请像 `admixtureplot` 那样检查长度并给出警告）。

## 4. 编写一个新绘图函数（配方）

```python
# geneview/<subpackage>/_myplot.py
import numpy as np
from .._core import styled_plot, color_cycle   # 按需再加 get_or_create_axes / resolve_colors


@styled_plot(figsize=(6, 4), subplot_kws={"facecolor": "w"})
def myplot(data, ax=None, color="#3B5488,#53BBD5",
           title=None, xlabel=None, ylabel=None, style=None, **kwargs):
    """一句话概述。

    Parameters
    ----------
    ...
    style : str、PlotStyle 或 None，可选
        要应用的绘图样式：已注册样式名（"nature"、"science"、"cell"）、
        PlotStyle 对象，或 None（使用当前激活样式）。
    ax : matplotlib axis，可选
        绘图所用坐标轴；省略时自动创建。

    Returns
    -------
    ax : matplotlib Axes
    """
    # 1. 校验输入（尽早抛出、给出清晰信息）。
    #    提示：轻量校验可放在绘图之前；若在此抛异常，装饰器会关闭自建的图形。

    # 2. 绘图。此时 `ax` 一定非 None，样式上下文已激活。
    colors = color_cycle(color)
    ax.plot(...)

    # 3. 标题 / 标签。
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax
```

随后像现有函数一样，在子包的 `__init__.py` 和顶层 `geneview/__init__.py` 中导出它。

### 约定

* **始终**接受 `ax=None` 与 `style=None`，并**始终**返回 `Axes`。
* 保持把 `**kwargs` 透传给底层的 matplotlib 调用（如 `scatter`）。
* 领域逻辑（染色体偏移、λ 膨胀系数、聚类、Venn 几何等）放在模块内——引擎刻意
  保持通用，对基因组学一无所知。
* 在 `geneview/tests/test_<name>.py` 下补充测试。

## 5. `apply_spines` 取值指南

装饰时唯一真正需要决策的，是“脊柱是否交由样式管理”。被迁移函数已确立范式：

| 函数            | `apply_spines` | 原因 |
| --------------- | -------------- | ---- |
| `manhattanplot` | `True`         | 标准带框图；geneview 默认隐藏 top/right（与旧的手写脊柱代码一致）。 |
| `qqplot` / `qqnorm` | `True`     | 同样是带框图契约；现在与激活样式保持一致。 |
| `admixtureplot` | `False`        | 自绘边框：对**四条**脊柱设置了自定义线宽。 |
| `venn`（`vennx`） | `False`      | 调用了 `ax.set_axis_off()`，脊柱可见性无意义。 |
| `karyoplot`     | `False` + `use_gca=True` | 保留 `plt.gca()` 历史回退与自定义坐标轴；不希望隐藏 top/right。 |

经验法则：标准 x/y 带框图用 **`True`**；自绘边框或关闭坐标轴时用 **`False`**。

## 6. 需要独立运行的内部绘图助手

部分模块把“公开入口”与“内部 `_draw_*` 工作函数”分离，后者会被单元测试直接调用
（如 `_draw_admixtureplot`）。装饰器只包裹**公开**函数，因此该工作函数在被单独
调用时仍需自行获取坐标轴。把这一步交给 `get_or_create_axes`，使用相同的默认值，
这样当装饰器已提供坐标轴时它就是一个空操作：

```python
def _draw_admixtureplot(..., ax=None):
    ax = get_or_create_axes(
        ax, figsize=(14, 2), apply_spines=False, 
        facecolor="w", constrained_layout=True
    )
    ...
```

## 7. 行为与兼容性说明

* **API 未变。** 所有被迁移函数的签名与返回值保持不变；`karyoplot` 只是**新增**了
  可选的 `style=` 关键字。
* **`qqplot`/`qqnorm` 的脊柱：** 此前它们保留全部四条脊柱（matplotlib 默认）；
  现在改为遵循激活样式、隐藏 top/right，与 `manhattanplot` 及期刊样式保持一致。
  这是有意为之的一致化行为，并非 API 变更。
* **图形尺寸保持不变。** 每个被装饰函数都传入其历史默认 `figsize`，因此独立
  （不传 `ax`）时的输出在视觉上与之前一致（manhattan `9x3`、qq `5x5`、
  admixture `14x2`、venn `7x7`）。
* **不泄漏图形。** 若装饰器创建图形后校验才抛异常，装饰器会关闭该图形。

## 8. 测试

运行被迁移模块的测试以及完整测试套件：

```bash
python -m pytest geneview/tests/test_manhattan.py geneview/tests/test_qq.py \
    geneview/tests/test_venn.py geneview/tests/test_admixture.py \
    geneview/tests/test_karyotype.py -q
python -m pytest geneview/tests/ -q
```

## 9. 文件结构

```
geneview/_core/
├── __init__.py       # 对外导出
├── canvas.py         # get_or_create_axes
├── colors.py         # color_cycle、resolve_colors
└── decorators.py     # styled_plot
```
