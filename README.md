# 遥感课程小作业

**Global Distribution of Lunar FeO Abundance and a Case Study of the Chang'e-4 Landing Site**

本项目是一项遥感科学课程的数据分析练习。使用 NASA Planetary Data System（PDS）公开归档的 Lunar Prospector（“月球勘探者号”）伽马射线光谱仪（GRS）元素丰度数据，利用 Python 对月球全球 FeO 丰度进行读取、可视化和面积加权统计，并以嫦娥四号着陆区为案例比较不同空间范围的 FeO 丰度。

**研究范围说明：** 本项目只对 **NASA Lunar Prospector GRS 数据**进行数值处理。嫦娥四号的公开着陆坐标用于地理定位和网格选择；**未下载或分析玉兔二号 VNIS 原始光谱数据**。着陆点对应的 FeO 值是轨道数据网格值，不能当作月球车的原位测量值。

## 1. 数据来源

- **数据平台：** NASA PDS Geosciences Node
- **任务 / 仪器：** Lunar Prospector / Gamma Ray Spectrometer（GRS）
- **数据集编号：** `LP-L-GRS-5-ELEM-ABUNDANCE-V1.0`
- **原始文件：** `lpgrs_high1_elem_abundance_2deg.tab`
- **记录数：** 11,306 条，每条 61 个数值字段
- **主要变量：** `W_FEO`（以 FeO 当量表示的质量分数，g/g）
- **转换：** `FeO (wt.%) = W_FEO × 100`
- **数据说明：** https://pds-geosciences.wustl.edu/missions/lunarp/grs_elem_abundance.html
- **文件下载地址：** https://pds-geosciences.wustl.edu/lunar/lp-l-grs-5-elem-abundance-v1/lp_9001/data/lpgrs_high1_elem_abundance_2deg.tab

原始表格中的网格经纬度边界用于确定空间位置和球面面积权重。
此数据已是 NASA PDS 提供的丰度反演产品，**本项目没有从原始伽马射线计数重新反演 FeO**。

嫦娥四号着陆位置使用 **177.5991°E、45.4446°S**，参考 Liu 等（2019）：https://doi.org/10.1038/s41467-019-12278-3 。

## 2. 运行环境

以下为我在 Windows、Visual Studio Code1.145 和 Anaconda `base` 环境中**实际成功复现**时使用的版本：

| 软件 / 库 | 版本 |
| --- | --- |
| Python | 3.10.9 |
| NumPy | 1.23.5 |
| pandas | 1.5.3 |
| Matplotlib | 3.7.0 |


直接安装：

```bash
python -m pip install numpy==1.23.5 pandas==1.5.3 matplotlib==3.7.0
```

## 3. 项目文件

建议将下面的文件放在**同一个目录**，与本地成功复现时的结构一致：

```text
Lunar-FeO-Distribution-Analysis/
├── data_preprocess.py                  # 数据读取与预处理
├── global_analysis.py                  # 全球分布图及统计
├── ce4_analysis.py                      # 嫦娥四号区域分析
├── README.md                            # 本说明文件
├── processed_feo_grid.csv              # 预处理结果（运行后生成）
├── global_statistics.csv               # 全球统计结果（运行后生成）
├── ce4_comparison_statistics.csv       # 区域统计结果（运行后生成）
├── ce4_neighborhood_10deg_cells.csv    # 10° 区域网格数据（运行后生成）
├── figure1_global_feo_map.png          # 全球 FeO 分布图
├── figure2_area_weighted_histogram.png # 全球面积加权直方图
├── figure3_ce4_regional_feo_map.png    # 嫦娥四号附近 FeO 图
├── figure4_ce4_global_comparison.png   # 不同区域平均值比较图
└── lpgrs_high1_elem_abundance_2deg.tab # NASA 原始数据：自行下载，担心版权不上传 GitHub，见数据来源最后一项
```

**重要：** 此处使用的是作者本地重命名后的脚本名称。

先从上面的官方文件目录下载，确保文件名完全一致，并与三个 `.py` 文件放在同一目录中，再从第 1 步运行。

## 4. 运行顺序

在 VS Code 中打开项目文件夹，在终端依次运行：

### Step 1 — 原始数据预处理

```bash
python data_preprocess.py
```

此脚本：

1. 读取 `lpgrs_high1_elem_abundance_2deg.tab`，检查原始表格的 **11,306 × 61** 数据结构。
2. 提取 `PIXEL_INDEX`、`MIN_LAT`、`MAX_LAT`、`MIN_LON`、`MAX_LON` 和 `W_FEO`。
3. 将 `W_FEO` 从 g/g 转换为 wt.%，同时计算各网格的经纬度中心。
4. 利用嫦娥四号坐标定位包含着陆点的网格，标记为 `CE4_LANDING_PIXEL`。
5. 对 FeO 恰好为零的记录添加 `ZERO_FEO_REVIEW` 标记，**不直接删除**。
6. 输出 `processed_feo_grid.csv`。

核对：着陆点对应 `PIXEL_INDEX = 1716`，`FeO = 7.2173 wt.%`。

### Step 2 — 全球分布与统计

```bash
python global_analysis.py
```

此脚本读取上一步生成的 `processed_feo_grid.csv`，完成：

1. 使用全部 11,306 个地理网格的原始边界绘制 **Figure 1：全球 FeO 分布图**，没有空间插值。
2. 根据每个网格在球面上的实际范围计算面积权重。权重与 `Δλ × [sin(φmax) − sin(φmin)]` 成正比，其中角度须换成弧度。
3. 计算全球面积加权平均值，以及按网格记录计算的算术平均值、中位数、百分位数等。
4. 绘制 **Figure 2：面积加权 FeO 丰度直方图**，并保存 `global_statistics.csv`。
5. 对唯一 FeO 零值记录进行保留 / 排除的全球平均值敏感性检查。

**说明：** Figure 1 是按经纬度画出的矩形地图，并非真正的等面积地图投影。现有 Figure 2 代码的直方图区间截至 30 wt.%；原始最大值为 30.749 wt.%，因此最末端高值未包含在该直方图的区间内，但仍参与全球平均值等统计计算。若用于严格的面积比例展示，应调整直方图区间覆盖全部记录。

### Step 3 — 嫦娥四号着陆区分析

```bash
python ce4_analysis.py
```

此脚本：

1. 使用嫦娥四号着陆点位置 `177.5991°E, 45.4446°S`，在 PDS 网格中定位包含该点的网格。
2. 计算着陆点与各网格**中心**之间的球面大圆角距离。
3. 分别筛选中心位于 **5°、10°、15°** 范围内的网格，并计算面积加权 FeO 平均值。
4. 生成 **Figure 3：嫦娥四号附近 FeO 分布图** 和 **Figure 4：全球与不同区域的 FeO 对比图**。
5. 输出 `ce4_comparison_statistics.csv` 和 `ce4_neighborhood_10deg_cells.csv`。

区域筛选按**网格中心**进行；落入范围的网格全部计入面积权重，不对边界网格做局部裁剪。5°、10°、15° 研究范围彼此重叠，因此不是独立样本；其中 10° 对应的月表弧长约为 303 km，远大于月球车行驶尺度。所以15° 只是做个补充，基本没啥用。

## 5. 复现结果核对

以下结果已使用 Python 3.10.9 在作者本地电脑上运行验证：

| 分析对象 | 网格数 | FeO（wt.%） |
| --- | ---: | ---: |
| 全球面积加权平均 | 11,306 | 7.405905 |
| 嫦娥四号着陆点所在网格 | 1 | 7.217300 |
| 着陆点周围 5° | 23 | 9.175638 |
| 着陆点周围 10° | 90 | 8.892866 |
| 着陆点周围 15° | 200 | 9.060852 |

其他全球统计指标：

- 网格算术平均值：**7.378923 wt.%**
- 网格中位数：**5.498156 wt.%**
- 网格第 5 / 95 百分位数：**2.953825 / 19.119000 wt.%**
- 原始 FeO 最小 / 最大值：**0 / 30.749 wt.%**
- FeO 为零的网格：**1 个**，已标记，默认保留在统计中


## 6. 研究局限与使用说明

- 区域平均值取决于网格选择方法、范围大小及数据反演误差；本项目未进行完整的误差传播分析。
- 本项目主要展示公开遥感数据的获取、清洗、绘图和统计流程，不对复杂月球地质过程作超出数据支持范围的推断。

## 7. 数据引用

1. NASA Planetary Data System, Lunar Prospector GRS Elemental Abundance Dataset, `LP-L-GRS-5-ELEM-ABUNDANCE-V1.0`. https://pds-geosciences.wustl.edu/missions/lunarp/grs_elem_abundance.html
2. Liu et al. (2019), *Nature Communications*, 嫦娥四号着陆位置研究。https://doi.org/10.1038/s41467-019-12278-3

本 README 为项目复现说明。原始 NASA 数据的来源需要下载。
