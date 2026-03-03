# ClusterMap V9 Workflow Documentation

## 概述

这是一个基于 tile 的空间转录组学分段聚类流程，使用分布式计算架构（SLURM）处理大规模空间转录组数据。

## 流程架构

### 三个阶段

```
┌─────────────────────────────────────────────────────────────────┐
│                         PREP 阶段 (共享)                         │
│  - DAPI 图像加载与旋转 (270°)                                      │
│  - 特征点坐标转换                                                 │
│  - 全局 spot 检测与预处理                                          │
│  - 图像分块 (tiling)                                              │
│  - 输出：manifest.json, out.pkl, global_spots.pkl, tiles/       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ├──────────────────────────────────┐
                              │                                  │
                              v                                  v
                    ┌──────────────────┐              ┌──────────────────┐
                    │  TAG1 (参数组合 1) │              │  TAG2 (参数组合 2) │
                    │   DRIVER → TILE   │              │   DRIVER → TILE   │
                    │   → RESCUE→STITCH │              │   → RESCUE→STITCH │
                    └──────────────────┘              └──────────────────┘
```

### 1. PREP 阶段（预处理）

**目的**: 全局预处理，为多个参数组合共享

**资源**: 60 CPUs, 480G 内存，8 小时

**关键步骤**:
1. 加载 DAPI 图像 (ch0X.tif)
2. 旋转 270°（适应坐标系）
3. 读取 goodPoints_max3d.csv（特征点）
4. 坐标转换（旋转校正）
5. 初始化 ClusterMap 模型
6. 运行 preprocess()（网格化、过滤）
7. 将图像分割为重叠 tiles（window_size=300, overlap=20%）
8. 保存全局 spots 和 per-tile 输入数据

**输出文件**:
```
WORK_ROOT/{POS}/_shared_prep/
├── manifest.json              # 元数据（n_tiles, window_size, rotation 等）
├── out.pkl                    # DataFrame[tile_id, img, spots, label_img]
├── global_spots.pkl           # 全局 spot 表（含 index 列）
├── max_rotated_dapi.tif       # DAPI 最大投影（用于 FIJI 验证）
└── tiles/
    ├── tile_00000.npz         # DAPI + label_img
    └── tile_00000_spots.csv   # spots（含全局 index）
```

### 2. TILE 阶段（并行分段）

**目的**: 对每个 tile 独立执行 ClusterMap 分段

**资源**: 8 CPUs, 64G 内存，24 小时/tile，60 个 tile 并行

**关键步骤**:
1. 加载 tile 的 DAPI 和 spots
2. 检查 spots 数量 >= reads_filter (5)
3. 初始化 ClusterMap（tile 级别）
4. 运行 segmentation()：
   - DAPI 网格化
   - 细胞核检测
   - 基因表达分配
   - clustermap 标签生成
5. 保存 tile 模型

**输出**: `tile_results/tile_00000.pkl` (ClusterMap 模型)

### 3. RESCUE 阶段（故障恢复）

**目的**: 重试失败/缺失的 tiles

**逻辑**:
- 检查哪些 tile_XXXXX.pkl 缺失或损坏
- 提交 RETRY 数组作业（最多 20 轮）
- 每轮重试失败的 tiles
- 全部成功或达到最大轮次后退出

**并发控制**: RETRY_CONCURRENCY=12（每轮 12 个并行）

### 4. STITCH 阶段（拼接）

**目的**: 将 tile 级别的聚类结果合并到全局坐标系

**资源**: 16 CPUs, 192G 内存，12 小时

**关键步骤**:
1. 加载 global_spots.pkl（全局坐标）
2. 初始化累加器模型（accumulator model）
3. 按顺序加载每个成功的 tile 模型
4. 调用 `model_acc.stitch(tile_model, out, tid)`
5. 通过全局 ID 写回 clustermap 分配
6. 输出最终结果

**最终输出**:
```
OUT_ROOT/{POS}/{TAG}_IC{IC}_ID{ID}_.../
├── model_final.pkl           # 完整 ClusterMap 模型
├── remain_reads_raw.csv      # 所有 spots（含 clustermap 列）
├── remain_reads.csv          # 过滤后（clustermap >= 0）
├── cell_center.csv           # 细胞中心坐标
└── .cm_done                  # 完成标记
```

## 目录结构

### 原始数据
```
BASE/01_data/{ROUND_DAPI}/{POS}/
├── ch00.tif, ch01.tif, ...   # DAPI 通道
└── genes.csv                  # 基因列表

BASE/02_registration/{POS}/
└── goodPoints_max3d.csv       # 特征点坐标
```

### 工作目录
```
WORK_ROOT/{POS}/
├── _shared_prep/              # PREP 输出（共享）
│   ├── manifest.json
│   ├── out.pkl
│   ├── global_spots.pkl
│   └── tiles/
│
└── {TAG}/                     # 参数组合特定
    ├── manifest.json -> symlink
    ├── out.pkl -> symlink
    ├── global_spots.pkl -> symlink
    ├── tiles -> symlink
    └── tile_results/
        └── tile_XXXXX.pkl
```

### 输出目录
```
OUT_ROOT/{POS}/{TAG}_IC{IC}_ID{ID}_IPF{IPF}_R{XY}_{Z}/
├── model_final.pkl
├── remain_reads_raw.csv
├── remain_reads.csv
├── cell_center.csv
└── .cm_done
```

## SLURM 作业流

```
cm_prep_{POS} (1 job, 60CPU/480G/8h)
    │
    └─afterok→ cm_drive_{POS}_{TAG} (1 job, 1CPU/2G/30m)
                 │
                 ├─→ cm_tile_{POS}_{TAG} (array: 0..N-1, 8CPU/64G/24h, 60 并发)
                 │      │
                 │      └─afterany→ cm_rescue_{POS}_{TAG} (检查缺失，提交 RETRY)
                 │                   │
                 │                   └─→ cm_retry_{POS}_{TAG}_{ROUND} (array, 最多 20 轮)
                 │
                 └─afterok→ cm_stitch_{POS}_{TAG} (1 job, 16CPU/192G/12h)
```

### 作业依赖
- **DRIVER** → 依赖 PREP 完成 (afterok)
- **TILES** → 由 DRIVER 提交
- **RESCUE** → 依赖 TILES 完成 (afterany)
- **STITCH** → 依赖 RESCUE 完成 (afterok)

## 参数文件

### params-ky8.tsv 格式
```
TAG     IC      ID      IPF     ICR_xy  ICR_z
AA1     0.080   4       0.003   50      10
AA2     0.100   4       0.003   50      10
...
```

**参数说明**:
- **TAG**: 参数组合标识符
- **IC** (input_cell_num_threshold): 细胞最小 spots 阈值（0.08 = 8%）
- **ID** (input_dapi_grid_interval): DAPI 网格间隔（4 像素）
- **IPF** (input_pct_filter): 过滤百分比（0.003 = 0.3%）
- **ICR_xy** (xy_radius): XY 方向细胞半径（50 像素）
- **ICR_z** (z_radius): Z 方向细胞半径（10 像素）

## 关键设计模式

### 1. 全局 Spot ID 管理
```python
# 所有阶段维护统一的 index 列
global_spots.index = global_spots['index'] = [0, 1, 2, ...]
# PREP → TILE → STITCH 全程保持一致
```

### 2. 共享 PREP
- 一个 Position 只运行一次 PREP
- 多个 TAG 共享同一个 PREP 输出
- 通过符号链接避免数据冗余

### 3. 原子写入
```python
# 所有关键文件使用临时文件 + mv
_atomic_joblib_dump(obj, path, compress=3)
# 写入 path.tmp → mv path.tmp path
```

### 4. 幂等性检查
```python
is_done_strict(out_dir, work_dir):
    - model_final.pkl 存在且非空
    - remain_reads*.csv 存在且非空
    - model_final.pkl 不旧于最新的 tile 结果
```

### 5. 错误恢复
- TILE 失败 → RESCUE 检测 → RETRY 重试
- 最多 20 轮重试
- 指数退避（2s, 4s, 8s, ...）

## Python 入口

### 主脚本
```
new_clustermaptest_dist_ky_fast_new.py
```

### 使用方式
```bash
# PREP 阶段
python new_clustermaptest_dist_ky_fast_new.py \
    --stage prep \
    -IP Position1 -IZ 42 -IC 0.08 -ID 4 -ICR "50,10" -IXY 2048 \
    -IDR 1 -IPF 0.003 -IEP T -IDir /path/to/data \
    -Igood_points_max3d goodPoints.csv \
    -OP /path/to/output \
    --work_dir /path/to/work \
    --window_size 300 --overlap_percent 0.2

# TILE 阶段
python new_clustermaptest_dist_ky_fast_new.py \
    --stage tile \
    --tile_id 0 \
    [其他参数同上] \
    --work_dir /path/to/work

# STITCH 阶段
python new_clustermaptest_dist_ky_fast_new.py \
    --stage stitch \
    [其他参数同上] \
    --work_dir /path/to/work \
    -OP /path/to/output
```

### ClusterMap 模块
```
ClusterMap/
├── clustermap.py          # 核心类：ClusterMap
├── preprocessing.py       # 预处理：网格化、过滤
├── stitch.py              # 拼接逻辑
├── utils.py               # 工具：get_img(), split()
├── Points2Cell.py         # 细胞调用
├── metrics.py             # 质量指标
└── postprocessing.py      # 输出格式化
```

## 输出文件格式

### remain_reads_raw.csv
```csv
index,gene_name,spot_location_1,spot_location_2,spot_location_3,gene,clustermap,cell_center_0,cell_center_1,cell_center_2,is_noise
0,ACTB,1024,2048,15,1,0,2048,1024,15,0
1,GAPDH,1030,2055,16,2,0,2048,1024,15,0
2,TUBB,500,1500,20,3,-1,-1,-1,-1,-1  # 未分配
```

### cell_center.csv
```csv
cell_barcode,column,row,z_axis
0,1024,2048,15
1,1030,2055,16
2,500,1500,20
```

### manifest.json
```json
{
  "n_tiles": 42,
  "window_size": 300,
  "overlap_percent": 0.2,
  "margin": 60,
  "img_z": 42,
  "left_rotation": 270,
  "reads_filter": 5,
  "useown_dapi_bi": false,
  "filter_value2": 1,
  "xy_radius": 50,
  "z_radius": 10,
  "cell_num_threshold": 0.08,
  "dapi_grid_interval": 4,
  "pct_filter": 0.003,
  "global_id_col": "index"
}
```

## 资源估算

### 单个 Position（单个 TAG）
| 阶段 | CPU 时 | 内存 | 时间 | 磁盘 |
|------|-------|------|------|------|
| PREP | 480 | 480G | 8h | ~50GB |
| TILES | 480 (60×8) | 64G×60 | 24h | ~20GB |
| RESCUE | ~60 | 480G | 24h | - |
| STITCH | 192 | 192G | 12h | ~10GB |
| **总计** | **1212** | - | **~2 天** | **~80GB** |

### 并发建议
- 单个 Position：1 个 pipeline
- 多个 Position：并行运行（SHARE_PREP_ACROSS_PARAMS=1 共享 PREP）
- MAX_LIVE_PIPELINES=0（无限制）或设置上限

## 日志与监控

### 日志文件
```
logs_cm_V9/
├── cm_prep_Position1.12345.out
├── cm_drive_Position1_AA1.12346.out
├── cm_tile_Position1_AA1.12347.out
├── cm_rescue_Position1_AA1.12348.out
└── cm_stitch_Position1_AA1.12349.out
```

### 提交记录
```
submitted.tsv
POS    TAG    params    prep_jobid    driver_jobid    tiles_jobid    rescue_jobid    stitch_jobid
Position1    AA1    IC0.08_ID4_IPF0.003_R50_10    -    12346    12347    12348    12349
```

## 常见问题

### Q: PREP 失败如何处理？
A: 检查 DAPI 文件路径和 goodPoints 文件是否存在。确认内存足够（480G）。

### Q: TILE 大量失败？
A: 检查 tile 输入文件（tiles/）是否完整。增加单 tile 内存（默认 64G）。

### Q: STITCH 内存不足？
A: 增加 STITCH 作业内存（默认 192G）。或减少并发 pipeline 数量。

### Q: 如何重用 PREP 结果？
A: 设置 `SHARE_PREP_ACROSS_PARAMS=1`（默认）。新 TAG 会自动 symlink 已有 PREP 输出。

## 相关文件

- `run_cm_sweep_multiPos-V9.sh`: 主调度脚本
- `new_clustermaptest_dist_ky_fast_new.py`: Python 入口
- `ClusterMap/`: 核心模块
- `params-ky8.tsv`: 参数组合定义
