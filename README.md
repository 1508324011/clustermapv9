# ClusterMap V9 快速入门指南

## 文件结构

```
clystermapv9/
├── README.md                              # 本文件
├── WORKFLOW.md                            # 详细流程文档
├── run_cm_sweep_multiPos-V9.sh            # 主调度脚本（多 Position 批量）
├── params-ky8.tsv                         # 参数组合定义文件
├── new_clustermaptest_dist_ky_fast_new.py # Python 入口脚本
└── ClusterMap/                            # ClusterMap 核心模块
```

## 快速开始

### 1. 准备数据

确保以下数据存在：

```bash
# DAPI 图像
/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/01_data/round001/Position1/ch00.tif
/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/01_data/round001/Position1/ch01.tif
...

# 基因列表
/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/01_data/genes.csv

# 特征点（registration 输出）
/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/02_registration/Position1/goodPoints_max3d.csv
```

### 2. 修改配置

编辑 `run_cm_sweep_multiPos-V9.sh`：

```bash
# 第 23 行起：基础路径
BASE="/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica"
WORK_ROOT="/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/03_segmentation_ky/clustermap_ky_tmp"
OUT_ROOT="/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/03_segmentation_ky/clustermap_ky"

# 第 31 行起：Position 列表
POSITIONS="Position1"  # 单个 Position
# 或
POSITIONS="Position1,Position109,Position217"  # 多个
# 或
POSITIONS="Position1-108"  # 范围

# 第 36 行起：参数文件
PARAM_TSV="./params-ky8.tsv"
```

### 3. 提交作业

```bash
# 进入目录
cd /gpfs/share/home/2201112433/zhui/clystermapv9

#  dry-run（预览，不实际提交）
bash run_cm_sweep_multiPos-V9.sh --dry-run 1

# 正式提交
bash run_cm_sweep_multiPos-V9.sh
```

### 4. 监控进度

```bash
# 查看提交的作业
cat logs_cm_V9/submitted.tsv

# 查看 SLURM 队列
squeue -u $USER

# 查看日志
tail -f logs_cm_V9/cm_stitch_Position1_AA1.*.out
```

### 5. 检查结果

```bash
# 完成标记
ls -la /gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/03_segmentation_ky/clustermap_ky/Position1/AA1_IC0.080_ID4_IPF0.003_R50_10/.cm_done

# 查看输出
ls -la /gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/03_segmentation_ky/clustermap_ky/Position1/AA1_*/
# - model_final.pkl
# - remain_reads.csv
# - cell_center.csv
```

## 单个 Position + 单个 TAG 示例

如果要手动运行单个 Position 和单个参数组合（不通过批量脚本）：

### 步骤 1: PREP（共享预处理）

```bash
POS="Position1"
TAG="AA1"
WORK_DIR="/path/to/work/${POS}/_shared_prep"
OUT_DIR="/path/to/output/${POS}"

python new_clustermaptest_dist_ky_fast_new.py \
    --stage prep \
    -IP "${POS}" \
    -IZ 42 \
    -IC 0.08 \
    -ID 4 \
    -ICR "50,10" \
    -IXY 2048 \
    -IDR 1 \
    -IPF 0.003 \
    -IEP T \
    -IDir "${BASE}/01_data/round001/${POS}" \
    -Igood_points_max3d "${BASE}/02_registration/${POS}/goodPoints_max3d.csv" \
    -OP "${OUT_DIR}" \
    --work_dir "${WORK_DIR}" \
    --window_size 300 \
    --overlap_percent 0.2
```

### 步骤 2: TILE（并行分段）

```bash
# 获取 tile 数量
N_TILES=$(cat ${WORK_DIR}/manifest.json | python -c "import sys,json; print(json.load(sys.stdin)['n_tiles'])")

# 提交数组作业（0 到 N_TILES-1）
for TID in $(seq 0 $((N_TILES-1))); do
    python new_clustermaptest_dist_ky_fast_new.py \
        --stage tile \
        --tile_id ${TID} \
        [其他参数同上] \
        --work_dir "${WORK_DIR}"
done
```

### 步骤 3: STITCH（拼接）

```bash
python new_clustermaptest_dist_ky_fast_new.py \
    --stage stitch \
    -IP "${POS}" \
    [其他参数同上] \
    --work_dir "${WORK_DIR}" \
    -OP "${OUT_DIR}/${TAG}_IC0.080_ID4_IPF0.003_R50_10"
```

## 参数调优

### 敏感参数（推荐测试多个组合）

| 参数 | 含义 | 推荐范围 | 影响 |
|------|------|----------|------|
| IC (cell_num_threshold) | 细胞最小 spots 比例 | 0.05 - 0.15 | 越小→细胞越多，可能过分割 |
| ID (dapi_grid_interval) | DAPI 网格间隔 | 3 - 6 | 越小→精度越高，速度慢 |
| IPF (pct_filter) | 噪声过滤比例 | 0.001 - 0.01 | 越小→保留更多 spots |
| ICR_xy (xy_radius) | 细胞 XY 半径 | 40 - 60 | 需匹配实际细胞大小 |
| ICR_z (z_radius) | 细胞 Z 半径 | 8 - 15 | 需匹配实际细胞大小 |

### 测试策略

1. **粗调**: 用 `params-ky8.tsv` 中的 5-10 个组合测试单个 Position
2. **评估**: 检查细胞数量、spots 分配率、细胞大小分布
3. **精选**: 选择 2-3 个最佳参数用于全部 Position

## 故障排查

### PREP 失败

```bash
# 检查 DAPI 文件
ls -lh ${BASE}/01_data/round001/${POS}/ch0*.tif

# 检查特征点
head ${BASE}/02_registration/${POS}/goodPoints_max3d.csv

# 查看日志
cat logs_cm_V9/cm_prep_${POS}.*.out
```

### TILE 大量失败

```bash
# 检查 PREP 输出
ls -lh ${WORK_ROOT}/${POS}/_shared_prep/tiles/

# 检查 manifest
cat ${WORK_ROOT}/${POS}/_shared_prep/manifest.json

# 查看失败 tile 日志
grep -l "FAILED" logs_cm_V9/cm_tile_*.out
```

### STITCH 卡住

```bash
# 检查 tile 完成情况
ls ${WORK_ROOT}/${POS}/${TAG}/tile_results/tile_*.pkl | wc -l

# 对比 manifest 的 n_tiles
cat ${WORK_ROOT}/${POS}/${TAG}/manifest.json | python -c "import sys,json; print(json.load(sys.stdin)['n_tiles'])"

# 内存不足？增加 STITCH 资源
# 编辑脚本：STITCH_CPUS=16, STITCH_MEM=192G, STITCH_TIME=12:00:00
```

## 性能优化

### 加速 PREP

```bash
# 减少 CPU（边际效应）
PREP_CPUS=48  # 默认 60

# 降低 DAPI 旋转精度（如果图像质量允许）
export CM_ROTATE_ORDER=0  # 默认 1
```

### 加速 TILE

```bash
# 增加并发
ARRAY_CONCURRENCY=80  # 默认 60

# 减少单 tile 时间（如果数据简单）
TILE_TIME=12:00:00  # 默认 24:00:00
```

### 节省磁盘

```bash
# 清理中间文件（STITCH 完成后）
rm -rf ${WORK_ROOT}/${POS}/${TAG}/tile_results/
rm ${WORK_ROOT}/${POS}/${TAG}/*.pkl
# 保留：manifest.json, global_spots.pkl（可重用）
```

## 输出解释

### remain_reads.csv

```csv
index,gene_name,spot_location_1,spot_location_2,spot_location_3,gene,clustermap,cell_center_0,cell_center_1,cell_center_2,is_noise
0,Actb,1024,2048,15,1,0,2048,1024,15,0
# clustermap=0 → 分配到细胞 0
# cell_center_* → 细胞 0 的中心坐标
```

### cell_center.csv

```csv
cell_barcode,column,row,z_axis
0,1024,2048,15
1,1030,2055,16
# 每个细胞的 ID 和中心坐标
```

### model_final.pkl

```python
import joblib
model = joblib.load('model_final.pkl')

# 查看细胞数量
len(model.cellid_unique)

# 查看 spots 分配
model.spots[model.spots['clustermap'] >= 0].shape[0]
```

## 进阶用法

### 多轮 RESCUE 配置

```bash
# 编辑脚本（约第 280 行）
MAX_RESCUE_ROUNDS=20  # 默认 20 轮
RETRY_CONCURRENCY=12  # 每轮 12 个并行
```

### 自定义并发限制

```bash
# 限制同时运行的 pipeline 数量
MAX_LIVE_PIPELINES=3  # 默认 0（无限制）

# 适用于内存紧张的情况
```

### 使用特定 Python 环境

```bash
# 编辑脚本（约第 18 行）
PY="/path/to/your/python"

# 确保环境包含：
# - numpy, pandas, scipy, scikit-learn
# - tifffile, joblib
# - ClusterMap 模块在 PYTHONPATH
```

## 相关文件

- `WORKFLOW.md`: 完整流程文档（架构、依赖、文件格式）
- `run_cm_sweep_multiPos-V9.sh`: 批量调度脚本
- `params-ky8.tsv`: 参数组合示例
- `new_clustermaptest_dist_ky_fast_new.py`: Python 主入口
- `ClusterMap/`: 核心算法模块

## 获取帮助

1. 查看 `WORKFLOW.md` 了解详细架构
2. 检查日志文件（`logs_cm_V9/*.out`）
3. 查看 `submitted.tsv` 追踪作业状态
4. 使用 `--dry-run 1` 预览作业提交流程
