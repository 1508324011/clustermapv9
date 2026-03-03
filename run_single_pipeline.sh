#!/bin/bash
# ClusterMap V9 - 单 Position 单 TAG 运行脚本
# 用法：bash run_single_pipeline.sh Position1 AA1

set -e

# ==================== 配置 ====================

# 基础路径
BASE="/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica"
WORK_ROOT="/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/03_segmentation_ky/clustermap_ky_tmp"
OUT_ROOT="/gpfs/share/home/2501111695/labShare/2501111695/data/202512_mouse_testis_r1_leica/03_segmentation_ky/clustermap_ky"

# Python 环境
PY="/gpfs/share/home/2501111695/labShare/2501111695/envs/starfinder/bin/python"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/new_clustermaptest_dist_ky_fast_new.py"

# 参数
POSITION="${1:-Position1}"
TAG="${2:-AA1}"

# 默认参数（AA1 组合）
IC="0.08"        # cell_num_threshold
ID="4"           # dapi_grid_interval
IPF="0.003"      # pct_filter
ICR_XY="50"      # xy_radius
ICR_Z="10"       # z_radius

# 固定参数
Z_NUM=42
XY_SIZE=2048
ROUND_DAPI="round001"
DAPI_ROUND_NUM=1
GENE_CSV="${BASE}/01_data/genes.csv"
WINDOW_SIZE=300
OVERLAP=0.2

# ==================== 路径 ====================

WORK_DIR="${WORK_ROOT}/${POSITION}/_shared_prep"
TAG_WORK_DIR="${WORK_ROOT}/${POSITION}/${TAG}"
OUT_DIR="${OUT_ROOT}/${POSITION}/${TAG}_IC${IC}_ID${ID}_IPF${IPF}_R${ICR_XY}_${ICR_Z}"
LOG_DIR="${SCRIPT_DIR}/logs_single"

mkdir -p "${LOG_DIR}"

# ==================== 辅助函数 ====================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/pipeline.log"
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log "ERROR: 文件不存在：$1"
        exit 1
    fi
}

check_done() {
    if [[ -f "${OUT_DIR}/.cm_done" ]]; then
        log "已完成：${OUT_DIR}/.cm_done 存在"
        exit 0
    fi
}

# ==================== 检查 ====================

log "=========================================="
log "ClusterMap V9 单 Pipeline 运行"
log "=========================================="
log "Position: ${POSITION}"
log "TAG: ${TAG}"
log "参数：IC=${IC}, ID=${ID}, IPF=${IPF}, ICR=${ICR_XY},${ICR_Z}"
log ""

# 检查输入文件
log "检查输入文件..."
DAPI_DIR="${BASE}/01_data/${ROUND_DAPI}/${POSITION}"
GOOD_POINTS="${BASE}/02_registration/${POSITION}/goodPoints_max3d.csv"

check_file "${GENE_CSV}"
check_file "${GOOD_POINTS}"

if [[ ! -d "${DAPI_DIR}" ]]; then
    log "ERROR: DAPI 目录不存在：${DAPI_DIR}"
    exit 1
fi

DAPI_COUNT=$(ls "${DAPI_DIR}"/ch0*.tif 2>/dev/null | wc -l)
if [[ ${DAPI_COUNT} -eq 0 ]]; then
    log "ERROR: 未找到 DAPI 文件 (ch0*.tif) in ${DAPI_DIR}"
    exit 1
fi
log "找到 ${DAPI_COUNT} 个 DAPI 文件"

# 检查是否已完成
check_done

log ""
log "工作目录：${WORK_DIR}"
log "输出目录：${OUT_DIR}"
log ""

# ==================== PREP 阶段 ====================

run_prep() {
    log "=========================================="
    log "阶段 1: PREP (预处理)"
    log "=========================================="
    
    if [[ -f "${WORK_DIR}/manifest.json" ]]; then
        log "PREP 已完成：manifest.json 存在"
        return 0
    fi
    
    log "提交 PREP 作业..."
    
    sbatch \
        --job-name="cm_prep_${POSITION}" \
        --partition=AMCAmpere \
        --cpus-per-task=60 \
        --mem=480G \
        --time=08:00:00 \
        --output="${LOG_DIR}/prep.%j.out" \
        --error="${LOG_DIR}/prep.%j.err" \
        <<EOF
#!/bin/bash
set -e

log() { echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$*" | tee -a "${LOG_DIR}/prep.log"; }

log "PREP 开始"
log "工作目录：${WORK_DIR}"

module load anaconda/3
source activate starfinder

log "运行 PREP..."
${PY} ${SCRIPT} \\
    --stage prep \\
    -IP "${POSITION}" \\
    -IZ ${Z_NUM} \\
    -IC ${IC} \\
    -ID ${ID} \\
    -ICR "${ICR_XY},${ICR_Z}" \\
    -IXY ${XY_SIZE} \\
    -IDR ${DAPI_ROUND_NUM} \\
    -IPF ${IPF} \\
    -IEP T \\
    -IDir "${DAPI_DIR}" \\
    -Igood_points_max3d "${GOOD_POINTS}" \\
    -OP "${OUT_DIR}" \\
    --work_dir "${WORK_DIR}" \\
    --window_size ${WINDOW_SIZE} \\
    --overlap_percent ${OVERLAP}

log "PREP 完成"

# 验证输出
if [[ -f "${WORK_DIR}/manifest.json" ]]; then
    log "✓ manifest.json 已生成"
    N_TILES=\$(python -c "import json; print(json.load(open('${WORK_DIR}/manifest.json'))['n_tiles'])")
    log "✓ Tiles 数量：\${N_TILES}"
else
    log "ERROR: manifest.json 未生成"
    exit 1
fi
EOF
    
    log "PREP 作业已提交"
}

# ==================== TILE 阶段 ====================

run_tiles() {
    log "=========================================="
    log "阶段 2: TILE (并行分段)"
    log "=========================================="
    
    # 创建 TAG 特定的工作目录（symlink）
    if [[ ! -d "${TAG_WORK_DIR}" ]]; then
        log "创建 TAG 工作目录 symlink..."
        mkdir -p "$(dirname "${TAG_WORK_DIR}")"
        ln -s "_shared_prep" "${TAG_WORK_DIR}"
    fi
    
    # 获取 tile 数量
    N_TILES=$(python -c "import json; print(json.load(open('${WORK_DIR}/manifest.json'))['n_tiles'])")
    log "Tiles 数量：${N_TILES}"
    
    # 检查是否全部完成
    DONE_COUNT=$(ls "${TAG_WORK_DIR}/tile_results"/tile_*.pkl 2>/dev/null | wc -l)
    if [[ ${DONE_COUNT} -ge ${N_TILES} ]]; then
        log "TILE 已完成：${DONE_COUNT}/${N_TILES}"
        return 0
    fi
    log "已完成：${DONE_COUNT}/${N_TILES}"
    
    log "提交 TILE 数组作业..."
    
    sbatch \
        --job-name="cm_tile_${POSITION}_${TAG}" \
        --partition=AMCAmpere \
        --cpus-per-task=8 \
        --mem=64G \
        --time=24:00:00 \
        --array=0-$((N_TILES-1))%60 \
        --output="${LOG_DIR}/tile.%A_%a.out" \
        --error="${LOG_DIR}/tile.%A_%a.err" \
        <<EOF
#!/bin/bash
set -e

TILE_ID=\${SLURM_ARRAY_TASK_ID}

${PY} ${SCRIPT} \\
    --stage tile \\
    --tile_id \${TILE_ID} \\
    -IP "${POSITION}" \\
    -IZ ${Z_NUM} \\
    -IC ${IC} \\
    -ID ${ID} \\
    -ICR "${ICR_XY},${ICR_Z}" \\
    -IXY ${XY_SIZE} \\
    -IDR ${DAPI_ROUND_NUM} \\
    -IPF ${IPF} \\
    -IEP T \\
    -IDir "${DAPI_DIR}" \\
    -Igood_points_max3d "${GOOD_POINTS}" \\
    -OP "${OUT_DIR}" \\
    --work_dir "${TAG_WORK_DIR}" \\
    --window_size ${WINDOW_SIZE} \\
    --overlap_percent ${OVERLAP}
EOF
    
    log "TILE 作业已提交 (array=0-$((N_TILES-1))%60)"
}

# ==================== RESCUE 阶段 ====================

run_rescue() {
    log "=========================================="
    log "阶段 3: RESCUE (故障恢复)"
    log "=========================================="
    
    N_TILES=$(python -c "import json; print(json.load(open('${WORK_DIR}/manifest.json'))['n_tiles'])")
    
    # 检查完成情况
    for ROUND in {0..5}; do
        DONE_COUNT=$(ls "${TAG_WORK_DIR}/tile_results"/tile_*.pkl 2>/dev/null | wc -l)
        MISSING=$((N_TILES - DONE_COUNT))
        
        if [[ ${MISSING} -eq 0 ]]; then
            log "✓ 所有 tiles 完成：${DONE_COUNT}/${N_TILES}"
            return 0
        fi
        
        log "第 ${ROUND} 轮检查：缺失 ${MISSING} 个 tiles"
        
        # 找出缺失的 tiles
        MISSING_TILES=""
        for TID in $(seq 0 $((N_TILES-1))); do
            if [[ ! -f "${TAG_WORK_DIR}/tile_results/tile_$(printf '%05d' ${TID}).pkl" ]]; then
                MISSING_TILES="${MISSING_TILES} ${TID}"
            fi
        done
        
        if [[ -z "${MISSING_TILES}" ]]; then
            break
        fi
        
        log "重试 tiles:${MISSING_TILES}"
        
        # 提交 RETRY
        sbatch \
            --job-name="cm_retry_${POSITION}_${TAG}_r${ROUND}" \
            --partition=AMCAmpere \
            --cpus-per-task=8 \
            --mem=64G \
            --time=12:00:00 \
            --array="${MISSING_TILES// /,}%12 \
            --output="${LOG_DIR}/retry_r${ROUND}.%A_%a.out" \
            --error="${LOG_DIR}/retry_r${ROUND}.%A_%a.err" \
            <<EOF
#!/bin/bash
set -e

TILE_ID=\${SLURM_ARRAY_TASK_ID}

${PY} ${SCRIPT} \\
    --stage tile \\
    --tile_id \${TILE_ID} \\
    -IP "${POSITION}" \\
    -IZ ${Z_NUM} \\
    -IC ${IC} \\
    -ID ${ID} \\
    -ICR "${ICR_XY},${ICR_Z}" \\
    -IXY ${XY_SIZE} \\
    -IDR ${DAPI_ROUND_NUM} \\
    -IPF ${IPF} \\
    -IEP T \\
    -IDir "${DAPI_DIR}" \\
    -Igood_points_max3d "${GOOD_POINTS}" \\
    -OP "${OUT_DIR}" \\
    --work_dir "${TAG_WORK_DIR}" \\
    --window_size ${WINDOW_SIZE} \\
    --overlap_percent ${OVERLAP}
EOF
        
        log "RETRY 轮次 ${ROUND} 已提交，等待 30 秒..."
        sleep 30
    done
    
    # 最终检查
    DONE_COUNT=$(ls "${TAG_WORK_DIR}/tile_results"/tile_*.pkl 2>/dev/null | wc -l)
    MISSING=$((N_TILES - DONE_COUNT))
    
    if [[ ${MISSING} -gt 0 ]]; then
        log "WARNING: 仍有 ${MISSING} 个 tiles 失败，继续 STITCH"
    fi
}

# ==================== STITCH 阶段 ====================

run_stitch() {
    log "=========================================="
    log "阶段 4: STITCH (拼接)"
    log "=========================================="
    
    if [[ -f "${OUT_DIR}/.cm_done" ]]; then
        log "STITCH 已完成：.cm_done 存在"
        return 0
    fi
    
    log "提交 STITCH 作业..."
    
    sbatch \
        --job-name="cm_stitch_${POSITION}_${TAG}" \
        --partition=AMCAmpere \
        --cpus-per-task=16 \
        --mem=192G \
        --time=12:00:00 \
        --output="${LOG_DIR}/stitch.%j.out" \
        --error="${LOG_DIR}/stitch.%j.err" \
        --dependency=singleton \
        <<EOF
#!/bin/bash
set -e

log() { echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$*" | tee -a "${LOG_DIR}/stitch.log"; }

log "STITCH 开始"

module load anaconda/3
source activate starfinder

log "运行 STITCH..."
${PY} ${SCRIPT} \\
    --stage stitch \\
    -IP "${POSITION}" \\
    -IZ ${Z_NUM} \\
    -IC ${IC} \\
    -ID ${ID} \\
    -ICR "${ICR_XY},${ICR_Z}" \\
    -IXY ${XY_SIZE} \\
    -IDR ${DAPI_ROUND_NUM} \\
    -IPF ${IPF} \\
    -IEP T \\
    -IDir "${DAPI_DIR}" \\
    -Igood_points_max3d "${GOOD_POINTS}" \\
    -OP "${OUT_DIR}" \\
    --work_dir "${TAG_WORK_DIR}" \\
    --window_size ${WINDOW_SIZE} \\
    --overlap_percent ${OVERLAP}

log "STITCH 完成"

# 验证输出
if [[ -f "${OUT_DIR}/model_final.pkl" ]]; then
    log "✓ model_final.pkl 已生成"
else
    log "ERROR: model_final.pkl 未生成"
    exit 1
fi

if [[ -f "${OUT_DIR}/remain_reads.csv" ]]; then
    COUNT=\$(wc -l < "${OUT_DIR}/remain_reads.csv")
    log "✓ remain_reads.csv: \${COUNT} 行"
fi

if [[ -f "${OUT_DIR}/cell_center.csv" ]]; then
    COUNT=\$(wc -l < "${OUT_DIR}/cell_center.csv")
    log "✓ cell_center.csv: \${COUNT} 个细胞"
fi

# 标记完成
touch "${OUT_DIR}/.cm_done"
log "✓ 完成标记：.cm_done"

log "=========================================="
log "Pipeline 完成！"
log "输出：${OUT_DIR}"
log "=========================================="
EOF
    
    log "STITCH 作业已提交"
}

# ==================== 主流程 ====================

run_prep
# 等待 PREP 完成
while [[ ! -f "${WORK_DIR}/manifest.json" ]]; do
    log "等待 PREP 完成..."
    sleep 30
done
log "PREP 完成"

run_tiles
# 等待 TILE 完成（简单等待，实际应检查作业状态）
log "等待 TILE 完成（约 24 小时）..."
log "TILE 作业提交后继续检查..."

run_rescue

run_stitch

log ""
log "=========================================="
log "Pipeline 提交完成！"
log "=========================================="
log "作业状态：squeue -u \$USER"
log "日志目录：${LOG_DIR}"
log "输出目录：${OUT_DIR}"
log ""
log "注意：TILE 和 STITCH 作业需要较长时间完成"
log "使用 tail -f ${LOG_DIR}/*.out 查看实时日志"
