#!/bin/bash -l
#$ -P paxlab
#$ -N edge_nfrags_unsup
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l h_rt=08:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/edge_nfrags_unsup.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/edge_nfrags_unsup.$JOB_ID.err

set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_PATH="analysis/src/build_tissue/build_tissue_barcodes_edge_nfrags_plots.R"
OUT_BARCODE_DIR="${OUT_BARCODE_DIR:-$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes}"
OUT_PLOT_DIR="${OUT_PLOT_DIR:-$PROJECT_ROOT/analysis/plots/variant_qc/edge_effect_nfrags}"

if [[ -f /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
fi
module load R

cd "${PROJECT_ROOT}"

mkdir -p "${OUT_BARCODE_DIR}" "${OUT_PLOT_DIR}"

echo "[start] unsupervised edge runs: axes=col,row"

for axis in col row; do
  export EDGE_AXIS=${axis}
  export FORCE_RECOUNT=0
  echo "[run] EDGE_AXIS=${EDGE_AXIS}"
  Rscript "${SCRIPT_PATH}"

  # Archive results for this axis
  axis_dir_barcode="${OUT_BARCODE_DIR}/axis_${axis}"
  axis_dir_plot="${OUT_PLOT_DIR}/axis_${axis}"
  mkdir -p "${axis_dir_barcode}" "${axis_dir_plot}"

    # Move object folders and summary into axis-specific folder
    for obj in deepseq_488B deepseq_489 lowseq_488B lowseq_489; do
        mv -f "${OUT_BARCODE_DIR}/${obj}" "${axis_dir_barcode}/" 2>/dev/null || true
    done
  mv -f "${OUT_BARCODE_DIR}/edge_effect_nfrags_thresholds.tsv" "${axis_dir_barcode}/edge_effect_nfrags_thresholds_axis_${axis}.tsv" 2>/dev/null || true

    # Move per-object plot folders
    for obj in deepseq_488B deepseq_489 lowseq_488B lowseq_489; do
        mv -f "${OUT_PLOT_DIR}/${obj}" "${axis_dir_plot}/" 2>/dev/null || true
    done

done

echo "[compare] choosing axis per dataset/tissue by larger number of edge cells"

# Build combined preference table
python3 - <<'PY'
import pandas as pd, glob
import os
base=os.environ.get('OUT_BARCODE_DIR','Data/01_inputs/barcodes/tissue_barcodes')
dfs={}
for axis in ('col','row'):
    p=f'{base}/axis_{axis}/edge_effect_nfrags_thresholds_axis_{axis}.tsv'
    try:
        dfs[axis]=pd.read_csv(p,sep='\t')
    except Exception:
        dfs[axis]=pd.DataFrame()

keys=set()
for d in dfs.values():
    for idx,row in d.iterrows():
        keys.add((row['dataset'],str(row['tissue'])))

out=[]
for dataset,tissue in sorted(keys):
    a_col=dfs['col']
    a_row=dfs['row']
    rcol=a_col[(a_col['dataset']==dataset)&(a_col['tissue']==tissue)]
    rrow=a_row[(a_row['dataset']==dataset)&(a_row['tissue']==tissue)]
    ecol=int(rcol['edge_cells'].values[0]) if len(rcol) else 0
    erow=int(rrow['edge_cells'].values[0]) if len(rrow) else 0
    prefer='col' if ecol>=erow else 'row'
    out.append({'dataset':dataset,'tissue':tissue,'edge_col':ecol,'edge_row':erow,'prefer_axis':prefer})

pd.DataFrame(out).to_csv(f'{base}/axis_choice_summary.tsv',sep='\t',index=False)
print('Wrote axis_choice_summary.tsv')
PY

echo "[publish] deploy preferred barcode lists and plots"

python3 - <<'PY'
import pandas as pd, shutil, os
base=os.environ.get('OUT_BARCODE_DIR','Data/01_inputs/barcodes/tissue_barcodes')
plotbase=os.environ.get('OUT_PLOT_DIR','analysis/plots/variant_qc/edge_effect_nfrags')
df=pd.read_csv(f'{base}/axis_choice_summary.tsv',sep='\t')
obj_map={
    ('deepseq','488B'):'deepseq_488B',
    ('deepseq','489'):'deepseq_489',
    ('lowseq','488B'):'lowseq_488B',
    ('lowseq','489'):'lowseq_489',
}
for idx,row in df.iterrows():
    a=row['prefer_axis']
    ds=row['dataset']; ts=row['tissue']
    obj=obj_map.get((ds,str(ts)))
    if obj is None:
        continue
    srcdir=f'{base}/axis_{a}'
    srcplot=f'{plotbase}/axis_{a}'

    # copy chosen barcode object folder
    s_obj=os.path.join(srcdir,obj)
    d_obj=os.path.join(base,obj)
    if os.path.exists(s_obj):
        if os.path.exists(d_obj):
            shutil.rmtree(d_obj)
        shutil.copytree(s_obj,d_obj)

    # copy chosen plot object folder
    s_plot=os.path.join(srcplot,obj)
    d_plot=os.path.join(plotbase,obj)
    if os.path.exists(s_plot):
        if os.path.exists(d_plot):
            shutil.rmtree(d_plot)
        shutil.copytree(s_plot,d_plot)
print('Published preferred files')
PY

echo "[done] unsupervised edge selection finished"
