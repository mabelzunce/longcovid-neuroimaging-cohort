#!/bin/bash
# ==========================================================
# 🔹 MRtrix3 DWI Multishell Preprocessing Pipeline (1 sujeto)
# ==========================================================
# Uso:
#   bash run_pipeline_subject.sh sub-CP0124
# Guarda resultados en: /Preprocessed/<SUJETO>/
# ==========================================================

set -e
trap 'echo "Error en $subj. Revisar log: $LOGFILE"; exit 1' ERR

MULTISHELL_DIR="/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/DTI/Nifti/multishell"
OUT_BASE="/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/DTI/Preprocessed"

export MRTRIX_NTHREADS=16
export OMP_NUM_THREADS=16

subj=$1
subj_dir="$MULTISHELL_DIR/$subj"
dwi_dir="$subj_dir/dwi"

if [[ ! -d "$dwi_dir" ]]; then
  echo "No existe carpeta DWI para $subj"
  exit 1
fi

echo "Procesando $subj ..."

# ===================== ARCHIVOS BASE =====================
DWI_AP=$(ls "$dwi_dir"/${subj}_acq-b1000b2000_dir-AP_dwi.nii.gz 2>/dev/null)
BVALS="${DWI_AP%.nii.gz}.bval"
BVECS="${DWI_AP%.nii.gz}.bvec"
B0_PA=$(ls "$dwi_dir"/${subj}_acq-b0_dir-PA_dwi.nii.gz 2>/dev/null)
JSON_PA=$(ls "$dwi_dir"/${subj}_acq-b0_dir-PA_dwi.json 2>/dev/null)

# Verificación de existencia detallada
missing=()
[[ ! -f "$DWI_AP" ]] && missing+=("DWI_AP ($DWI_AP)")
[[ ! -f "$BVALS" ]] && missing+=("BVALS ($BVALS)")
[[ ! -f "$BVECS" ]] && missing+=("BVECS ($BVECS)")
[[ ! -f "$B0_PA" ]] && missing+=("B0_PA ($B0_PA)")
[[ ! -f "$JSON_PA" ]] && missing+=("JSON_PA ($JSON_PA)")

if [[ ${#missing[@]} -gt 0 ]]; then
  echo "⚠️  Archivos faltantes para $subj:"
  for f in "${missing[@]}"; do echo "   - $f"; done
  echo "   ➤ Sujeto $subj omitido por archivos incompletos."
  exit 0
fi

# Leer readout time desde JSON
READOUT=$(grep -o '"TotalReadoutTime": *[0-9.]*' "$JSON_PA" | grep -o '[0-9.]*')
if [[ -z "$READOUT" ]]; then
  echo "⚠️  No se encontró TotalReadoutTime para $subj — saltando."
  exit 0
fi

# ===================== DIRECTORIOS Y LOG =====================
OUT_DIR="$OUT_BASE/$subj"
mkdir -p "$OUT_DIR"
LOGFILE="$OUT_DIR/pipeline_log.txt"

echo "===========================================" | tee "$LOGFILE"
echo "🧠 Preprocesamiento DTI - $subj" | tee -a "$LOGFILE"
echo "Fecha: $(date)" | tee -a "$LOGFILE"
echo "Readout time: $READOUT s" | tee -a "$LOGFILE"
echo "===========================================" | tee -a "$LOGFILE"

start_time=$(date +%s)

# Archivos intermedios
DWI_AP_MIF="$OUT_DIR/ap.mif"
DWI_PA_B0_MIF="$OUT_DIR/pa_b0.mif"
DWI_AP_MIF_DEN="$OUT_DIR/ap_den.mif"
DWI_AP_MIF_DEN_DEG="$OUT_DIR/ap_den_deg.mif"
DWI_PA_B0_MIF_DEN="$OUT_DIR/pa_b0_den.mif"
DWI_PA_B0_MIF_DEN_DEG="$OUT_DIR/pa_b0_den_deg.mif"
OUT_TOPUP_EDDY="$OUT_DIR/dwi_preproc_APPA.mif"
OUT_MASK="$OUT_DIR/dwi_preproc_APPA_mask.mif"
OUT_BIAS="$OUT_DIR/dwi_preproc_APPA_biascorr.mif"
GRAD_MATRIX="$OUT_DIR/grad_checked.txt"
PE_AP="ap"

# ===================== 1) Importación =====================
echo "[1/8] Importando DWI..." | tee -a "$LOGFILE"
mrconvert "$DWI_AP" -fslgrad "$BVECS" "$BVALS" -force "$DWI_AP_MIF"
mrconvert "$B0_PA" -force "$DWI_PA_B0_MIF"

# ===================== 2) Denoise & Degibbs =====================
echo "[2/8] Denoising + Gibbs..." | tee -a "$LOGFILE"
dwidenoise "$DWI_AP_MIF" "$DWI_AP_MIF_DEN" -nthreads 16 -force
mrdegibbs "$DWI_AP_MIF_DEN" "$DWI_AP_MIF_DEN_DEG" -nthreads 16 -force
dwidenoise "$DWI_PA_B0_MIF" "$DWI_PA_B0_MIF_DEN" -nthreads 16 -force
mrdegibbs "$DWI_PA_B0_MIF_DEN" "$DWI_PA_B0_MIF_DEN_DEG" -nthreads 16 -force

# ===================== 3) TOPUP + EDDY =====================
echo "[3/8] TOPUP + EDDY..." | tee -a "$LOGFILE"
dwifslpreproc "$DWI_AP_MIF_DEN_DEG" "$OUT_TOPUP_EDDY" \
  -pe_dir $PE_AP \
  -rpe_pair -se_epi "$DWI_PA_B0_MIF_DEN_DEG" \
  -readout_time "$READOUT" \
  -eddyqc_text "$OUT_DIR/eddy_qc_text" \
  -eddy_options "--slm=linear --repol --very_verbose" \
  -force 2>&1 | tee -a "$LOGFILE"

if [[ ! -f "$OUT_TOPUP_EDDY" ]]; then
  echo "dwifslpreproc falló para $subj, saltando..." | tee -a "$LOGFILE"
  exit 1
fi

# ===================== 4) Máscara =====================
echo "[4/8] Creando máscara..." | tee -a "$LOGFILE"
dwi2mask "$OUT_TOPUP_EDDY" "$OUT_MASK" -nthreads 16 -force

# ===================== 5) Bias Field =====================
echo "[5/8] Corrigiendo bias field..." | tee -a "$LOGFILE"
dwibiascorrect ants "$OUT_TOPUP_EDDY" "$OUT_BIAS" -mask "$OUT_MASK" -nthreads 16 -force

# ===================== 6) Gradiente QC =====================
echo "[6/8] Chequeando gradientes..." | tee -a "$LOGFILE"
dwigradcheck "$OUT_BIAS" -export_grad_mrtrix "$GRAD_MATRIX" -force

# ===================== 7) Tensor + Mapas =====================
echo "[7/8] Calculando tensor y mapas..." | tee -a "$LOGFILE"
TENSOR="$OUT_DIR/tensor.mif"
FA="$OUT_DIR/fa.mif"
MD="$OUT_DIR/md.mif"
RD="$OUT_DIR/rd.mif"
AD="$OUT_DIR/ad.mif"

dwi2tensor "$OUT_BIAS" "$TENSOR" -mask "$OUT_MASK" -grad "$GRAD_MATRIX" -nthreads 16 -force
tensor2metric "$TENSOR" -fa "$FA" -adc "$MD" -rd "$RD" -ad "$AD" -nthreads 16 -force

# ===================== 8) Exportar NIfTI =====================
echo "[8/8] Exportando a NIfTI..." | tee -a "$LOGFILE"
for metric in fa md rd ad; do
  mrconvert "$OUT_DIR/${metric}.mif" "$OUT_DIR/${metric}.nii.gz" -force
done

# ===================== FINAL =====================
end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "-------------------------------------------" | tee -a "$LOGFILE"
echo "✅ Sujeto $subj procesado exitosamente." | tee -a "$LOGFILE"
echo "Duración total: $(date -ud "@$elapsed" +'%H:%M:%S')" | tee -a "$LOGFILE"
echo "-------------------------------------------" | tee -a "$LOGFILE"



