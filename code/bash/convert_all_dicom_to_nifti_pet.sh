#!/bin/bash

# Directorios principales
INPUT_DIR="/media/sol/Expansion/PET/Fleni/fgd_martin_pitossi-20251031T165436Z-1-001/fgd_martin_pitossi"
OUTPUT_DIR="/media/sol/Expansion/PET/Fleni/fgd_martin_pitossi-20251031T165436Z-1-001/fdg_martin_pitossi_nifti"
mkdir -p "$OUTPUT_DIR"

for subj_dir in "$INPUT_DIR"/*; do
    if [ -d "$subj_dir" ]; then
        subj_name=$(basename "$subj_dir")
        echo "Procesando sujeto: $subj_name"

        subj_out="$OUTPUT_DIR/$subj_name"
        mkdir -p "$subj_out"

        # Verificar si hay archivos con cabecera DICM
        dicom_check=$(find "$subj_dir" -type f -exec grep -Il "DICM" {} \; | head -n 1)

        if [ -n "$dicom_check" ]; then
            echo "Cabecera DICOM detectada, convirtiendo..."
        else
            echo " No se encontró 'DICM', se intentará conversión igual (puede ser DICOM válido sin esa marca)..."
        fi

        # Convertir con dcm2niix (usa profundidad 9 para buscar dentro de subcarpetas)
        dcm2niix -d 9 -z y -o "$subj_out" "$subj_dir"

        # Comprobar si se generó algún NIfTI
        nii_count=$(find "$subj_out" -type f -name "*.nii*" | wc -l)
        if [ "$nii_count" -gt 0 ]; then
            echo "Conversión exitosa ($nii_count archivos generados) para $subj_name"
        else
            echo "No se generaron archivos NIfTI para $subj_name"
        fi

        echo "---------------------------------------------"
    fi
done

echo "Conversión completa. Resultados en: $OUTPUT_DIR"

