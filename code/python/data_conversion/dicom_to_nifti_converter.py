import os
import argparse
import shutil
import subprocess
import logging
#!/usr/bin/env python3
"""
Convert DICOM files in a directory tree to NIfTI, preserving structure.

Default dicom root (as requested):
/home/martin/data/UNSAM/CovidProject/ConsorcioNacional/Datos/fgd_martin_pitossi/

Output root can be provided; default will be created next to the input root.
This script will try to use dcm2niix (recommended). If not found and
dicom2nifti is installed, it will use that as a fallback.

Usage:
    python3 dicom_to_nifti.py --dicom-root /path/to/dicom_root --output-root /path/to/output_root
"""


try:
    import pydicom
except Exception:
    pydicom = None

try:
    import dicom2nifti
except Exception:
    dicom2nifti = None

# Default path from the user's prompt
DEFAULT_DICOM_ROOT = "/home/martin/data/UNSAM/CovidProject/ConsorcioNacional/Datos/fgd_martin_pitossi/"
DEFAULT_NIFTI_ROOT = "/home/martin/data/UNSAM/CovidProject/ConsorcioNacional/Datos/fgd_martin_pitossi_nifti/"

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def is_probably_dicom_file(path):
    """Quick check: extension or try to read header with pydicom if available."""
    if not os.path.isfile(path):
        return False
    ext = os.path.splitext(path)[1].lower()
    if ext in (".dcm", ""):  # many DICOMs have no extension
        if pydicom:
            try:
                # read minimal header (no pixel data)
                pydicom.dcmread(path, stop_before_pixels=True, force=True)
                return True
            except Exception:
                return False
        else:
            # fallback: check for 'DICM' at byte offset 128 (not always present)
            try:
                with open(path, "rb") as f:
                    f.seek(128)
                    header = f.read(4)
                    return header == b"DICM"
            except Exception:
                return False
    return False


def is_dicom_dir(dirpath, max_files_check=20):
    """Return True if directory contains at least one probable DICOM file."""
    try:
        entries = os.listdir(dirpath)
    except Exception:
        return False
    files_checked = 0
    for name in entries:
        if files_checked >= max_files_check:
            break
        fp = os.path.join(dirpath, name)
        if os.path.isfile(fp):
            files_checked += 1
            if is_probably_dicom_file(fp):
                return True
    return False


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def convert_with_dcm2niix(input_dir, out_dir, extra_args=None):
    args = ["dcm2niix", "-z", "y", "-o", out_dir, "-f", "%p_%s"]
    if extra_args:
        args.extend(extra_args)
    args.append(input_dir)
    logging.info("Running: %s", " ".join(args))
    res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        logging.error("dcm2niix failed for %s: %s", input_dir, res.stderr.strip())
        return False
    logging.debug(res.stdout)
    return True


def convert_with_dicom2nifti(input_dir, out_dir):
    try:
        logging.info("Using dicom2nifti for %s -> %s", input_dir, out_dir)
        dicom2nifti.convert_directory(input_dir, out_dir, compression=True, reorient=False)
        return True
    except Exception as e:
        logging.error("dicom2nifti failed for %s: %s", input_dir, e)
        return False


def main(dicom_root, output_root):
    dicom_root = os.path.abspath(dicom_root)
    output_root = os.path.abspath(output_root)

    if not os.path.exists(dicom_root):
        logging.error("DICOM root does not exist: %s", dicom_root)
        return 1

    ensure_dir(output_root)

    use_dcm2niix = shutil.which("dcm2niix") is not None
    if use_dcm2niix:
        logging.info("Found dcm2niix on PATH; will use it for conversion.")
    else:
        logging.info("dcm2niix not found. Will try dicom2nifti if installed.")

    processed = 0
    skipped = 0

    for root, dirs, files in os.walk(dicom_root):
        # skip hidden directories
        if os.path.basename(root).startswith("."):
            continue

        if not files:
            continue

        if not is_dicom_dir(root):
            skipped += 1
            continue

        rel = os.path.relpath(root, dicom_root)
        dest_dir = os.path.join(output_root, rel)
        ensure_dir(dest_dir)

        success = False
        if use_dcm2niix:
            success = convert_with_dcm2niix(root, dest_dir)
        else:
            success = convert_with_dicom2nifti(root, dest_dir)

        if success:
            logging.info("Converted: %s -> %s", root, dest_dir)
            processed += 1
        else:
            logging.error("Failed to convert: %s", root)

    logging.info("Done. Processed: %d directories. Skipped (non-DICOM): %d", processed, skipped)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch convert DICOM folders to NIfTI preserving folder structure.")
    parser.add_argument("--dicom-root", default=DEFAULT_DICOM_ROOT, help="Root folder containing DICOM subfolders")
    parser.add_argument("--output-root", default=DEFAULT_NIFTI_ROOT, help="Root folder where NIfTI outputs will be written (preserves structure).")
    args = parser.parse_args()

    if args.output_root:
        out_root = args.output_root
    else:
        parent = os.path.dirname(os.path.abspath(args.dicom_root.rstrip("/")))
        out_root = os.path.join(parent, os.path.basename(args.dicom_root.rstrip("/")) + "_nifti")

    raise_code = main(args.dicom_root, out_root)
    if raise_code:
        raise SystemExit(raise_code)