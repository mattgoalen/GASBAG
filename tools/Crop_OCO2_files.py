"""
Crop_OCO2_files.py - file to crop OCO-2 data files to a given set of latitude limits, saves cropped files in 
                    "cropped_files" directory in the same location as the original file to be cropped.

Example usage:  python Crop_OCO2_files.py -y 2025 -m 5 -d 16 -f L1bSc 51,52
                # crop all L1bSc files from the 16th May 2025 to data between 51 and 52 degrees latitude
                python Crop_OCO2_files.py -y 2025 -f L2Met 51,52
                # crop all L2Met files from 2025 to data between 51 and 52 degrees latitude
                python Crop_OCO2_files.py -i /mnt/t/mg13/data/OCO2/L1bSc/2025/05/oco2_L1bScND_57802a_250514_B11213r_250606195253.h5 51,52
                # crop specified file (denoted by -i) to data between 51 and 52 degrees latitude
                python Crop_OCO2_files.py -m 05 -D /mnt/t/mg13
                # crop all OCO2 files from the month of May found within the data directory /mnt/t/mg13

-i -input_file      : Optional, specify an input file to be cropped
-y --year           : Optional, if input file unspecified, year of interest used to select file/s
-m --month          : Optional, if input file unspecified, month of interest used to select file/s
-d --day            : Optional, if input file unspecified, day of interest used to select file/s
-D --data_directory : Optional, if input file unspecified, top level data directory in which to search for files, defaults to /mnt/t/mg13
-f --file_type      : Optional, type of OCO-2 data to be cropped, if unspecified defaults to all types (*)
lat_limits          : Required, latitude limits used for cropping: lat_min,lat_max
"""

import argparse
import os
import glob
import numpy as np
import scipy
import h5py

def lat_lims(s):
    try:
        l_min, l_max = map(int, s.split(','))
        return l_min, l_max
    except:
        raise argparse.ArgumentTypeError("Must have two latitude limits: min,max")

# --------------------------
# Main program
# --------------------------

parser = argparse.ArgumentParser(description="Crop OCO-2 HDF5 files using h5py")
parser.add_argument("lat_limits", type=lat_lims)
parser.add_argument("-i", "--input_file", default=False)
parser.add_argument("-y", "--year", type=int, default=0)
parser.add_argument("-m", "--month", type=int, default=0)
parser.add_argument("-d", "--day", type=int, default=0)
parser.add_argument("-D", "--data_directory", default="/mnt/t/mg13")
parser.add_argument("-f", "--file_type", choices=["L1bSc", "L2Met", "*"], default="*")
args = parser.parse_args()

lat_min, lat_max = args.lat_limits
DATA_DIR = args.data_directory

year = "*" if args.year == 0 else f"{args.year}"
month = "*" if args.month == 0 else f"{args.month:02}"
day = "*" if args.day == 0 else f"{args.day:02}"
file_type = args.file_type

# File discovery
if not args.input_file:
    file_list = glob.glob(
        os.path.join(DATA_DIR, "**", f"oco2_{file_type}*_{year[-2:]}{month}{day}_*.h5"),
        recursive=True
    )
    file_list = [i for i in file_list if "cropped" not in i]
else:
    file_list = [args.input_file]


# --------------------------
# Loop over files
# --------------------------

for f in file_list:

    with h5py.File(f, "r") as h5in:

        # --- Read geometry group ---
        geom = h5in["SoundingGeometry"]
        lat = geom["sounding_latitude"][:]
        lon = geom["sounding_longitude"][:]

        sounding_dim, frame_dim = lat.shape

        # determine correct orientation
        lat_min, lat_max = (lat_min, lat_max) if lat[0,0] < lat[-1,-1] else (lat_max, lat_min)

        # find crop indices
        (sounding_0, frame_0) = [i[0] for i in np.where(lat - lat_min > 0)]
        (sounding_1, frame_1) = (np.array(lat.shape)-1 -
                                 np.array([i[0] for i in np.where(lat[::-1, ::-1] - lat_max < 0)]))

        min_sounding = int(min(sounding_0, sounding_1))
        max_sounding = int(max(sounding_0, sounding_1))

        # --------------------------
        # Create cropped output file
        # --------------------------

        out_dir = os.path.join(os.path.dirname(f), "cropped_files")
        os.makedirs(out_dir, exist_ok=True)

        new_name = os.path.basename(f).split(".")[0] + f"_cropped_{lat_min}{lat_max}.h5"
        out_path = os.path.join(out_dir, new_name)

        with h5py.File(out_path, "w") as h5out:

            # Copy all groups/datasets but slice sounding/frame dims
            def copy_group(g_in, g_out, sounding_dim, frame_dim,
                        min_sounding, max_sounding):

                for key, item in g_in.items():

                    # ---------------------------
                    # Case 1: Group → recurse
                    # ---------------------------
                    if isinstance(item, h5py.Group):
                        g_sub = g_out.create_group(key)
                        copy_group(
                            item, g_sub,
                            sounding_dim, frame_dim,
                            min_sounding, max_sounding,
                        )
                        continue

                    # ---------------------------
                    # Case 2: Dataset → crop
                    # ---------------------------
                    if isinstance(item, h5py.Dataset):
                        data = item[...]

                        # Identify axes with sizes equal to sounding_dim or frame_dim
                        axes = data.shape
                        sounding_axes = [i for i, s in enumerate(axes) if s == sounding_dim]
                        frame_axes    = [i for i, s in enumerate(axes) if s == frame_dim]

                        # If neither dimension matches, copy unchanged
                        if not sounding_axes and not frame_axes:
                            dset_out = g_out.create_dataset(key, data=data, compression="gzip")
                            for a,v in item.attrs.items(): dset_out.attrs[a] = v
                            continue

                        # Build full-slicing list
                        slicer = [slice(None)] * data.ndim

                        # Apply sounding-dim crop wherever it appears
                        for ax in sounding_axes:
                            slicer[ax] = slice(min_sounding, max_sounding)

                        cropped = data[tuple(slicer)]

                        # Write out cropped
                        dset_out = g_out.create_dataset(key, data=cropped, compression="gzip")

                        # Copy attributes
                        for a, v in item.attrs.items():
                            dset_out.attrs[a] = v

            copy_group(
                h5in, h5out,
                sounding_dim=sounding_dim,
                frame_dim=frame_dim,
                min_sounding=min_sounding,
                max_sounding=max_sounding,
            )

    print(f"Saved cropped file: {out_path}")
