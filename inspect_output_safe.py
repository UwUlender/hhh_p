
import h5py
import time
import os

def print_structure(name, obj):
    print(name)
    if isinstance(obj, h5py.Dataset):
        print(f"  Shape: {obj.shape}, Dtype: {obj.dtype}")
        # Print attributes if any
        for key, val in obj.attrs.items():
            print(f"  Attr {key}: {val}")

def inspect_file(filename):
    print(f"Inspecting {filename}...")
    try:
        # Try waiting a bit or opening in read only mode strictly
        with h5py.File(filename, 'r', swmr=True) as f: # Use SWMR mode if possible
             f.visititems(print_structure)
    except Exception as e:
        print(f"Error inspecting file: {e}")
        try:
             # Fallback to standard open
             with h5py.File(filename, 'r') as f:
                f.visititems(print_structure)
        except Exception as e2:
             print(f"Error inspecting file (fallback): {e2}")

if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    inspect_file("output.h5")
