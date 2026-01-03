
import h5py

def print_structure(name, obj):
    print(name)
    if isinstance(obj, h5py.Dataset):
        print(f"  Shape: {obj.shape}, Dtype: {obj.dtype}")
        # Print attributes if any
        for key, val in obj.attrs.items():
            print(f"  Attr {key}: {val}")

def inspect_file(filename):
    try:
        with h5py.File(filename, 'r') as f:
            print(f"Inspecting {filename}...")
            f.visititems(print_structure)
    except Exception as e:
        print(f"Error inspecting file: {e}")

if __name__ == "__main__":
    inspect_file("output.h5")
