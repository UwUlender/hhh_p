
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

def process_output(filename, output_dir="plots"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Processing {filename}...")
    try:
        # Use file locking false just in case
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
        
        with h5py.File(filename, 'r') as f:
            # Get mesh
            x = f['mesh/x_coords'][:]
            
            # Find all steps
            steps = []
            for key in f['data'].keys():
                if key.startswith('step_'):
                    steps.append(int(key.split('_')[1]))
            
            steps.sort()
            if not steps:
                print("No steps found.")
                return

            last_step = steps[-1]
            print(f"Found {len(steps)} steps. Last step: {last_step}")
            
            # Get data for last step
            group = f[f'data/step_{last_step}']
            
            n_e = group['n_0'][0, :] # Assuming 1D row 0
            n_Ar_plus = group['n_1'][0, :]
            n_Ar_meta = group['n_2'][0, :]
            n_eps = group['n_eps'][0, :]
            phi = group['phi'][0, :]
            
            # Plot Densities
            plt.figure(figsize=(10, 6))
            plt.plot(x, n_e, label='Electron Density (n_e)')
            plt.plot(x, n_Ar_plus, label='Ar+ Density')
            plt.plot(x, n_Ar_meta, label='Ar Metastable Density')
            plt.xlabel('Position (m)')
            plt.ylabel('Density (m^-3)')
            plt.title(f'Species Densities at Step {last_step}')
            plt.legend()
            plt.grid(True)
            plt.yscale('log')
            plt.savefig(f"{output_dir}/densities_step_{last_step}.png")
            plt.close()
            
            # Plot Potential
            plt.figure(figsize=(10, 6))
            plt.plot(x, phi, label='Potential')
            plt.xlabel('Position (m)')
            plt.ylabel('Potential (V)')
            plt.title(f'Potential at Step {last_step}')
            plt.legend()
            plt.grid(True)
            plt.savefig(f"{output_dir}/potential_step_{last_step}.png")
            plt.close()
            
            # Plot Electron Energy
            plt.figure(figsize=(10, 6))
            plt.plot(x, n_eps, label='Electron Energy Density')
            # Calculate mean energy if n_e > 0
            # mean_energy = n_eps / n_e
            # Safe division
            with np.errstate(divide='ignore', invalid='ignore'):
                 mean_energy = np.where(n_e > 1e10, n_eps / n_e, 0)
                 
            plt.plot(x, mean_energy, label='Mean Electron Energy (eV)')
            plt.xlabel('Position (m)')
            plt.ylabel('Energy (eV)')
            plt.title(f'Electron Energy at Step {last_step}')
            plt.legend()
            plt.grid(True)
            plt.savefig(f"{output_dir}/energy_step_{last_step}.png")
            plt.close()
            
            print(f"Plots saved to {output_dir}/")
            
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    process_output("output.h5")
