#!/usr/bin/env python3
"""
Validate PETSc solver configuration in YAML files
Checks that all solver types are valid and lowercase
"""

import yaml
import sys
from pathlib import Path

# Valid PETSc types (must be lowercase)
VALID_KSP_TYPES = {
    'gmres', 'fgmres', 'bcgs', 'cg', 'bicg', 'tfqmr', 
    'richardson', 'chebyshev', 'preonly', 'minres'
}

VALID_PC_TYPES = {
    'pbjacobi', 'bjacobi', 'asm', 'jacobi', 'ilu', 'sor', 
    'lu', 'none', 'cholesky', 'icc', 'mg', 'hypre', 'gamg'
}

def validate_config(config_file):
    """Validate a single configuration file"""
    print(f"Validating: {config_file}")
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        if not config:
            print(f"  ⚠️  Warning: Empty or invalid YAML file")
            return True
        
        # Check if solver section exists
        if 'solver' not in config:
            print(f"  ℹ️  Info: No solver section (will use defaults)")
            return True
        
        solver = config['solver']
        errors = []
        warnings = []
        
        # Check KSP type
        if 'ksp_type' in solver:
            ksp_type = solver['ksp_type']
            if not isinstance(ksp_type, str):
                errors.append(f"ksp_type must be a string, got {type(ksp_type)}")
            elif ksp_type != ksp_type.lower():
                errors.append(f"ksp_type '{ksp_type}' must be lowercase (use '{ksp_type.lower()}')")
            elif ksp_type.lower() not in VALID_KSP_TYPES:
                warnings.append(f"ksp_type '{ksp_type}' may not be valid. Common types: {', '.join(sorted(list(VALID_KSP_TYPES)[:5]))}")
        
        # Check PC type
        if 'preconditioner' in solver:
            pc_type = solver['preconditioner']
            if not isinstance(pc_type, str):
                errors.append(f"preconditioner must be a string, got {type(pc_type)}")
            elif pc_type != pc_type.lower():
                errors.append(f"preconditioner '{pc_type}' must be lowercase (use '{pc_type.lower()}')")
            elif pc_type.lower() not in VALID_PC_TYPES:
                if pc_type.upper() == 'PBP':
                    errors.append(f"preconditioner 'PBP' is NOT a valid PETSc type! Use 'pbjacobi' instead")
                else:
                    warnings.append(f"preconditioner '{pc_type}' may not be valid. Common types: {', '.join(sorted(list(VALID_PC_TYPES)[:5]))}")
        
        # Report results
        if errors:
            print(f"  ❌ ERRORS found:")
            for err in errors:
                print(f"     - {err}")
            return False
        elif warnings:
            print(f"  ⚠️  WARNINGS:")
            for warn in warnings:
                print(f"     - {warn}")
            print(f"  ✓  No critical errors (but check warnings)")
            return True
        else:
            print(f"  ✓  Valid configuration")
            return True
            
    except Exception as e:
        print(f"  ❌ ERROR: Failed to parse file: {e}")
        return False

def main():
    """Validate all YAML configuration files"""
    print("=" * 70)
    print("PETSc Solver Configuration Validator")
    print("=" * 70)
    print()
    
    # Find all YAML files in config directory
    config_dir = Path(__file__).parent / 'config'
    if not config_dir.exists():
        print(f"Error: Config directory not found: {config_dir}")
        sys.exit(1)
    
    yaml_files = list(config_dir.glob('*.yaml')) + list(config_dir.glob('*.yml'))
    
    if not yaml_files:
        print("No YAML files found in config directory")
        sys.exit(1)
    
    print(f"Found {len(yaml_files)} YAML configuration file(s)\n")
    
    all_valid = True
    for yaml_file in sorted(yaml_files):
        if not validate_config(yaml_file):
            all_valid = False
        print()
    
    print("=" * 70)
    if all_valid:
        print("✓ All configuration files are valid!")
        print("=" * 70)
        sys.exit(0)
    else:
        print("❌ Some configuration files have errors!")
        print("=" * 70)
        print("\nPlease fix the errors and re-run validation.")
        print("See docs/PETSC_SOLVER_GUIDE.md for valid solver types.")
        sys.exit(1)

if __name__ == '__main__':
    main()
