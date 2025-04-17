#!/usr/bin/env python
"""
CLI tool for cdSirius Compound Discoverer integration
"""
import os
import sys
import shutil
import argparse
import json
import subprocess
from pathlib import Path


def find_python_executable():
    """Find the Python executable that's running this script."""
    return sys.executable


def find_sirius_executable():
    """Try to locate the Sirius executable installation."""
    default_paths = [
        r"C:\Program Files\sirius\sirius.exe",
        r"C:\Program Files\sirius-gui\sirius.exe"
    ]
    
    for path in default_paths:
        if os.path.exists(path):
            return path
    
    return None


def setup_cd_integration(cd_path, sirius_path, python_path=None):
    """
    Set up the integration with Compound Discoverer.
    
    Args:
        cd_path: Path to Compound Discoverer installation
        python_path: Path to Python executable
        sirius_path: Path to Sirius executable
    """
    if not python_path:
        python_path = find_python_executable()
    
    # Make sure the paths exist
    if not os.path.exists(cd_path):
        print(f"Error: Compound Discoverer path does not exist: {cd_path}")
        return False
    
    if not os.path.exists(python_path):
        print(f"Error: Python executable path does not exist: {python_path}")
        return False
    
    if not os.path.exists(sirius_path):
        print(f"Warning: Sirius executable path does not exist: {sirius_path}")
        print("The node.json will be updated with this path anyway, but please verify it's correct.")
    
    # Create the Scripts/cdSirius directory if it doesn't exist
    scripts_path = os.path.join(cd_path, "Tools", "Scripts", "cdSirius")
    os.makedirs(scripts_path, exist_ok=True)
    
    # Get package directory
    package_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Copy node.json and icon files
    files_to_copy = ["node.json", "IMG_16x16.png", "IMG_32x32.png"]
    source_dir = os.path.dirname(package_dir)  # parent of cdsirius package directory
    
    for filename in files_to_copy:
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(scripts_path, filename)
        
        if os.path.exists(source_file):
            shutil.copy2(source_file, target_file)
            print(f"Copied {filename} to {target_file}")
        else:
            print(f"Warning: Could not find {filename} in package directory")
    
    # Update node.json with Python path
    node_json_path = os.path.join(scripts_path, "node.json")
    if not os.path.exists(node_json_path):
        print("Error: Failed to find node.json in the target directory")
        return False
    
    try:
        with open(node_json_path, 'r') as f:
            node_data = json.load(f)
        
        # Update Python path in ScriptProcessorArguments and script path
        if "ScriptProcessorArguments" in node_data:
            # Update the ExecutablePath
            if "ExecutablePath" in node_data["ScriptProcessorArguments"]:
                node_data["ScriptProcessorArguments"]["ExecutablePath"] = python_path
            
            # Update the ExecutableCommandLineArguments to point to the bootstrap script
            if "ExecutableCommandLineArguments" in node_data["ScriptProcessorArguments"]:
                bootstrap_path = os.path.join(scripts_path, "cdSirius.py").replace("\\", "\\\\")
                cmd_args = node_data["ScriptProcessorArguments"]["ExecutableCommandLineArguments"]
                # Replace just the script path part while preserving %NODEARGS%
                parts = cmd_args.split()
                if len(parts) > 0:
                    node_data["ScriptProcessorArguments"]["ExecutableCommandLineArguments"] = f"{bootstrap_path} %NODEARGS%"
        
        # Update the Sirius Program Path in Parameters if provided
        if sirius_path and "Parameters" in node_data:
            for param in node_data["Parameters"]:
                if param.get("Name") == "Sirius Program Path":
                    param["Default"] = sirius_path
                    print(f"Updated Sirius Program Path in node.json to {sirius_path}")
                    break
        
        # Write updated node.json
        with open(node_json_path, 'w') as f:
            json.dump(node_data, f, indent=2)
        
        print(f"Updated Python path in node.json to {python_path}")
        
        # Create cdSirius.py in the target directory that imports from the package
        bootstrap_script = os.path.join(scripts_path, "cdSirius.py")
        with open(bootstrap_script, 'w') as f:
            f.write(f'''#!/usr/bin/env python
"""
Bootstrap script that imports from the installed cdsirius package
"""
import sys
import os
from cdsirius import cdSirius

if __name__ == "__main__":
    sys.exit(cdSirius.main())
''')
        print(f"Created bootstrap script at {bootstrap_script}")
        
        print("\nCompound Discoverer integration completed successfully!")
        print("\nFinal steps:")
        print("1. Launch Compound Discoverer 3.3")
        print("2. Navigate to Help -> License Manager")
        print("3. Run 'Scan for Missing Features'")
        print("4. Restart Compound Discoverer")
        
        return True
        
    except Exception as e:
        print(f"Error updating node.json: {str(e)}")
        return False


def main():
    parser = argparse.ArgumentParser(description="cdSirius Compound Discoverer Integration Tool")
    
    parser.add_argument("cd_path", 
                        help="Path to Compound Discoverer installation directory")
    
    parser.add_argument("sirius-path",
                        help="Path to Sirius executable (e.g., C:\\Program Files\\sirius\\sirius.exe)")
    
    parser.add_argument("--python-path", 
                        help="Path to Python executable (default: current Python)")
    
    
    
    args = parser.parse_args()
    
    success = setup_cd_integration(args.cd_path, args.sirius_path, args.python_path)
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()