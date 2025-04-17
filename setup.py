from setuptools import setup

setup(
    # Note: Both pysirius and pyeds need to be installed separately from GitHub
    # pysirius: pip install git+https://github.com/sirius-ms/sirius-client-openAPI@v0.0.0#subdirectory=client-api_python/generated
    # pyeds: pip install git+https://github.com/thermofisherlsms/pyeds
    
    # All dependencies needed for the package
    install_requires=[
        "pandas",
        # pyeds is installed separately from GitHub
        "rdkit-pypi",
        "molmass",
        # pysirius is not included here as it must be installed from GitHub
        "numpy",
        "matplotlib",
        "scipy",
        "openpyxl",
        "requests"
    ],
    # Post-install message to notify users about external dependencies
    # This will appear in the terminal after installation
    classifiers=[
        "Framework :: Setuptools Plugin",
    ],
    # This line enables the post-install message
    entry_points={
        "distutils.setup_keywords": [
            "install_requires=setuptools.dist:check_requirements",
        ],
    },
)

# Attempt to display a post-install message to inform users about dependencies
try:
    import sys
    import os
    
    # If this is being installed (not just imported)
    if any(arg.startswith('install') for arg in sys.argv):
        print("\n" + "="*80)
        print("IMPORTANT: cdSirius requires the following packages which must be installed separately:")
        print("1. pySirius - install with the command:")
        print("   pip install git+https://github.com/sirius-ms/sirius-client-openAPI#subdirectory=client-api_python/generated")
        print("\n2. pyEDS - install with the command:")
        print("   pip install git+https://github.com/thermofisherlsms/pyeds")
        print("="*80 + "\n")
except Exception:
    # Fail silently if there's an issue with the message
    pass