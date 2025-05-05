# cdSirius
**cdSirius** is an implementation of the [Sirius suite of computational mass spectrometry tools](https://bio.informatik.uni-jena.de/software/sirius/) maintained by the Böcker lab at University of Jena within [Compound Discoverer v3.3](https://www.thermofisher.com/ch/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/compound-discoverer-software.html) from ThermoFisher Scientific.  This software is intended to allow annotation of metabolites and other small molecules with molecular formulas, 2D structures, and compound classes based on acquisition of high-resolution, accurate mass MS/MS through data-dependent LC-MS/MS analysis of complex samples.  Note that cdSirius will only work with high-resolution MS and MS/MS data, and has only been tested with data-dependent analysis (DDA) results.  It may work with DIA data, but no guarantees are given.

Implementation of Sirius within Compound Discoverer uses the scripting node functionality within CD as an extensible interface between independent software programs.  Within this implementation, Sirius runs on the host PC and is called as a background service during post-processing through the Sirius API.  Mass spectral data for discrete compounds detected by Compound Discoverer are passed to Sirius for processing, and results are reported back to CD and persisted to the result file for viewing and interpretation.  Linkages among CD compounds and Sirius results are maintained within the resulting tables.

When using cdSirius results in a publication, please be sure to cite the work that enabled creation of this resource.  Visit the Sirius development group's site referenced above for detailed citation information, and use the primary citation as follows:

Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Alexander A. Aksenov, Alexey V. Melnik, Marvin Meusel, Pieter C. Dorrestein, Juho Rousu and Sebastian Böcker. [SIRIUS 4: Turning tandem mass spectra into metabolite structure information](https://doi.org/10.1038/s41592-019-0344-8). _Nature Methods_ 16, 299–302, 2019.
## Dependencies
Installation of cdSirius requires a number of dependencies to be installed on the host PC.  _Important_: It is recommended to install all of these dependencies as administrator so that all users will have access.  This is particularly important for the python packages.  Also, the packages should be installed in the base python environment.

> **⚠️ IMPORTANT:** cdSirius requires the following packages which must be installed separately from GitHub:
> 
> 1. **pySirius** - install with:
> ```
> pip install git+https://github.com/sirius-ms/sirius-client-openAPI@v0.0.0#subdirectory=client-api_python/generated
> ```
>
> 2. **pyEDS** - install with:
> ```
> pip install git+https://github.com/thermofisherlsms/pyeds
> ```
> 
> These packages are not available on PyPI and must be installed before using cdSirius.

Other required dependencies:
* Fully licensed installation of Compound Discoverer v3.3, SP 3
* [Sirius](https://v6.docs.sirius-ms.io/install/) minimum program version 6.0.1
* [Python v.3.11](https://www.python.org/downloads/release/python-3110/).  Note that more recent versions of Python may work but have not been tested.  _Note:_ it is recommended to install Python at `C:/Program Files/Python311/python.exe`.
* [RDKit](https://pypi.org/project/rdkit-pypi/) cheminformatics python package (for implementation of future functionality in cdSirius)
* [pandas](https://pypi.org/project/pandas/) data science python package
* [molmass](https://github.com/cgohlke/molmass) python package for molecular formula manipulation
## Sirius user account
To run Sirius with CSI:FingerID and CANOPUS functionality, you will need a Sirius user account.  You can create a user account as described in the [Sirius Wiki](https://v6.docs.sirius-ms.io/account-and-license/).  The username and password will be used in the cdSirius node for authentication.  Do not re-use a sensitive password for this user account, as the password will not be encrypted and will be visible in plain text within the CD method editor. _Note_: If you are an academic user, you will qualify for a free account, but you must use your institutional email address when you register your account.  
## Installation

cdSirius is available on PyPI and can be installed using pip. This method makes it easier to manage dependencies and updates.

1. Ensure you have Python 3.11 installed. Recommended path is `C:/Program Files/Python311/python.exe`

2. Install the required GitHub dependencies first:
   ```
   pip install git+https://github.com/sirius-ms/sirius-client-openAPI@v0.0.0#subdirectory=client-api_python/generated
   pip install git+https://github.com/thermofisherlsms/pyeds
   ```

3. Install the package from PyPI:
   ```
   pip install cdsirius
   ```

   This will install cdSirius and other required dependencies:
   - pandas
   - rdkit-pypi
   - molmass
   - numpy
   - matplotlib
   - scipy
   - openpyxl
   - requests

4. Set up Compound Discoverer integration using the included CLI tool:
   ```
   cdsirius-setup "C:\Program Files\Thermo\Compound Discoverer 3.3" "C:\Program Files\sirius\sirius.exe"
   ```
   
   This tool will:
   - Create the necessary directory structure in your Compound Discoverer installation
   - Copy required files to the Compound Discoverer Scripts folder
   - Update paths in node.json to point to your Python installation
   - Create a bootstrap script to connect Compound Discoverer to the installed package

   The paths to your Compound Discoverer and Sirius installations are required. You can also specify:
   - The Python executable path:
     ```
     cdsirius-setup "C:\Program Files\Thermo\Compound Discoverer 3.3" "C:\Program Files\sirius\sirius.exe" --python-path "C:\path\to\python.exe"
     ```

5. Complete the installation:
   - Launch Compound Discoverer 3.3 SP3
   - Navigate to Help -> License Manager
   - Run "Scan for Missing Features"
   - Restart Compound Discoverer

## Using cdSirius within a Compound Discoverer workflow
The cdSirius node is a post-processing node that can be appended to an existing full processing workflow, or it can be included in a "reprocessing" workflow to retrospectively add Sirius results to the cdResult file.  Either way, you will find the new Sirius node within the Workflow Editor Node menu, in the _10. Post-Processing_ sub-menu:
   
   <img width="194" alt="image" src="https://github.com/user-attachments/assets/c5a6a15b-ec71-4182-91bf-ca05b04d28fd" />

After adding the node to the workflow, the processing configuration dialogue is available for editing.  Default parameters are provided for most settings:

   ![image](https://github.com/user-attachments/assets/30da7f94-f87f-4e95-aaae-25222286860a)

   **Figure 3.** Sirius node parameter configuration.

### cdSirius parameter settings
The range of possible settings for Sirius is very large and the corresponding job configurations can become quite complicated.  The settings available within cdSirius represent a subset of possible parameters, chosen based on their general applicability and typical use cases.  A complete guide for Sirius job parameters is beyond the scope of this program, but extensive documentation is available for Sirius [elsewhere](https://v6.docs.sirius-ms.io/methods-background/).
1.  **Sirius Program Settings:**  These are global settings for the Sirius program service.
   - <ins>Sirius Program Path</ins>: Here you can set the program path for the Sirius executable.  The default should be suitable for most installations, but can be changed for non-standard installations.
   - <ins>Save Sirius Result</ins>: Setting this to `True` will enable the .sirius workspace to be persisted as a permanent file, saved to the same directory with the cdResult file.  This is useful if e.g. you plan to re-open and analyze data using the Sirius GUI at a later time.
   - <ins>Save Sirius Predicted Fingerprints</ins>: This setting will toggle saving tab-separated ASCII files (.txt) containing the predicted fingerprints for compounds processed in Sirius as well as the fingerprint definition key used by Sirius.  These files will be saved to the same directory as the cdResult file.
   - <ins>Sirius Username</ins>: Your username as configured on the BrightGiant web portal as described above
   - <ins>Sirius Password</ins>: Your password as configured on the BrightGiant web portal as described above (**Note**: Do not re-use sensitive passwords here as this information is not encrypted or hidden)
2.  **Compound Selection Settings:**  These settings modify which Compounds within the CD processing environment are selected for submission to Sirius
   - <ins>Checked Feature Status Handling</ins>: This is a switch that allows for down-selection of only compounds of interest when cdSirius is used in reprocessing of existing cdResult files.  Selecting "Checked" will pass only "checked" compounds to Sirius for calculation. The default of "All" will pass all non-background compounds (which have not been pre-filtered within the workflow by peak quality thresholds) to the Sirius service.
   - <ins>Peak Quality Threshold</ins>: This threshold allows a more stringent selection of chromatographic peaks for Sirius processing.  See the CD documentation for an explanation of peak quality metrics.
   - <ins>Maximum MW</ins>: Sirius computation becomes extremely slow for large molecules, so it is best to limit the upper MW range to only those of interest for a particular analysis.
3.  **Molecular Formula Prediction Settings:**  This set of parameters controls the "base" molecular formula calculation functions within Sirius and is necessary for all further processing
   - <ins>Maximum Formula Candidates</ins>: Adjust to increase or decrease the allowable formula candidates that can be considered by Sirius.
   - <ins>MS1 Mass Accuracy [ppm]</ins>: The known mass accuracy threshold (in ppm) for your MS1 data.  _This is a critical parameter and must be set accordingly.  If your instrument is equipped with EasyIC, it is suggested that you use it_
   - <ins>MS2 Mass Accuracy [ppm]</ins>: The known mass accuracy threshold (in ppm) for your MS2 data.  _This is also critical and is often a less accurate measure than for MS1, especially when using lower resolutions (e.g. 15K).  It is not recommended to use EasyIC for Orbitrap MS2 with resolutions < 60K_
   - <ins>Filter by Isotope Pattern</ins>: Enabling this parameter will use isotope pattern measurement as a pre-filter to exclude formulas that are inconsistent with the measured isotope pattern, regardless of MS/MS tree score.
   - <ins>Enforce Lipid Detection Filtering</ins>: This setting enables an internal Sirius algorithm that attempts to detect fragmentation patterns characteristic of lipids.  When detected, the corresponding lipid-like molecular formula will be prioritized as a candidate.
   - <ins>Perform Bottom-Up Formula Search</ins>: This setting allows for the use of a "bottom-up" formula candidate selection strategy, which uses combinations of known formulas for potential sub-fragments to build candidate molecular formulas.  It is less restrictive than searching a database for candidate formulas but is also less computationally-intensive than a true _de novo_ formula prediction strategy.
   - <ins>_De novo_ Formula Generation Threshold</ins>: Below this _m/z_, all molecular formulas are calculated using a _de novo_ approach, which maximizes the chance to observe novel formulas (which are not present in any databases).  Above this _m/z_, formula candidates are predicted using the "bottom-up" strategy if enabled above.
   - <ins>Formula Elemental Constraints</ins>: Use this string to specify which elements should be considered for _de novo_ formula prediction.  Numbers in brackets represent maximum possible element counts.  Elements without numbers are given unlimited maximum counts.  **Note:** Do not include B, Cl, Br, S, or Se in this list, as those elements are detected automatically using observed isotope patterns in the MS1 spectra.
4.  **Structure Prediction Settings:**  These settings control the CSI:FingerID structure prediction toolset within Sirius.
   - <ins>Predict Structures</ins>: Toggle enabling CSI:FingerID database search
   - <ins>PubChem as Fallback</ins>: cdSirius uses the US EPA DSSTox database as its default molecular structure database for searching compound structure candidates.  When this parameter is set to "True", Sirius will search PubChem for structure candidates in the event that no viable structure candidates were found within the target database.  In future releases of cdSirius, database choices beyond DSSTox will be available.
5.  **Compound Class Prediction Settings:**  This parameter set controls the CANOPUS classification algorithm.
   - <ins>Predict Compound Classes</ins>: When this parameter is set to "True", the CANOPUS implementation of the ClassyFire algorithm is used to predict compound classes from molecular fingerprints.
6.  **De Novo Structure Prediction Settings:**  These settings control the MSNovelist toolset for database-independent molecular structure prediction.  **Caution:** MSNovelist is quite computationally intensive.  Be careful when using this toolset with full (unfiltered) compound sets from Compound Discoverer, as the job compute times can become extremely long.
   - <ins>Predict de Novo Structures</ins>: Enables or disables MSNovelist processing.
   - <ins>De Novo Structure Candidates Limit</ins>: This setting throttles the number of possible _de Novo_ structure candidates considered for each molecular formula, for each Compound
7.  **General:**
   - <ins>Archive Datafiles</ins>: If set to True, the temporary table files and JSON-based response file produced by the cdSirius node for ingestion by Compound Discoverer are persisted in a folder within the same directory as the corresponding cdResult file.


