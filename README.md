# WAVE

(By Ting Chen, [Fengnan Gao](https://gaofn.xyz/ "Fengnan's Homepage") and [Yan-Wen Tan](https://phys.fudan.edu.cn/f7/50/c7605a63312/page.htm "Yan-Wen's faculty page"))

Implementation of WAVE (**W**asserstein distance **A**nalysis in steady-state **V**ariations in smFR**E**T) to detect and locate non-equilibrium transition positions in FRET trajectories, which is first introduced in Chen, Gao and Tan (2023).[^1]

[^1]: Chen, T., Gao, F. and Tan, Y-W. (2023) Transition Time Determination of Single-Molecule FRET Trajectories via Wasserstein Distance Analysis in Steady-state Variations in smFRET (WAVE).

## Main script description

1. **`TestNonequilibrium3`**: Data analysis script, the core of WAVE.

2. **`MakeNonequilibriumHMMPoissonData`**: This script generates simulated data (FRET trajectories) via HMM.

3. **`NonequilibriumCompare`**: This script evaluates the analysis result after `TestNonequilibrium3` has finished processing the simulated data generated by `MakeNonequilibriumHMMPoissonData`.

## User guide

### 1. Run `TestNonequilibrium3` individually

Run `TestNonequilibrium3` script in Matlab, and select a folder with structure as follows:

***
**Folder_name** (a folder containing all FRET trajectories)  

-->**XXX.txt** (records intensity-time trajectory, containing 2 columns data with same length. The first column is the intensity-time trajectory of donor channel, and the second column is the intensity-time trajectory of acceptor channel)  
------**E** (a folder, contain FRET efficiency and region information of all trajectories in FRET folder)  
--------->**XXX.txt Efficiency.txt** (corresponding to XXX.txt, records FRET efficiency-time  trajectory in the first column, with the same length as intensity-time trajectory. In FRET region: [FRETbegin, FRETend], the recorded FRET efficiency~=0, while in crosstalk region [Crosstalkbegin, Crosstalkend] and background region [Backgroundbegin, Backgroundend], the recorded FRET efficiency=0)  
--------->**XXX.txt Region.txt** (corresponding to XXX.txt, records region boundaries. This file contains 6 elements, like follows [FRETbegin, FRETend, Crosstalkbegin, Crosstalkend, Backgroundbegin, Backgroundend])
***

We advise the user study the `test` folder for this specific structure of organizing data.

**Note**: Set carefully the parameters in `TestNonequilibrium3` before running, especially `Changeframe`, which represents the time point at which the conditions change, depending on the experimental design.

### 2. Test `TestNonequilibrium3` with simulated data

Follow the steps below:

1. Create an empty folder `new_folder` to save simulated data.  
2. Run `MakeNonequilibriumHMMPoissonData` with the newly created empty folder `new_folder`.  
3. Run `TestNonequilibrium3` with the folder `new_folder` containing simulated data.  
4. Run `NonequilibriumCompare` to evaluate the analysis result.

**Note**: The evaluation criteria in `NonequilibriumCompare` should match the kinetic rate coefficients set in `MakeNonequilibriumHMMPoissonData`. Read the comments in these scripts for details.

## Additional information

### Output of `TestNonequilibrium3`

**BNESTAnanysis** is the output of `TestNonequilibrium3`. It is a `cell` table, with the following format：

- 1st row: record the name of each FRET trajectory (.txt file)
- 2nd row: record the length of FRET trajectory after BNEST
- 3rd row: record the last state before BNEST
- 4th row: record the first state after BNEST
- 5th row: record the trajectory length after condition change and before BNEST
- 6th row: record the Length of last segment before BNEST
- 7th row: record the FRET trajectory after BNEST
- 8th row: record the fitting trajectory after BNEST

### Function classification

1. **STaSI part:** `MDLMulit.m`, `StategroupingMulti.m`, `change_point_detection`.
2. **HMM part:** `HMM_Algorithm2.m`, `HMM_Core.m`, `HMM_Data_Preparation.m`,
`HMM_FitTrace_Get_Distribution.m`, `HMM_MainSimple.m`, `HMM_Max_Likelihood.m`
`Tcalculation.m`
3. **MWD  analysis part:** `CountMatrixCalculation.m`, `FindMaximumWassersteinDistance3.m`
`FindMaximumWassersteinDistanceHMMSequence.m`, `TcalculationC.m`, `FindMaximumWassersteinDistanceOutOrder.m`
4. **Forward/Backward algorithm part:** `FindBestTracePart2.m`, `MDLCompare.m`
5. **Data simulation and evaluation part:** `MakeHMMPoissonData.m`, `postcalculationTimeBin.m`, `precalculationTimeBin.m`
