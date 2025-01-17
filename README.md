# stoichiometry_analysis
Original author: Yuri Quintana (yuriqp@gmail.com), 2016.

Working author: Jonas Aufdermauer (jonas.aufdermauer@gmail.com), since 2022.

Version status: v1.0

How to run/install (Windows 11):

Option 1: Run dist/analyze_stoichiometry.exe (no installation necessary)

Option 2: Manual installation:

1. Install all required packages with pip install -r requirements.txt
2. Run python stoichiometry_main.py


How to use the script for batch analysis:

1. Select a folder for the calibration (eg. test_dataset/cal_nup96)
2. Select a folder for the data (eg. test_dataset/data).
3. To batch analyze all biological replicates at once, checkmark analyze subfolders. 
4. For the test dataset, the threshold should be 3e-5.
5. Press analyze. This may take a while for batch analysis (10-20min)
6. After the analysis, the analyzed images can be observed in the Particle Detection tab.
7. For the numbers (stoichiometry, number of particles kinetics and so on), the Brightness tab features an Export as csv button.


If single replicates should be analyzed:

1. Select a folder for the calibration (eg. test_dataset/cal_nup96)
2. Uncheck analyze subfolders.
3. Select a folder for the data (eg. test_dataset/data/1).
4. Go to the Particle Detection tab (sometimes you have to go twice back and forth) and there you can press detect, then localize and then several constraints can be applied on the particles like distance, threshold and stacking.
5. Go to Brightness tab and there you can either export first frame intensities, monomer intensities or first frame intensities in separate lines (they are the raw uncalibrated brightness, that divided by the calibration eg. nup96 intensity, yields the stoichiometry via the ratiometric approach.