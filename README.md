# ACVU_scripts

01. Run "01downsizedMovie.py"
Select experiment folder
Hatching timepoint can be set to 0 if not needed
Let the script run overnight to generate all the movie for all channels needed

02. For each worm, open the coolLED movie with imagej
In a txt file called skin.txt, annotate hatching time, and timepoints in which skin appears (L2 and L3)
Timepoints are as shown by imagej - 1 (compare skin files already present)

02. Run "02markGoand.py"
Select the _analyzedImages folder
Annotate the gonad position starting from 5 hours before L2 (time shown on top-left corner of the image)
NOTE: once clicked, the figure automatically shows the next image, mouse wheels manually goes to next/previous timepoint, right click removes the gonad position of current timepoint
NOTE: the gonad pos at the current timepoint is not saved, so you need to move to next or previous timepoint before clicking the save button

03. Run "03cropImages.py"
Select experiment folder
Let the script run overnight
NOTE: When done, the raw data can be deleted!

04 Run "04markCells.py"
Manually annotate the timepoint in which cells divide and the ID of the AC by looking at images by eye
Annotate data in "birthOrderToOutcome.txt" file (compare txt files already existing)

05 Run "10histoBirthOrderVSOutcome"
to generate figure of 1st born/2nd born = AC
NOTE: folders with experimental data are hardcoded. Change it to include more data if running more experiments.

06 Run "04markCells.py"
Label cells and background in every timepoint possible, ending at L3 molt
Save data

07 Run "09cellsMovie.py"
experimental folder is hardcoded!

08 Run "15quickAndDirtyFluoCompute.py"
experimental folder is hardcoded!

09 Run "16quickAndDirtyFluoPlot", "16quickAndDirtyRatioPlot"
to visualize data
Run "histoTimeToDecide_q&d" to look at time to decision
All these scripts are hard coded

10 To make life easier, run "CreatSingleFileWithAllData"
to generate flitered trajectories and store raw and filtered fluorescence data in a single list
NOTE: The object is not pickled anywhere, manually use pickle.dump to store it anywhere you like

11 Use "q&dAbsPlot_test", "q&dRatioHist_test", "q&dRatioPlot_test", "q&dTdec_test"
to make all the plots discussed above in a much faster way with the singleFile created at point 10
NOTE: you need to specify the location of the file to load!