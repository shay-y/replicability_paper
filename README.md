Prerequisites:
--------------

(1) Install R and Rstudio:
    Follow the instruction in RStudio's website https://www.rstudio.com/products/RStudio/#Desktop

(2) Run R and install the following packages (needed only once) with the command:

    `install.packages(c('rJava','XLConnect','xtable','nlme','multcomp','knitr','dplyr','lubridate'))`
    
    rJava installation may run into errors if Java is not updated.

(3) Install MikTex by running the following:  
    `install.packages('installr');  
    library(installr);   
    install.miktex()`
    or by following the directions from: http://miktex.org/download
    
Reproducing the Analyses and the Simulation:
--------------------------------------------

(1) Clone this repository or download the .zip from from Github or from 
    http://replicability.tau.ac.il/images/Replicability_Paper.zip and extract the files to a folder named `Replicability_paper`.  

(2) The folder will serve as the working directory for R and knitr processes. Open the file `Replicability_Paper.Rproj` to begin a new RStudio project in the folder.

(4) In Rstudio, open one of the three `.Rnw `files:
      
      * `1_IMPC_datasets_analysis.Rnw` for IMPC datasets (Males and Females) analyses.
      * `2_other_datasets_analysis.Rnw` for the other six datasets analyses.
      * `3_simulation.Rnw` for the simulation run.

(5) Click "compile PDF" on the top of the open script window to run the script and compile a pdf with a report resulted from the analysis.
    The `.pdf` file will output to the same directory, along with a TeX file and a log file of the Latex compiling process.

Troubleshooting:
----------------
    
* Make sure that RStudio uses __knitr__ for waeving the document (the default). Follow the steps:
    enter 'Tools' menu in Rstudio's toolbar >> Global Options... >> Sweave >> Weave Rnw files using: >> choose 'knitr'.   
* If prompted, install the "underscore.sty" MikTeX package, using the MikTeX package manager (accessed from the start menu, in windows).

Notes:
------

* The "wahlsten datasets" document compiling may take between 10 minutes to 1 hour, depending on the system resources.

* You may run the analysis on a new dataset saved as .xlsx file with the following structure:
    * each row represent an observation of a specific strain in a specific lab and contains all endpoints
    * first column is named "lab" and contains the lab variable (TEXT, will be converted to factor). 
    * second column is named "strain" and contains the strain variable (TEXT, will be converted to factor).
    * the other columns contain measurments in various endpoints (each endpoint in each column, named in the first row, all numeric).  
    
* In datasets folder there are the eight datasets analysed in the paper.
  Foo.csv includes the numbers for set.seed() function in the simulation script (this makes the simulation fully reproducible)

* Each file contains its own instructions, commented at the beginning of the files.

Simulation Notes:
-----------------

(1) The simulation code uses set.seed() function in order to reproduce the exact final results. without providing seed to the simulation,
    it will result with slightly different results, due to different random data generations.

(2) The simulation is set for a multi lab data with 3 labs and 3 strain, and one "feature" single lab.
    You may Change the parameters. In the code: " Nlab=3 " and "Nstrain=3 " in lines 23 and 24.

(1) The simulation is set for 1000 iterations. It takes about 22 hours. To change number of iterations look for "Rounds=1000 " in line 28.


Support:
--------
* For further technical questions regarding running the code, Please contact the author of the code Iman Jaljule at imanjljule@gmail.com 