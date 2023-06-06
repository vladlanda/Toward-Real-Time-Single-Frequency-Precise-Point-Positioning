

## Genrating TEC RMS maps

1. Use "DMD_RMS_TEC.ipynb" to download the relevant ionex files, for various agencies (WHU) into "ionex..." folders
3. Run "DMD_RMS_TEC.ipync" file in order to generate DMD RMS TEC predictions (at code section :"RUNNING" change the "
days_type variable in order to utilize different solar activities periuds


## Generating NEU errors
1. Install "gLab" on you pc and make sure you can run it from terminal
2. Use "gLab_executer.ipynb" in order to execute the possition estimation using gLAB tool
3. To do so you will have to place the nessesary file, such as sattelite observation(.14o), orbits(.sp3) and clocks(.clk_30s) 
corsponding to a choosen date and place it inside ./gLAB/files/disturbance (or quiet)/<date>
###see 27_10_2014 date example.
4. The generated files will be in the output directory. 
5.Use "gLab_output.ipynb" file inorder to generate the NEU graphs.
6.You can use the "copy_tec_to_glab_folder.ipynb" script, to help you move the generated filse into gLab folder.