Before running the code, make sure you have all the neccessary python scripts in the same folder together with the txt files with the data you want to analyse. Then to run the code:

1) Edit BEHR_DIR in HR_driver_const_counts.py with the directory path to where you have BEHR installed. If you don't have it yet, you can install it from: http://hea-www.harvard.edu/astrostat/behr/
2) Open the terminal
3) Run ciao
4) Go to the directory with your data txt files and code scripts.
5) In terminal write "python PlottingDataWithHR.py binsize N", where binsize = size of bins in seconds for lightcurve plot, N - photons to the right and left to consider for HR plot. Eg. "python PlottingDataWithHR.py 1000 8"
Extras) Can also change the divide energy between soft and hard photons by writing "python PlottingDataWithHR.py binsize N divide". Eg. "python PlottingDataWithHR.py 1000 8 2100"

This code outputs png images with plots and possibly a "log.txt" file where all the errors while running BEHR are stored. If an error occurs running BEHR, then only 3 plots without HR will be plotted and stored in png image.
