# matlabcodes

This repo conatins the matlab code of the final year project we worked on.The code here is for the modules of automated cell counting and 
tumor growth modelling based on the probability of contracting cancer and number of days for simulaitng cancer cell growth.

tumorui.m - provides an ui for doctor where he can key in the various parameters required for simulation and he can choose to visualise
a 2-D or 3-D version of tumor growth .

cellcountui.m - provides an UI for doctor to upload the growth model for automated cell counting.Also contains the
function that does the cell counting.We use various binarization techniques and weiner filter. Finally we use watershed segmentation which 
consider gray regions as altitude of relief.
