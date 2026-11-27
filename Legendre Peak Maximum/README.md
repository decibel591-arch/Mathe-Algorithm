# Cubic peak maximum fit

In the limit a peak maximum is quadratic about the maximum.  But it is more practical, a better approach, to fit to a cubic and solve for the maximum:  Assume the peak is asymmetric.  

<img width="567" height="409" alt="PeakMaxFit" src="https://github.com/user-attachments/assets/87fedfe8-f9f7-4041-a57f-afd074e8995f" />. 

[Download ROOT Code file](https://github.com/user-attachments/files/35517078/LgndrPkMx.txt) for peak fit.  

Uses data file from Legendre Smooth.   

Start Root in a terminal

load code

    .L LgndrPkMx.txt

command

     LgndrMax("trial.dta", 5, 143)

first argument is file with data, second argument defines the fit window and third the data centre. 

[Download Python Script](https://github.com/user-attachments/files/35520460/LgndrMax.zip) and data.

