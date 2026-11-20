# Sliding Legendre window

Uses a sliding window of a Legendre polynomial fit.  


Start Root in a terminal

load code

    .L LgndrSmth.txt

issue command

     LgndrSld("trial.dta", 3, 5)

first argument is file (for the example trial.dta) with data, second is maximum polynomial degree (0lg -1) and last the size of smoothing window (2*wnd-1).  

<img width="555" height="400" alt="LegSmth" src="https://github.com/user-attachments/assets/ce2b481e-fe91-4a1a-942c-59be83d42366" />

The Figure show that although the smoothing window is large the sharpness of feaures is retained.

<img width="555" height="400" alt="TopHat" src="https://github.com/user-attachments/assets/8baec037-1981-4934-aa79-43de64466e41" />

This is not the case for a sliding Topcoat: The peaks start to broaden and reduce in intensity. 


