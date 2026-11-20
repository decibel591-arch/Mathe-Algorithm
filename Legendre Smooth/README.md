# Sliding Legendre window

Uses a sliding window of a Legendre polynomial fit.

Start Root
load code
  .L LgndrSmth
issue command
   LgndrSld("trial.dta", 3, 5)
first argument is file with data second is maximum polynomial degree (0lg -1) last the size of smoothing window (2*wnd-1)

