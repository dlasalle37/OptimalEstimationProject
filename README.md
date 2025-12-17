Downloading locally should allow you to run this without any problems. Contents are described below

# Main Directory
Scripts beginning with `main_` are the scripts where the simulations are. The UKF and EKF are sims combined. GMEKF is separate (and not done yet). They are also separated by cases. Use these to generate the figures shown in the report.
  
  Some options set at top of these scripts: 
  
    `ntrials`: how many trials to run.
    
    `noiselevel`: "high" or "low". Sets which pair of process/measurement noise to use. I may or may not include this in the final project paper.
    
    `save_figs`: "true" or "false", just determines whether or not to save the generated figures.

Scripts beginning with `test_` are just things I used to make sure functions are working as intended.
Each other .m file here is a function used within the scripts.

# Extras
This directory just contains the code I used to generate the jacobians of the system model and the measurement model. The .c/.txt files in here are created by the .m scripts.

# Plots
This directory contains all the saved figures. By turning `save_figs` to true in the `main_` scripts, these will be overwritten when running.


