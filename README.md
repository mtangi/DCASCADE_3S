# DCASCADE_3S
The repository contains the code, input and output data used in the paper M. Tangi, S. Bizzi, R. Schmitt, A. Castelletti "Balancing sediment connectivity and energy production via optimized reservoir sediment management strategies" 

----

The data are divided into five folders

    Folder "DCASCADE_input" contains the input used for the simulations, including the 3S River network and the hydrological datasets;
    Folder "DCASCADE_functions" contains all the general functions called by D-CASCADE;
    Folder "Step_1_Baseline_sediment_budget" contains the script and function used in Step 1: Baseline sediment budget definition, as described in the    paper;
    Folder "Step_2_Reservoir_siting" contains the script and function used in Step 2: Reservoir siting assessment, as described in the paper;
    Folder "Step_3_Reservoir_sediment_management" contains the script and function used in Step 3: Reservoir sediment management, as described in the paper.
    
---

Each of the folders relative to a research step contains;
   a single script "Main_script..." reporting the general operations necessary to obtain the output discussed in the paper;
   "Input" folders containing all the additional inputs used specifically for the step considered;
   "Function" folders reporting the functions used only in the step considered;
   "Result" folders containing a Matlab data repository file (.mat) with the raw data obtained by the analysis conducted in the step;
   "Plot" folders containing scripts to plot the figures included in the paper.
   
---

Step 3: Reservoir sediment management is divided into 2 part, with 2 dedicated folders:
    Folder "Part_1_flushing_design_optimization" contains the operations to run the multi-objective optimization of the reservoir flushing designs;
    Folder "Part_2_sensitivity_analysis" contains the operations to run a sensitivity analysis on the Pareto Optimal designs derived in Part 1.
    
---

The code to run the BORG Multi-objective evolutionary algorithm used in the paper is avaiable on request from http://borgmoea.org/

---
