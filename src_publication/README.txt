We recommend setting up a virtual environment with the 

- required versions of python packages (requirements.txt)
- a free academic licence for the gurobi solver



MPM_main.py
    -     main file to generate the multi-phase environment-coupled model of leaf metabolism and to perform analysis
    
plot_linker_fluxes_and_pareto_frontier.py   
        - plot linker fluxes and pareto frontier as in Fig.2
        
plot_budgets.py
        - plot budgets as in Fig.3

plot_parameter_scan_results.py
        - plot parameter scan results as in Fig.5



Parameters folder:

	model_parameters.xml
	- phase, linker and balance constraints

	pathwayprioritizer.csv
	- sorted list of metabolic pathways present in the model

	rxnpathwaydict.csv
	- matches reactions and pathways present in the model

	Temperature_and_Humidity_Curves.csv
	- parameter file for T and RH

	Temperature_dependent_WaterVapour.csv
	- parameter for gas-diffusion equation


Functions folder:

        - contains helping functions for main file



