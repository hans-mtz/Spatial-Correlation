# Spatial Correlation

To reproduce:

- Simulations, run `make matlab`
- Report of simulations, run `make report`

## Monte Carlo Simulations with Different Rho Values

The makefile includes targets for running Monte Carlo simulations (`mainMC.m`) with different values of the spatial correlation parameter rho. Here are the available commands:

### Run All Rho Values
Run Monte Carlo simulations for all predefined rho values (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
```bash
make matlab-rho
```

### Run Specific Rho Value
Run Monte Carlo simulation for a specific rho value:
```bash
make matlab-rho-single RHO=0.3
make matlab-rho-single RHO=0.7
```

### Run Individual Simulations by Log File
You can also run simulations by targeting specific log files:
```bash
make Theta/mainMC-r0.3.log
make Theta/mainMC-r0.7.log
```

### Clean Up Log Files
To remove all rho-specific log files:
```bash
make clean-ml
```

The log files will be created in the `Theta/` directory with names like `mainMC-r0.3.log`, `mainMC-r0.7.log`, etc., making it easy to identify which simulation corresponds to which rho value.

## Bandwidth and Critical Value Selection

`mainCV.m` runs the code for bandwidth and critical value selection. The script is located in the `Theta` folder. 

To review the code online, open script with this link [mainCV.m](Theta/mainCV.m) . The functions used in `mainCV.m` are located in the [functions](Theta/functions/) folder. 

To run the `mainCV.m` script in your machine, download the `Theta` folder. Once you are in the file, change the working directory to the `Theta` folder, and add the `functions` folder to the path. You can then run the file by sections from MATLAB or the whole script. You can do this in MATLAB with the following commands:

```matlab
cd Theta % Change to the Theta directory
addpath('functions'); % Add the functions directory to the path
mainCV; % Run the mainCV script
```
The script will print to the console, the critical value for a choice of bandwidth, number of PCs of the tensor product of a 10x10 B triangular B-Splines, rho, and theta, simulating 100 times. 

The script will also print the power of the test for a value of beta candidate.