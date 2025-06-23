# Spatial Correlation

To reproduce:

- Simulations, run `make matlab`
- Report of simulations, run `make report`

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