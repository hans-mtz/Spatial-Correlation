# Spatial Correlation

To reproduce
    - Simulation and getting results for tables 1-6 run `make matlab`
    - Paper run `make paper`

## Bandwidth and Critical Value Selection

`mainCV.m` runs the code for bandwidth and critical value selection. The script is located in the `Theta` folder. 

To review the code online, open script with this link [mainCV.m](Theta/mainCV.m) . The functions used in `mainCV.m` are located in the [functions](Theta/functions/) folder. 

To run the `mainCV.m` script in your machine, download the `Theta~ folder. Once you are in the file, change the working directory to the `Theta` folder, and add the `functions` folder to the path. You can then run the file by sections or the whole script. You can do this in MATLAB with the following commands:

```matlab
cd Theta % Change to the Theta directory
addpath('functions'); % Add the functions directory to the path
mainCV; % Run the mainCV script
```
The script is not currently printing output to the console.