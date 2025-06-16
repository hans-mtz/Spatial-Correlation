# Spatial Correlation

To reproduce
    - Simulation and getting results for tables 1-6 run `make matlab`
    - Paper run `make paper`

## Bandwidth and Critical Value Selection

To review the code, open ![mainCV.m](Theta/mainCV.m)

To run the `mainCV.m` script, you need to change the working directory to the `Theta` folder, and add the `functions` folder to the path. You can do this in MATLAB with the following commands:

```matlab
cd('Theta'); % Change to the Theta directory
addpath('functions'); % Add the functions directory to the path
mainCV; % Run the mainCV script
```
