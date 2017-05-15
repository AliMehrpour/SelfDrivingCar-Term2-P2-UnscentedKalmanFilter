# Unscented Kalman Filter
This Project is second project of second term of Udacity Self-Driving Car Nanodegree program. The goal is to apply Unscented Kalman Filter to data from LIDAR and RADAR sensors by C++.

## Content
- `scr` contains source code:
  - `main.cpp` - loads data, calls a function to run the Unscented Kalamn Filter(UKF), calls a function to calculate RMSE
  - `data_source.cpp` - load input data into measurement and ground truth vectors
  - `ukf.cpp` - initializes the filter, calls the predict function, calls the update function
  - `tools.cpp` - a function to calculate RMSE and the Jacobian matrix
- `data` a directory with input files, provided by Udacity
- `result` a directory with output files, log files

## Results
In order to see how sensor data impact accuracy of Kalman Filter, I ran the algorithm on two data input (porvided by Udacity). Here are the RMSEs comparion of expected and above calculations:

### Data 1
|   RMSE    | px | py | vx | vy |
| --------- | -- | -- | -- | -- |
| Threshold | 0.09 | 0.09 | 0.65 | 0.65 |
| Radar + Lidar | 0.0772625 | 0.0817972 | 0.589278 | 0.574289 |
| Radar | 0.104241 | 0.111766 | 0.57746 | 0.572511 |
| Lidar | 0.0924323 | 0.0831921 | 0.650915 | 0.608865 |


### Data 2
|   RMSE    | px | py | vx | vy |
| --------- | -- | -- | -- | -- |
| Threshold | 0.20 | 0.20 | 0.55 | 0.55 |
| Radar + Lidar | 0.19322 | 0.189588 | 0.300693 | 0.435957 |
| Radar | 0.249102 | 0.707066 | 1.63246 | 0.932646 |
| Laser | 0.319107 | 2.76411 | 4.85391 | 5.10034 |


## Lessons Learned
- Considering both sensor data (LIDAR and RADAR) give much better accuracy. You can see in RSME table that when we are considering both sensor data, we have better and acceptable RMSE.
- Lidar sensor give more accurate measurement than radar
- UKF get better accuracy than EKF.

## Build Instructions
1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake ../src && make` 
   * On windows, you may need to run: `cmake ../src -G "Unix Makefiles" && make`
4. Run it: `./UKF path/to/input.txt path/to/output.txt true|false true|false`. You can find
   some sample inputs in 'data/'.
    - eg. `./UKF ../data/sample-laser-radar-measurement-data-1.txt output.txt true|false true|false`
