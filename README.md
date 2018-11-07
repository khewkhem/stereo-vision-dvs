# stereo-vision-dvs

* C++ implementation of a belief-propagation event-based stereo matching algorithm (originally implemented in MATLAB [EMP](https://github.com/harryxz/EMP)). 
* program takes in as input a specified csv file (file location prompted by user, e.g. 'data/box1.csv') containing DVS events and outputs a csv file containing stereo-matched events with associated disparity.

## Getting Started

### Prerequisites

* compilation depends on the Armadillo library, which can be installed by following instructions on  [Getting started with Armadillo a C++ Linear Algebra Library on Windows, Mac and Linux](https://solarianprogrammer.com/2017/03/24/getting-started-armadillo-cpp-linear-algebra-windows-mac-linux/)

### Installing

* enter 'make process_csv' to compile the program

## Running the tests

* run main program with './process_csv'
* when prompted enter the name of the csv file containing events (e.g. 'data/box1.csv', 'data/walking2.csv')
* then when prompted enter the desired name of output csv file containing stereo-matched events (e.g. 'box1_out.csv', 'walking2_out.csv')

## Built With

* [Armadillo](http://arma.sourceforge.net/) - Linear algebra library

## Authors

* **Seth Siriya** - *Contributing author* - [khewkhem](https://github.com/khewkhem)
* **Alexander Matthies** - *Contributing author*
* **Nico Hertel** - *Contributing author*

## Acknowledgments

* Thanks to Lukas Everding for supervising the project, and Matthias Emde for assisting with hardware setup.
* Also thanks to authors of original MATLAB implementation of algorithm [EMP](https://github.com/harryxz/EMP), which this C++ implementation is based on, and for providing event data.
