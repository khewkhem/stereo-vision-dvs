# stereo-vision-dvs
What does the program do?
- Takes in as input a specified csv file (file location prompted by user, e.g. 'data/box1.csv') and outputs a csv file containing stereo-matched events with disparity.

Compiling programs
 - enter 'make process_csv' to compile the program

Sidenote
 - compilation depends on the Armadillo library, which can be installed by following instructions on the following website
 - (https://solarianprogrammer.com/2017/03/24/getting-started-armadillo-cpp-linear-algebra-windows-mac-linux/)

How to run program?
 - run main program with './process_csv'
 - when prompted enter the name of the csv file containing events (e.g. 'data/box1.csv', 'data/rotate.csv')
 - then when prompted enter the desired name of output csv file containing stereo-matched events (e.g. 'box1_out.csv', 'rotate_out.csv')
