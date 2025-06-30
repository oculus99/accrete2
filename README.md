# accrete2

This code is modified classical Dole accretion model.
Under development.
Code has some command line parameters, csv file output end experimental Pov-ray file output.


Compilation on Linux :

gcc accreteb.c -lm -o accreteb


Usage:
  ./accreteb [-mass <value>] [-migrate] [-coeff <value>] [-ecosphere <value>]
  -mass <value>      : Stellar mass ratio (double, default 1.0)
  -migrate           : Enable planetary migration (no argument needed)
  -coeff <value>     : Migration coefficient (double, default 0.1)
  -ecosphere <value> : Ecosphere radius (double, default 1.0)
  -h, --help         : Display this help message


To render output POV-ray file, you must have installed Pov-Ray and povray path in your system path.
