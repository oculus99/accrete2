# accrete2

This code is modified classical Dole accretion model.
Under development.
Code has some command line parameters, csv file output end experimental Pov-ray file output.


Compilation on Linux :

gcc accreteb.c -lm -o accreteb


<code>
  Program name: ./accreteb

Usage:
  ./accreteb [-mass <value>] [-migrate] [-coeff <value>] [-ecosphere <value>]
  -mass <value>      : Stellar mass ratio (double, default 1.0)
  -dust_density_coeff <value>      :  dust amount factor. Doles default 0.002 
  -gasdust <value>      : gas per dust, default 50. Sun maybe 70 
  -alpha <value>      : dust exp(-a) density by distance modifier , default 5
  -nanna <value>    : dust exp(-a) density by distance modifier , default 3
  -innermost_planet <value>      : innermost planet, default 0.3 AU
  -disk_radius <value>      :  dust disk radius, default 50 AU
  -migrate           : Enable planetary migration (no argument needed)
  -coeff <value>     : Migration coefficient (double, default 0.1)
  -ecosphere <value> : Ecosphere radius (double, default 1.0)
  -h, --help         : Display this help message

</code>


To render output POV-ray file, you must have installed Pov-Ray and povray path in your system path.
