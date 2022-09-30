# Electron-phonon coupling
 Calculation of the electron-phonon coupling in metals within the simplest jelly model
 
Using given material DOS the code evaluates:

1) electron-phonon coupling factor (within the jelly model)

2) electron heat capacity

3) chemical potential

4) electron effective mass in effective-one-band approximation

The jelly model is too simplistic and *is not* recomended for scientific use and research, rather for simple estimations, illustrative and educational perpuses.

## Disclaimer

Although we endeavour to ensure that the code and results delivered are correct, no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions. We shall not be liable for any damage arising from the use of this code or its parts or any results produced with it, or from any action or decision taken as a result of using this code or any related material.

This code is distributed as is for non-commercial peaceful purposes only, such as research and education. It is explicitly prohibited to use the code, its parts, its results or any related material for military-related and other than peaceful purposes. By using this code or its materials, you agree with these terms and conditions.

## Instructions

The main file of the code is Analyzing_DOS_main.f90. Can be compiled using Make.bat file (written for Windows; can be analogously created for Unix).
To execute, requires input file INPUT_PARAMETERS.txt - the file contains explanation inside on how to set material parameters.
Also required to have a file with the electronic density of states (DOS), saved in the directory DOS. New materials may be added by analogy.

The code is not optimised and rather slaggish, not meant for efficient calculations, rather for simplistic estimations, given that the model behind it is itself very simplistic (for details see Section 4 in Ref.[1]).

Output files are self-explanatory, and each has first two lines commenting on the values and units in the printed out columns.

## How to cite

The use of the code is at your own risk. Should you choose to use it, appropriate citations are mandatory:

[1] S.A. Gorbunov, N.A. Medvedev, R.A. Rymzhanov, P.N. Terekhin, A.E. Volkov, Nuclear Instrum. Meth. B 326 (2014) 163-168
