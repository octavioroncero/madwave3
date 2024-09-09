In this example we study the process:
$\text{H} + \text{DH}(v=0, j=0) \rightarrow \text{HD}(v', j') + \text{H}$
## Installation
Make sure you have compiled MadWave3. For that, go to the `BIN` directory in the *madwave3* repository and run `make`.

## Compute PES
MadWave3 requires to represent the PES into a grid in Jacobi coordinates.

1. Compile pot.out:
	- Inspect the `Makefile_pot` file. It can be found in the path `{madwave3}/BIN/` and should, in general, be copied from that directory.
	- Add the path to madwave3's `SRC` folder in the variable `SRC_DIR`.
	- Include all the files of your PES in the `POT_DIR` and `POT` variables. (These variables are defined for the present example).
	- Run `make -f Makefile_pot` (in case of some error run `make -f Makefile_pot clean` for a clean compilation).
2.  Inspect the `input.dat` file which defines the grid.
```
&INPUTGRIDBASE
 NPUN1=    128, RMIS1=     0.1d0     , RFIN1=     9.d0     , NPUN1MIN= 64,
 NPUN2=        128, RMIS2= 0.1d0 RFIN2=  9.d0    ,
 NANGU= 60, NANGPLOT= 2,
 NELECMAX= 1,
 JTOT=     0, IPARITY=  1
 IOMMIN=   0, IOMMAX=   0,
 J0=       0, INC= 1
 JINI=     0, JMAX=     39,
 NVINI=     0, NVMAX=     5,
 NVREF=     0, JREF=      0, IOMREF=    0, IELECREF=  1,
 /
```
Use 256 grid points in each radial coordinate `r, R` and 80 angular points.
Use a single electronic state.
```
!-----------> masses(in amu) and potential(in eV) cut data
system='H+DH'
xm0=2.014101779d0
xm1=1.007825035d0
xm2=1.007825035d0
VcutmaxeV=2.5d0
radcutmaxeV=2.5d0
rotcutmaxeV=5.d0
/
```
Specify the masses of each of the atoms. It also sets the order of the atoms (atom 1 will be deuterium and 2 and 3 protium).

3. Create the folders `pot` and `func`.
4. Execute `mpiexec -n 1 ./pot.out`. Note this code can only be called with one processor.
The output is written into the `sal.000` file. After the print of the namelist we get a list of the angle grid being used:
```
          Angular grid to plot pes and wvp every   1
           1   177.72256551502545     
           2   174.77233758453366     
           3   171.80468702752137     
           4   168.83310004270163     
           5   165.86001997818192     
           6   162.88621751861174     
           7   159.91201146051998     
           8   156.93755731164626     
           9   153.96293990907643     
           ...
```
Next, the rotational states of the reactants AB fragment  are computed:
```
  ** reactant AB^j(k) (01 fragments) wvf energies expanded in all elec. states **

 ---------> Ielectronic =            1
            Be (eV) =    5.6650653314570679E-003   r_e(Angs.)=   0.74112056028013995     
            Vmin(eV)=   -6.3049003333886253E-006     Vinfty(eV)=    4.7482552217269580     

   ** j=            0  **
     0     0    0.2343133E+00    0.2343351E+00    0.1000000E+01    0.1182163E+00
     1     1    0.6845299E+00    0.6846296E+00    0.1000000E+01    0.3545251E+00
     2     2    0.1112677E+01    0.1112914E+01    0.9999999E+00    0.5904677E+00
     3     3    0.1519214E+01    0.1519632E+01    0.9999996E+00    0.8259493E+00
     4     4    0.1904469E+01    0.1905094E+01    0.9999994E+00    0.1061014E+01
     5     5    0.2268630E+01    0.2269477E+01    0.1000005E+01    0.1295864E+01
   ** j=            1  **
   ...
```
since `iprod=2` the rotational states of the products are also computed:
```
  ** product BC^j(k) (02 fragments) wvf energies expanded in all elec. states **

   ** j=            0  **
     0    0.2343351E+00    0.1000000E+01
     1    0.6846296E+00    0.1000000E+01
     2    0.1112915E+01    0.1000000E+01
     3    0.1519632E+01    0.1000000E+01
     4    0.1905094E+01    0.1000000E+01
     5    0.2269468E+01    0.1000000E+01
     6    0.2612774E+01    0.1000000E+01
     7    0.2934917E+01    0.1000000E+01
     8    0.3235670E+01    0.1000000E+01
     9    0.3514481E+01    0.1000000E+01
    10    0.3770190E+01    0.1000000E+01
   ** j=            1  **
   ...
```
Next, the potential energy surface is computed for each of the grid's angle:
```

          only in idproc=0
        substracting Eref(vref)=   0.22631680258674705
 iang=           1
  writting file = pot/pot.001.dat
 iang=           2
  writting file = pot/pot.002.dat
  ...
```
and the number of grid points per angle:
```
  No. of points per angle
           1        4389
           2        4438
           3        4575
           4        4848
```
Finally, the minimum of the potential energy surface is displayed:
```
          Vmax (eV)=    2.5000000000000000       while Vcutmax(eV)=    2.5000000000000000     

          Vmin (eV)=  -0.23705248565602827

       at r1=rpeq (angstroms) =    3.0433070866141736
          r2=Rgran (angstroms)=    2.7629921259842525
          gam (degrees)       =    2.2774344849745574

     ** end potential calculation **
```

## Compute reaction probabilities (for J=0)
We will compute the reaction probabilities for the reactive process:
$\text{H} + \text{DH}(v=0, j=0) \rightarrow \text{HD}(v', j') + \text{H}$

1. In the calculations for total angular momentum J=0 we start creating a folder named `J000`. This folder should be at the same level than `pot` or `func`.
2. Copy the `input.dat` file into this folder.
3. Execute MadWave3: `mpirun -np 8 /path/to/madwave/mad3.out`. Note that size of the angular grid should be divisible by the number of processes (8 in this example).
After the execution the following files are generated:
- sal.{iproc}: Output files. The file sal.000 is the main output file.
- gaussE: Initial wavepacket in the energy domain. First column is the energy and second the value of the wavefunction.
- gaussR: Initial wavepacket in the position domain. First column is the $R$ coordinate and second the value of the wavefunction.
- distriS2reac.elec01: Total inelastic propability. Represents energy vs probability.
- distriS2reac.v{v'}.e1: Inelastic probabilities for the reactants in the vibrational state v'. The first column is the energy and the subsequent columns are the probabilities for increasing values of $j'=0, 1, \ldots$
- S2prod.v{v0}.J{j0}.k{i}: Reaction probabilities in the $ith$ iteration. To see the last results take the largest value of i. First column is the energy, second column the total reaction probability for all the possible product channels. Third is the total probability for reactants and products channels (it should be unitiy). The fourth is the total reaction probability for the reaction $01+2 \rightarrow 02 + 1$. The subsequent columns present the reaction probability for each of the final vibrational product states.
- potr1r2.e1.1: Contains the PES used in the calculation. First column is $r$ and the second is $R$. Subsequent columns are the values of the PES for the various angles in the grid.
- Cvjprod.{i}: Cumulative reaction probabilites.
- prodwv.v00.Omg00

## Reaction probabilties for $J>0$
The input for the calculations for $J>0$ are similar to those for $J=0$. The inputs for $J=5$ and $J=10$ are provided in this repository.

The main difference can be found in the namelist `INPUTGRIDBASE` where the values of JTOT, IPARITY and IOMMAX have been set.
```
&INPUTGRIDBASE
 ...
 JTOT=     5, IPARITY=  -1
 IOMMIN=   0, IOMMAX=   5,
 ...
 /
```

```
&INPUTGRIDBASE
 ...
 JTOT=    10, IPARITY=   1
 IOMMIN=   0, IOMMAX=  10,
 ...
 /
```
In order to compute the cumulative reaction probabilities once the calculations for the different $J$ have finished we create a new folder at the level of the J's folders named Omg0 and we copy there the J folders with their reactive probabilities. We will also create a CRP folder at the same level than Omg0 with a single file `CalculatedJ.dat` with two columns: the J value and the effective rotational constant in eV that optimizes the J-shifting interpolation. This should be the directory structure at this point:
```
project/
├── func
│   └── ...
├── input.dat
├── J000
│   ├── S2prod.v00.J000.k00008
│   └── ...
├── J005
│   ├── S2prod.v00.J005.k00008
│   └── ...
├── J010
│   ├── S2prod.v00.J010.k00008
│   └── ...
├── Omg0
│   ├── J000
│   │   └── S2prod.v00.J000.k00008
│   ├── J005
│   │   └── S2prod.v00.J000.k00008
│   └── J010
│       └── S2prod.v00.J000.k00008
├── CRP
│   ├── CalculatedJ.dat
├── pes
├── pot
│   └── ...
```

To calculate the cumulative reaction probabilities we navigate to the CRP folder and execute the `crp.out` program located in the BIN directory in the madwave3 repository. The state-to-state reactive cross sections are stored in the CRP.vf{vf}.ef1 files. Similarly, the inelastic reaction probability can be computed through the `cip.out` program, which will output the CIP.vf{vf}.ef1 files.

In both cases the first column is the collision energy, second is $K_0^2$ in $\AA^2$ and the following are the cumulative probabilities over $J$, $p$, $\Omega$ and $\Omega_0$. Each of the columns corresponds to a final $j'$ value, from 0 to its maximum selected value.

