# Gasphase chemistry
## The input file
The gas chemistry input file generally follows the chemkin input format. The following section briefly describes the setting up of the gasphase mechanism input file. For a complete example, the users may refer to the home pages of GRI Mech (http://combustion.berkeley.edu/gri-mech/) or the mechanisms hosted by Lawrence Livermore national lab (https://combustion.llnl.gov/mechanisms)

The general structure of the input file is shown below.
```
! GRI-Mech Version 3.0 7/30/99  CHEMKIN-II format
! See README30 file at anonymous FTP site unix.sri.com, directory gri;
! WorldWideWeb home page http://www.me.berkeley.edu/gri_mech/ or
! through http://www.gri.org, under 'Basic  Research',
! for additional information, contacts, and disclaimer
ELEMENTS
O  H  C  N  AR
END
SPECIES
H2      H       O       O2      OH      H2O     HO2     H2O2
C       CH      CH2     CH2(S)  CH3     CH4     CO      CO2
HCO     CH2O    CH2OH   CH3O    CH3OH   C2H     C2H2    C2H3
C2H4    C2H5    C2H6    HCCO    CH2CO   HCCOH   N       NH
NH2     NH3     NNH     NO      NO2     N2O     HNO     CN
HCN     H2CN    HCNN    HCNO    HOCN    HNCO    NCO     N2
AR      C3H7    C3H8    CH2CHO  CH3CHO
END
!THERMO
! Insert GRI-Mech thermodynamics here or use in the default file
!END
REACTIONS
2O+M<=>O2+M                              1.200E+17   -1.000        .00
H2/ 2.40/ H2O/15.40/ CH4/ 2.00/ CO/ 1.75/ CO2/ 3.60/ C2H6/ 3.00/ AR/  .83/
...
...
O+CO(+M)<=>CO2(+M)                       1.800E+10     .000    2385.00
   LOW/ 6.020E+14     .000    3000.00/
H2/2.00/ O2/6.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/3.50/ C2H6/3.00/ AR/ .50/
...
...
H+CH2(+M)<=>CH3(+M)                      6.000E+14     .000        .00
     LOW  /  1.040E+26   -2.760   1600.00/
     TROE/   .5620  91.00  5836.00  8552.00/
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/ .70/
...
...
2HO2<=>O2+H2O2                           1.300E+11     .000   -1630.00
 DUPLICATE
2HO2<=>O2+H2O2                           4.200E+14     .000   12000.00
 DUPLICATE
END 
```
Any line starting with **!** is treated as a comment line. The actual mechanism specification begins with the **ELEMENT** keyword. Instead of **ELEMENT**, **ELEM** may also be used. The element specification ends with the keyword **END**, which is not mandatory. However, there shall be only one instance of **ELEMENT** or **ELEM** keyword. This part is entirely ignored by **ReactionEngine**. i.e. it does not read and stores the element information provided in the mechanism input file. The species specification follows the element specification. The species specification starts with the keyword **SPECIES** or **SPEC** and ends with **END**. The **END** keyword is not mandatory. The species names may be specified either in uppercase or lowercase letters. The reactions follow the species definition and start with the keyword **REACTIONS** or **REAC** and ends with the keyword **END**. The default values of the activation energy are **CAL/MOL**. If the activation energy units are different from **CAL/MOL**, that must follow the keyword **REACTION** or **REAC**. The activation energy units may be specified in **CAL/MOLE, KCAL/MOLE, JOULES/MOLE, KELVINS**, or **EVOLTS**. An example specification is shown below.

```
REACTION JOULES/MOL
```

The units of pre-exponential factors by default are **cm-mol-s-K**. Instead of **cm-mol-s-K**, **cm-molecule-s-K** may be used. An example specification is shown below
```
REACTION JOULES/MOL MOLECULES
```



## Specifying a reaction
An example of reaction stoichiometry is as follows.
```
H+O2+H2O<=>HO2+H2O                       11.26E+18    -.760        .00
```
There are two parts to this reaction specification. The first part is the stoichiometric equation, and the second part is the reaction parameters. The stoichiometric equation represents the chemical transformation happening, which may be reversible or irreversible.  A reversible reaction is represented by *<=>* as shown in the above example. Instead of *<=>*, *=* may also be used for denoting a reversible reaction. An irreversible reaction, on the other hand, is represented by the symbol **=>**. An example is shown below.
```
O2+C3H7=>OH+C2H5+CH2O                   2.410E+13     .000       .00
```
The *+* symbol separates the reactant species. There may be space between the different participating species. The reaction's stoichiometry may be specified either in the uppercase or lowercase letter. The reaction parameters follow the stoichiometric equation. These are the Arrhenius parameters. The first is the pre-exponential factor, the second is the temperature exponent, and the third is the activation energy. All the three parameters must be specified for the mechanism parser to work.

The rate of progress of the reaction $i$,  $r_i$  for an Arrhenius reaction is calculated according to 
```math
r_i = k_{fi} \prod_{k=1}^{N_g} [X_k]^{\nu_{ki}^{\prime}} - k_{ri} \prod_{k=1}^{N_g} [X_k]^{\nu_{ki}^{\prime\prime}}
```
where
```math
k_{fi} = A_{fi} T^\beta \exp(-E/RT)
```
The reverse reaction rate constant for a reversible reaction is calculated from the thermodynamic information. i.e
```math
k_{ri} = \frac{k_{fi}}{K_{ci}}
```
The equilibrium constant $K_{ci}$ is related to $K_{pi}$ according to 
```math
K_{ci} = K_{pi} \left( \frac{p_\mathrm{atm}}{RT} \right)^{\sum_{k=1}^{N_g}\nu_{ki}}
```
The equilibrium constant $K_{pi}$ is calculated from the Gibb's free energy change of the reaction. $N_g$ is the number of gasphase species and $\nu_{ki} = \nu_{ki}^{\prime\prime}-\nu_{ki}^{\prime}$

## Auxiliary information
The chemkin gasphase mechanism format allows several auxiliary information to alter the rate of a reaction. The ones that are supported by **ReactionEngine** is described below. 

### Third body reactions
An example of a third body reaction is shown below.
```
2O+M<=>O2+M                              1.200E+17   -1.000        .00
H2/ 2.40/ H2O/15.40/ CH4/ 2.00/ CO/ 1.75/ CO2/ 3.60/ C2H6/ 3.00/ AR/  .83/
```
Here **M** is the third body. By default the third body collision efficiency is treated as 1. Any species that contribute differently (i.e efficiency is not 1) as a third body must be specified in a separate line that follows the reaction specification line. The species and its enhanced third body efficiency are separated by '/'. For instance, in the above example, the enhanced third body efficiency of H2 is 2.4. 


In the case of third body reactions, the Arrhenius rate is modified by the concentrations of the third body species. The rate of progress then becomes
```math
r_i = \sum_{k=1}^{N_g}\left( \alpha_{ki}[X_k] \right) \left( k_{fi} \prod_{k=1}^{N_g} [X_k]^{\nu_{ki}^{\prime}} - k_{ri} \prod_{k=1}^{N_g} [X_k]^{\nu_{ki}^{\prime\prime}} \right)
```

### Reaction order
The order of a reaction with respect to a species may be altered by specifying the **FORD** and **RORD** keywords. 
An example is shown below.
```
2H2+O2<=>2H2O        3.0e13    0.0e-16    0.0e-16
FORD/H2 2.0/O2 1.0/
RORD/H2O 1.0/
```
In the above example, the order of the forward reaction w.r.t $H_2$ is 2. By default the order of the reaction w.r.t the concentration of $H_2$ is 2 due to the value of the stoichiometric coefficient. The additional specification of order 2 w.r.t $H_2$ will make the reaction 4th order w.r.t. H$_2$. In **ReactionEngine** the rate for the above reaction translates into the following 
```math
r = A_f \exp\left(\frac{-E}{RT}\right) [H_2]^2[H_2]^2[O_2] - \frac{k_f}{K_c} [H_2O][H_2O]
```

### Pressure dependent reactions

***ReactionEngine*** implements three different models for the rate calculation of pressure dependent reactions. The pressure
dependent reactions auxiliary information line must follow the reaction specification. For fall off reactions, the **LOW** keyword
follows the reaction information, and the parameters correspond $A_0, \beta_0$, and $E_0$. For chemically activated bi-molecular
reactions, the auxiliary information for the **HIGH** keyword must be provided. The parameters correspond to  $A_{\infty}, \beta_{\infty}$, and $E_{\infty}$.

With the high pressure and low pressure parameters, the rate constant is calculated as follows.

```math
k_0 = A_0 T^{\beta_0} \exp(-E_0/RT)
```
```math
k_{\infty} = A_{\infty} T^{\beta_{\infty}} \exp(-E_{\infty}/RT)
```

```math
k = k_{\infty} \left( \frac{Pr}{1+Pr} \right)F
```
The reduced pressure $Pr$ is given by
```math
Pr = \frac{k_0[M]}{k_{}\infty}
```
where $[M]$ is the concentration of the mixture, including third body efficiencies. For the Lindemaan form $F=1$

#### Troe reactions

In the case of Troe reactions, in addition to the high pressure and low pressure parameters, the parameters following the **TROE** keyword must be
provided. The parameters are assumed to be in the following order $a, T^{***}, T^*$, and $T^{**}$. The fourth parameter is optional. 

The value of $F$ in Troe form is calculated as follows
```math
log F = \left[ 1+ \left( \frac{log Pr + c}{n-d(log Pr + c)} \right)^2 \right]^{-1} log F_\mathrm{cent}
```
```math
c = -0.4-0.67log F_\mathrm{cent}
```
```math
n = 0.75 -1.27 log F_\mathrm{cent}
```
```math
d = 0.14
```
```math
F_\mathrm{cent} = (1-\alpha) \exp(-T/T^{***}) + \alpha \exp( -T/T^*) +\exp(-T^{**}/T)
```

#### Sri reactions
In the case of SRI reactions, in addition to the low and high pressure limit parameters, the parameters following the **SRI** keyword must be provided. Either 
3 or 5 parameters must be provided. The parameters are assumed in the following order $a, b, c, d$ and $e$. If only first three are specified, then $d=1$ and $e=0$. The F value in the case of Sri fall off reaction is calculated as follows

```math
x = 1/(1+log^2 Pr)
```

```math
F = d\left[ a \exp(-b/T)+\exp(-T/c) \right]^{x} T^e
```
