## Running the code
On the Julia REPL 
```julia
julia>using GasphaseChemistry
julia>run("gaschem/gaschem.xml","lib/")
```

## Input file
The method takes two arguments *input\_path* and *lib\_dir*. The input\_files specifies the input XML file and
lib\_dir speficies the path to the data files folder where *therm.dat* is present. The structure of the XML input file is shown below.

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<gaschem>
        <molefractions>H2=0.25, H2O=0.25, O2=0.25, N2=0.25</molefractions>
        <T>1173.0</T>
        <p>1e5</p>
        <gas_mech>h2o2.dat</gas_mech>
</gaschem>
```

The meaning of different XML tags is explained below.

- <gaschem> : tag specifying the model
- <molefractions> : mole fractions of the gas-phase species. Instead of mass fractions, mole fractions may also be 
- <T>: Temperature in K
- <p> : Pressure in Pa
- <gas_mech> : tag specifying the name of the mechanism input file

## Output
The code does not generate any file output.  An example of terminal output that transport_properties generates is shown below.
```
T, p Conditions 

T(K): 	1173.0
p(Pa): 	100000.0
          H2 	 -1.5707e-01 	 -3.1663e-04
           O2 	 -1.5728e-01 	 -5.0327e-03
          H2O 	 -1.3727e-03 	 -2.4730e-05
            H 	 +1.8271e-02 	 +1.8416e-05
            O 	 +5.2291e-10 	 +8.3663e-12
           OH 	 +2.8130e-01 	 +4.7842e-03
          HO2 	 +1.7313e-02 	 +5.7144e-04
         H2O2 	 +0.0000e+00 	 +0.0000e+00
           N2 	 +0.0000e+00 	 +0.0000e+00
 Sum of sources: -5.421010862427522e-19

```

 