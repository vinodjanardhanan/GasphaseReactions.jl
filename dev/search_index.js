var documenterSearchIndex = {"docs":
[{"location":"run/#Running-the-code","page":"Execution","title":"Running the code","text":"","category":"section"},{"location":"run/","page":"Execution","title":"Execution","text":"On the Julia REPL ","category":"page"},{"location":"run/","page":"Execution","title":"Execution","text":"julia>using GasphaseChemistry\njulia>run(\"gaschem/gaschem.xml\",\"lib/\")","category":"page"},{"location":"run/#Input-file","page":"Execution","title":"Input file","text":"","category":"section"},{"location":"run/","page":"Execution","title":"Execution","text":"The method takes two arguments input_path and lib_dir. The input_files specifies the input XML file and lib_dir speficies the path to the data files folder where therm.dat is present. The structure of the XML input file is shown below.","category":"page"},{"location":"run/","page":"Execution","title":"Execution","text":"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<gaschem>\n        <molefractions>H2=0.25, H2O=0.25, O2=0.25, N2=0.25</molefractions>\n        <T>1173.0</T>\n        <p>1e5</p>\n        <gas_mech>h2o2.dat</gas_mech>\n</gaschem>","category":"page"},{"location":"run/","page":"Execution","title":"Execution","text":"The meaning of different XML tags is explained below.","category":"page"},{"location":"run/","page":"Execution","title":"Execution","text":"<gaschem> : tag specifying the model\n<molefractions> : mole fractions of the gas-phase species. Instead of mass fractions, mole fractions may also be \n<T>: Temperature in K\n<p> : Pressure in Pa\n<gas_mech> : tag specifying the name of the mechanism input file","category":"page"},{"location":"run/#Output","page":"Execution","title":"Output","text":"","category":"section"},{"location":"run/","page":"Execution","title":"Execution","text":"The code does not generate any file output.  An example of terminal output that transport_properties generates is shown below.","category":"page"},{"location":"run/","page":"Execution","title":"Execution","text":"T, p Conditions \n\nT(K): \t1173.0\np(Pa): \t100000.0\n          H2 \t -1.5707e-01 \t -3.1663e-04\n           O2 \t -1.5728e-01 \t -5.0327e-03\n          H2O \t -1.3727e-03 \t -2.4730e-05\n            H \t +1.8271e-02 \t +1.8416e-05\n            O \t +5.2291e-10 \t +8.3663e-12\n           OH \t +2.8130e-01 \t +4.7842e-03\n          HO2 \t +1.7313e-02 \t +5.7144e-04\n         H2O2 \t +0.0000e+00 \t +0.0000e+00\n           N2 \t +0.0000e+00 \t +0.0000e+00\n Sum of sources: -5.421010862427522e-19\n","category":"page"},{"location":"run/","page":"Execution","title":"Execution","text":"","category":"page"},{"location":"gchem/#Gasphase-chemistry","page":"Theory","title":"Gasphase chemistry","text":"","category":"section"},{"location":"gchem/#The-input-file","page":"Theory","title":"The input file","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The gas chemistry input file generally follows the chemkin input format. The following section briefly describes the setting up of the gasphase mechanism input file. For a complete example, the users may refer to the home pages of GRI Mech (http://combustion.berkeley.edu/gri-mech/) or the mechanisms hosted by Lawrence Livermore national lab (https://combustion.llnl.gov/mechanisms)","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The general structure of the input file is shown below.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"! GRI-Mech Version 3.0 7/30/99  CHEMKIN-II format\n! See README30 file at anonymous FTP site unix.sri.com, directory gri;\n! WorldWideWeb home page http://www.me.berkeley.edu/gri_mech/ or\n! through http://www.gri.org, under 'Basic  Research',\n! for additional information, contacts, and disclaimer\nELEMENTS\nO  H  C  N  AR\nEND\nSPECIES\nH2      H       O       O2      OH      H2O     HO2     H2O2\nC       CH      CH2     CH2(S)  CH3     CH4     CO      CO2\nHCO     CH2O    CH2OH   CH3O    CH3OH   C2H     C2H2    C2H3\nC2H4    C2H5    C2H6    HCCO    CH2CO   HCCOH   N       NH\nNH2     NH3     NNH     NO      NO2     N2O     HNO     CN\nHCN     H2CN    HCNN    HCNO    HOCN    HNCO    NCO     N2\nAR      C3H7    C3H8    CH2CHO  CH3CHO\nEND\n!THERMO\n! Insert GRI-Mech thermodynamics here or use in the default file\n!END\nREACTIONS\n2O+M<=>O2+M                              1.200E+17   -1.000        .00\nH2/ 2.40/ H2O/15.40/ CH4/ 2.00/ CO/ 1.75/ CO2/ 3.60/ C2H6/ 3.00/ AR/  .83/\n...\n...\nO+CO(+M)<=>CO2(+M)                       1.800E+10     .000    2385.00\n   LOW/ 6.020E+14     .000    3000.00/\nH2/2.00/ O2/6.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/3.50/ C2H6/3.00/ AR/ .50/\n...\n...\nH+CH2(+M)<=>CH3(+M)                      6.000E+14     .000        .00\n     LOW  /  1.040E+26   -2.760   1600.00/\n     TROE/   .5620  91.00  5836.00  8552.00/\nH2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/ .70/\n...\n...\n2HO2<=>O2+H2O2                           1.300E+11     .000   -1630.00\n DUPLICATE\n2HO2<=>O2+H2O2                           4.200E+14     .000   12000.00\n DUPLICATE\nEND ","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"Any line starting with ! is treated as a comment line. The actual mechanism specification begins with the ELEMENT keyword. Instead of ELEMENT, ELEM may also be used. The element specification ends with the keyword END, which is not mandatory. However, there shall be only one instance of ELEMENT or ELEM keyword. This part is entirely ignored by ReactionEngine. i.e. it does not read and stores the element information provided in the mechanism input file. The species specification follows the element specification. The species specification starts with the keyword SPECIES or SPEC and ends with END. The END keyword is not mandatory. The species names may be specified either in uppercase or lowercase letters. The reactions follow the species definition and start with the keyword REACTIONS or REAC and ends with the keyword END. The default values of the activation energy are CAL/MOL. If the activation energy units are different from CAL/MOL, that must follow the keyword REACTION or REAC. The activation energy units may be specified in CAL/MOLE, KCAL/MOLE, JOULES/MOLE, KELVINS, or EVOLTS. An example specification is shown below.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"REACTION JOULES/MOL","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The units of pre-exponential factors by default are cm-mol-s-K. Instead of cm-mol-s-K, cm-molecule-s-K may be used. An example specification is shown below","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"REACTION JOULES/MOL MOLECULES","category":"page"},{"location":"gchem/#Specifying-a-reaction","page":"Theory","title":"Specifying a reaction","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"An example of reaction stoichiometry is as follows.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"H+O2+H2O<=>HO2+H2O                       11.26E+18    -.760        .00","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"There are two parts to this reaction specification. The first part is the stoichiometric equation, and the second part is the reaction parameters. The stoichiometric equation represents the chemical transformation happening, which may be reversible or irreversible.  A reversible reaction is represented by <=> as shown in the above example. Instead of <=>, = may also be used for denoting a reversible reaction. An irreversible reaction, on the other hand, is represented by the symbol =>. An example is shown below.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"O2+C3H7=>OH+C2H5+CH2O                   2.410E+13     .000       .00","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The + symbol separates the reactant species. There may be space between the different participating species. The reaction's stoichiometry may be specified either in the uppercase or lowercase letter. The reaction parameters follow the stoichiometric equation. These are the Arrhenius parameters. The first is the pre-exponential factor, the second is the temperature exponent, and the third is the activation energy. All the three parameters must be specified for the mechanism parser to work.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The rate of progress of the reaction i,  r_i  for an Arrhenius reaction is calculated according to ","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"r_i = k_fi prod_k=1^N_g X_k^nu_ki^prime - k_ri prod_k=1^N_g X_k^nu_ki^primeprime","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"where","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"k_fi = A_fi T^beta exp(-ERT)","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The reverse reaction rate constant for a reversible reaction is calculated from the thermodynamic information. i.e","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"k_ri = frack_fiK_ci","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The equilibrium constant K_ci is related to K_pi according to ","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"K_ci = K_pi left( fracp_mathrmatmRT right)^sum_k=1^N_gnu_ki","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The equilibrium constant K_pi is calculated from the Gibb's free energy change of the reaction. N_g is the number of gasphase species and nu_ki = nu_ki^primeprime-nu_ki^prime","category":"page"},{"location":"gchem/#Auxiliary-information","page":"Theory","title":"Auxiliary information","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The chemkin gasphase mechanism format allows several auxiliary information to alter the rate of a reaction. The ones that are supported by ReactionEngine is described below. ","category":"page"},{"location":"gchem/#Third-body-reactions","page":"Theory","title":"Third body reactions","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"An example of a third body reaction is shown below.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"2O+M<=>O2+M                              1.200E+17   -1.000        .00\nH2/ 2.40/ H2O/15.40/ CH4/ 2.00/ CO/ 1.75/ CO2/ 3.60/ C2H6/ 3.00/ AR/  .83/","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"Here M is the third body. By default the third body collision efficiency is treated as 1. Any species that contribute differently (i.e efficiency is not 1) as a third body must be specified in a separate line that follows the reaction specification line. The species and its enhanced third body efficiency are separated by '/'. For instance, in the above example, the enhanced third body efficiency of H2 is 2.4. ","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"In the case of third body reactions, the Arrhenius rate is modified by the concentrations of the third body species. The rate of progress then becomes","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"r_i = sum_k=1^N_gleft( alpha_kiX_k right) left( k_fi prod_k=1^N_g X_k^nu_ki^prime - k_ri prod_k=1^N_g X_k^nu_ki^primeprime right)","category":"page"},{"location":"gchem/#Reaction-order","page":"Theory","title":"Reaction order","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The order of a reaction with respect to a species may be altered by specifying the FORD and RORD keywords.  An example is shown below.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"2H2+O2<=>2H2O        3.0e13    0.0e-16    0.0e-16\nFORD/H2 2.0/O2 1.0/\nRORD/H2O 1.0/","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"In the above example, the order of the forward reaction w.r.t H_2 is 2. By default the order of the reaction w.r.t the concentration of H_2 is 2 due to the value of the stoichiometric coefficient. The additional specification of order 2 w.r.t H_2 will make the reaction 4th order w.r.t. H_2. In ReactionEngine the rate for the above reaction translates into the following ","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"r = A_f expleft(frac-ERTright) H_2^2H_2^2O_2 - frack_fK_c H_2OH_2O","category":"page"},{"location":"gchem/#Pressure-dependent-reactions","page":"Theory","title":"Pressure dependent reactions","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"***ReactionEngine*** implements three different models for the rate calculation of pressure dependent reactions. The pressure dependent reactions auxiliary information line must follow the reaction specification. For fall off reactions, the LOW keyword follows the reaction information, and the parameters correspond A_0 beta_0, and E_0. For chemically activated bi-molecular reactions, the auxiliary information for the HIGH keyword must be provided. The parameters correspond to  A_infty beta_infty, and E_infty.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"With the high pressure and low pressure parameters, the rate constant is calculated as follows.","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"k_0 = A_0 T^beta_0 exp(-E_0RT)","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"k_infty = A_infty T^beta_infty exp(-E_inftyRT)","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"k = k_infty left( fracPr1+Pr right)F","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The reduced pressure Pr is given by","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"Pr = frack_0Mk_infty","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"where M is the concentration of the mixture, including third body efficiencies. For the Lindemaan form F=1","category":"page"},{"location":"gchem/#Troe-reactions","page":"Theory","title":"Troe reactions","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"In the case of Troe reactions, in addition to the high pressure and low pressure parameters, the parameters following the TROE keyword must be provided. The parameters are assumed to be in the following order a T^*** T^*, and T^**. The fourth parameter is optional. ","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"The value of F in Troe form is calculated as follows","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"log F = left 1+ left( fraclog Pr + cn-d(log Pr + c) right)^2 right^-1 log F_mathrmcent","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"c = -04-067log F_mathrmcent","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"n = 075 -127 log F_mathrmcent","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"d = 014","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"F_mathrmcent = (1-alpha) exp(-TT^***) + alpha exp( -TT^*) +exp(-T^**T)","category":"page"},{"location":"gchem/#Sri-reactions","page":"Theory","title":"Sri reactions","text":"","category":"section"},{"location":"gchem/","page":"Theory","title":"Theory","text":"In the case of SRI reactions, in addition to the low and high pressure limit parameters, the parameters following the SRI keyword must be provided. Either  3 or 5 parameters must be provided. The parameters are assumed in the following order a b c d and e. If only first three are specified, then d=1 and e=0. The F value in the case of Sri fall off reaction is calculated as follows","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"x = 1(1+log^2 Pr)","category":"page"},{"location":"gchem/","page":"Theory","title":"Theory","text":"F = dleft a exp(-bT)+exp(-Tc) right^x T^e","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GasphaseReactions","category":"page"},{"location":"#GasphaseReactions","page":"Home","title":"GasphaseReactions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GasphaseReactions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GasphaseReactions]","category":"page"},{"location":"#GasphaseReactions.calculate_molar_production_rates!-Tuple{Any, ReactionCommons.GasMechDefinition, Any}","page":"Home","title":"GasphaseReactions.calculate_molar_production_rates!","text":"calculatemolarproductionrates!(ms::MixtureState, gd::GasmechDefinition, thermoobj)    \n\nms: State object\ngd: Struct of the type GasMechDefinition \nthermo_obj: ThermoObj structure\n\nreturn value is in mol/m3\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.clean_up!-Union{Tuple{T}, Tuple{Array{T}, Array{T}}} where T<:AbstractString","page":"Home","title":"GasphaseReactions.clean_up!","text":"cleanup!(species::Array{T}, allspecies::Array{T}) where T <: AbstractString\n\nIn the case of fall-off reactions, the reactant and product species may contain \n\n\"(\" as part of the species name and this occurs at the end of the name. This needs to be cleaned up \n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.collect_species!-Union{Tuple{T}, Tuple{Array{T}, T}} where T<:AbstractString","page":"Home","title":"GasphaseReactions.collect_species!","text":"collectspecies!(gasphasespecies::Array{T}, data_string::T) \n\nThis function is not for extermal calls \n\nfunction for reading the species present in the gasphase mechanism from the SPECIES .. END block\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.compile_gaschemistry-Tuple{T} where T<:AbstractString","page":"Home","title":"GasphaseReactions.compile_gaschemistry","text":"compilegaschemistry(inputfile::T) \n\ninputfile : inputfile including the path \n\nFunction for reading the gasphase input file. The function returns the     definition of gasphase mechanism \n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.delete_tbs!-Tuple{Any, Vararg{Any}}","page":"Home","title":"GasphaseReactions.delete_tbs!","text":"delete_tbs(tbs, args...) function to delete the third body species from the reactant and product list \n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.fall_off_species-Tuple{AbstractString}","page":"Home","title":"GasphaseReactions.fall_off_species","text":"falloffspecies(rxn_string::AbstractString)\n\nFunction to find the species acting as third body in a fall-off\n\nreaction. The species could be (+M) or (+Sp), where sp is any species present in the mechanism\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.modify_coefficients!-Tuple","page":"Home","title":"GasphaseReactions.modify_coefficients!","text":"modify_coefficients!(args...)\n\nChemkin gasphase chemistry allows non unity stoichiometric coefficients.\n\nThis function identifies the coefficients of any species that is not 1 and     modifies the reactant and product species list accordingly\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.parse_reaction-Union{Tuple{T}, Tuple{T, Array{T}, Any}} where T<:AbstractString","page":"Home","title":"GasphaseReactions.parse_reaction","text":"function parsereaction(datastring::AbstractString)\n\nfunction for reading the reactions present in the REACTIONS... END block\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.parse_third_body_collision_data-Union{Tuple{T}, Tuple{Array{T}, T}} where T<:AbstractString","page":"Home","title":"GasphaseReactions.parse_third_body_collision_data","text":"parsethirdbodycollisiondata(gasphasespecies::Array{T}, datastring::T)\n\nFunction for parsing the third body collision data \n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.pressure_dependent_params-Tuple{AbstractString}","page":"Home","title":"GasphaseReactions.pressure_dependent_params","text":"parsepressuredependentmodelparams(gasphase_reactions::AbstractString) Function for parsing the pressure dependent model parameters TROE and SRI reactions\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.rate_constant-Tuple{Real, ReactionCommons.Arrhenius}","page":"Home","title":"GasphaseReactions.rate_constant","text":"Calculate Arrhenius rate constant\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.read_gaschem!-Union{Tuple{T}, Tuple{Array{T}, T}} where T<:AbstractString","page":"Home","title":"GasphaseReactions.read_gaschem!","text":"readgaschem!(gasphasespecies::Array{T}, inputfile::T, libdir::T) where T <: AbstractString\n\ngasphase_species : list of species present in the mechanism \ninput_file : input file including the path \n\nFunction for reading the gasphase mechanism file.     Only chemkin input file is supported     \n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.reduced_pressure-NTuple{7, Any}","page":"Home","title":"GasphaseReactions.reduced_pressure","text":"reduced_pressure(id, k, T, conc, low, high)    \n\ncalculate the  reduced pressure\nid: reaction id \nk : forward reaction rate constant \nT : mixture temperature\ntb_conc : third body concentration\nlow : struct Arrhenius \nhigh : struct Arrhenius \nunimol : whether unimolecular or bimolecular boolean\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.third_body_collision-Tuple{Int64, Array{Float64}, Dict{Int64, Dict{Int64, Float64}}}","page":"Home","title":"GasphaseReactions.third_body_collision","text":"thirdbodycollision(rxnid::Int64,conc::Array{Float64},tbcall::Dict{Int64, Dict{Int64, Float64}})\n\nFunction for the calculation of third body collision effciencies.   If third body effciency is not specified for a species then that        is assumed as 1. If no collision effciencies are specified then        [M] is nothing but the concentration of the mixture \n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.unit_conversion_A-Tuple{Array{ReactionCommons.Arrhenius}, Any}","page":"Home","title":"GasphaseReactions.unit_conversion_A","text":"unitconversionA(params::Array{Arrhenius}, factor)\n\nFunction for converting pre-exponential factor to units of moles    the units are in cm-mol-s-K\n\n\n\n\n\n","category":"method"},{"location":"#GasphaseReactions.unit_conversion_E-Tuple{Array{ReactionCommons.Arrhenius}, Any}","page":"Home","title":"GasphaseReactions.unit_conversion_E","text":"unitconversionE(params::Array{Arrhenius}, factor)\n\nFunction for converting the activation energy to SI units J/Mol \n\n\n\n\n\n","category":"method"}]
}
