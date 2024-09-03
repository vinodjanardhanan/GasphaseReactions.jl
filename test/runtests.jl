using GasphaseReactions
using Test
using RxnHelperUtils
using IdealGas
using ReactionCommons
using LightXML, Printf


@testset "GasphaseReactions.jl" begin
    
    if Sys.isapple() || Sys.islinux()
        lib_dir = "lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end

    @testset "Testing h2o2 chemistry " begin
        input_file = "gaschem.xml"
        input_file = joinpath("h2o2", "gaschem.xml")
    
        xmldoc = parse_file(input_file)
        xmlroot = root(xmldoc)

        local mech_file = get_text_from_xml(xmlroot, "gas_mech")

        mech_file = get_path(lib_dir, mech_file)
        gas_mech_def = compile_gaschemistry(mech_file)
        
        gasphase = gas_mech_def.gm.species

        thermo_file = lib_dir*"therm.dat"
        thermo_obj =  IdealGas.create_thermo(gasphase,thermo_file)

        mole_fracs = get_molefraction_from_xml(xmlroot,thermo_obj.molwt, gasphase)
        local T = get_value_from_xml(xmlroot,"T")
        local p = get_value_from_xml(xmlroot,"p")
        
        conc = similar(mole_fracs)
        source = zeros(length(conc))
        g_all = similar(mole_fracs)    
        Kp = zeros(length(gas_mech_def.gm.reactions))    
        rxn_rate = zeros(length(Kp))
        # ms = MixtureState(T,p,mole_fracs,conc, source, rxn_rate)
        ms = GasphaseState(T,p,mole_fracs,conc,rxn_rate,source, g_all, Kp)    
        calculate_molar_production_rates!(ms, gas_mech_def, thermo_obj)
        sum_source = 0.0
        println("T, p Conditions ")
        println("\nT(K): \t",ms.T)
        println("p(Pa): \t", ms.p)
        for i in 1:length(ms.source)        
            @printf("%12s \t %+.4e \t %+.4e\n ", gasphase[i], ms.source[i], ms.source[i]*thermo_obj.molwt[i])
            sum_source += ms.source[i]*thermo_obj.molwt[i]
        end
        println("Sum of sources: ", sum_source)
        @test abs(sum_source) < 1e-9
    
    end

    @testset "Testing GRI Mech " begin
        input_file = "gaschem.xml"
        input_file = joinpath("gri", "gaschem.xml")
    
        xmldoc = parse_file(input_file)
        xmlroot = root(xmldoc)

        local mech_file = get_text_from_xml(xmlroot, "gas_mech")

        mech_file = get_path(lib_dir, mech_file)
        gas_mech_def = compile_gaschemistry(mech_file)
        
        gasphase = gas_mech_def.gm.species

        thermo_file = lib_dir*"therm.dat"
        thermo_obj =  IdealGas.create_thermo(gasphase,thermo_file)

        mole_fracs = get_molefraction_from_xml(xmlroot,thermo_obj.molwt, gasphase)
        local T = get_value_from_xml(xmlroot,"T")
        local p = get_value_from_xml(xmlroot,"p")
        
        conc = similar(mole_fracs)
        source = zeros(length(conc))
        g_all = similar(mole_fracs)    
        Kp = zeros(length(gas_mech_def.gm.reactions))    
        rxn_rate = zeros(length(Kp))
        # ms = MixtureState(T,p,mole_fracs,conc, source, rxn_rate)
        ms = GasphaseState(T,p,mole_fracs,conc,rxn_rate,source, g_all, Kp)    
        calculate_molar_production_rates!(ms, gas_mech_def, thermo_obj)
        sum_source = 0.0
        println("T, p Conditions ")
        println("\nT(K): \t",ms.T)
        println("p(Pa): \t", ms.p)
        for i in 1:length(ms.source)        
            @printf("%12s \t %+.4e \t %+.4e\n ", gasphase[i], ms.source[i], ms.source[i]*thermo_obj.molwt[i])
            sum_source += ms.source[i]*thermo_obj.molwt[i]
        end
        println("Sum of sources: ", sum_source)
        @test abs(sum_source) < 1e-8
    
    end

            
end
