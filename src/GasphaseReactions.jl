module GasphaseReactions

using LightXML
using Printf
using IdealGas
using ReactionCommons
using RxnHelperUtils

include("Constants.jl")
   
export compile_gaschemistry, calculate_molar_production_rates!

"""
compile_gaschemistry(input_file::T) 
-   input_file : input_file including the path 

Function for reading the gasphase input file. The function returns the
    definition of gasphase mechanism 
"""
function compile_gaschemistry(input_file::T) where T <: AbstractString
    gasphase_species = Array{T,1}()
    gasphase_mech = read_gaschem!(gasphase_species, input_file)
    reaction_map = species_rxn_map(gasphase_species,gasphase_mech)       
    gmd = GasMechDefinition(gasphase_mech,reaction_map)
    return gmd
end




"""
read_gaschem!(gasphase_species::Array{T}, input_file::T, lib_dir::T) where T <: AbstractString
-   gasphase_species : list of species present in the mechanism 
-   input_file : input file including the path 
Function for reading the gasphase mechanism file.
    Only chemkin input file is supported     
"""
function read_gaschem!(gasphase_species::Array{T}, input_file::T) where T <: AbstractString
    # g_mech_file = lib_dir*input_file        
    # g_mech_file = get_path(lib_dir, input_file)
    
    reactions = Array{GasphaseReaction,1}()
    low_rxns = Dict{Int64,Arrhenius}()
    troe_rxns = Dict{Int64,Troe}()
    high_rxns = Dict{Int64,Arrhenius}()
    sri_rxns = Dict{Int64,Sri}()
    lt_rxns = Dict{Int64, LandauTeller}()
    rev_rxns = Dict{Int64,Arrhenius}()
    third_body_reactions = Dict{Int64,Dict{Int64,Float64}}() #Array{Order,1}()
    ford_params = Dict{Int64,Dict{Int64,Float64}}()
    rord_params = Dict{Int64,Dict{Int64,Float64}}()

    Ea_factor = 1.0
    A_factor = 1.0


    open(input_file) do io
        read_elements = false
        read_species = false
        read_rxns = false        
        rxn_id = 0
        
        while !eof(io)
            data_string = readline(io)
            if length(data_string) != 0 && SubString(data_string[1:1])!="!"
                keyword_check = strip(uppercase(split(data_string)[1]))
                if keyword_check == "ELEMENTS" || keyword_check == "ELEM"
                    read_elements = false                
                elseif keyword_check == "SPECIES" || keyword_check == "SPEC"
                    read_elements = true
                    read_species = false
                elseif keyword_check == "REACTIONS" || keyword_check == "REAC"
                    read_species = true
                    read_rxns = false    
                    unit_check = split(uppercase(data_string))
                    # The default units for Arrhenius parameters are cm-mol-s-K and cal/mol
                    energy_units = length(unit_check) > 1 ? unit_check[2] : "CAL/MOL"       
                    pre_exponential_units = length(unit_check) > 2 ? unit_check[3] : "MOLES"
                    if length(unit_check) > 0
                        Ea_factor = convert2si(string(energy_units))
                    end
                    if pre_exponential_units == "MOLECULES"
                        A_factor = 1.0/6.023e23
                    end
                end                
                if read_species == false && read_elements == true
                    # get the species participating in the reactions
                    collect_species!(gasphase_species, data_string)                          
                end
                if read_rxns == false && read_species == true && read_elements == true                     
                    # check if the datastring is a reaction string        
                             
                    if reaction_string(data_string)     
                    
                        rxn_id = rxn_id + 1                        
                        push!(reactions, parse_reaction(data_string, gasphase_species, rxn_id)) #Returns a Struct of Gasphase reactions 
                    else
                        if upstrip(data_string) == "REACTIONS" || upstrip(data_string) == "REAC" 
                            continue
                        elseif low(data_string)
                            low_rxns[rxn_id] = Arrhenius(arrhenius_params(strip(data_string))...)
                        elseif high(data_string)
                            high_rxns[rxn_id] = Arrhenius(arrhenius_params(strip(data_string))...)
                        elseif troe(data_string)
                            troe_rxns[rxn_id] = Troe(pressure_dependent_params(strip(data_string))...)                            
                        elseif sri(data_string)
                            sri_rxns[rxn_id] = Sri(pressure_dependent_params(strip(data_string))...)                                                       
                        elseif rev(data_string)                            
                            rev_rxns[rxn_id] = Arrhenius(parse_rev_rxn_params(strip(data_string))...)
                        elseif dup(data_string)
                            continue
                        elseif landau_teller(data_string)
                            lt_rxns[rxn_id] = LandauTeller(parse_lt_larams(strip(data_string))...)
                        elseif ford(data_string)
                            ford_params[rxn_id] = parse_order(gasphase_species,data_string)
                        elseif rord(data_string)
                            rord_params[rxn_id] = parse_order(gasphase_species,data_string)        
                        elseif upstrip(data_string) != "END"
                            third_body_reactions[rxn_id] = parse_third_body_collision_data(gasphase_species,data_string)                            
                        end
                    end
                end                
            end
        end

    end    
    
    # Conver the activation energy to Joules units if not specified in Jouls
    if Ea_factor != 1.0        
        rxn_params = [rxn.params for rxn in reactions]
        unit_conversion_E(rxn_params, Ea_factor)
        if !isempty(low_rxns) unit_conversion_E(collect(values(low_rxns)), Ea_factor) end 
        if !isempty(high_rxns) unit_conversion_E(collect(values(high_rxns)), Ea_factor) end 
        if !isempty(rev_rxns) unit_conversion_E(collect(values(rev_rxns)), Ea_factor) end
    end
    # Convert the pre-exponent to moles if specified in molecules 
    if A_factor != 1.0
        rxn_params = [rxn.params for rxn in reactions]
        unit_conversion_A(rxn_params, A_factor)
        if !isempty(low_rxns) unit_conversion_A(collect(values(low_rxns)), A_factor) end 
        if !isempty(high_rxns) unit_conversion_A(collect(values(high_rxns)), A_factor) end 
        if !isempty(rev_rxns) unit_conversion_A(collect(values(rev_rxns)), A_factor) end 
    end

    aux_data = AuxiliaryData(third_body_reactions,low_rxns,troe_rxns,high_rxns,sri_rxns,rev_rxns, lt_rxns,ford_params,rord_params)    
    GasphaseMechanism(gasphase_species,reactions,aux_data)    
end



"""
unit_conversion_E(params::Array{Arrhenius}, factor)
- Function for converting the activation energy to SI units J/Mol 
"""
function unit_conversion_E(params::Array{Arrhenius}, factor)
    for p in params
        p.E *= factor
    end
end


"""
unit_conversion_A(params::Array{Arrhenius}, factor)
- Function for converting pre-exponential factor to units of moles 
    the units are in cm-mol-s-K
"""
function unit_conversion_A(params::Array{Arrhenius}, factor)
    for p in params
        p.k0  *= factor
    end
end



"""
collect_species!(gasphase_species::Array{T}, data_string::T) 
-  This function is not for extermal calls 
function for reading the species present in the gasphase mechanism from the SPECIES .. END block
"""
function collect_species!(gasphase_species::Array{T}, data_string::T) where T <: AbstractString
    for sp in split(uppercase(data_string))
        if sp != "SPECIES" && sp != "SPEC" && sp != "END"
            push!(gasphase_species,sp)
        end
    end
end


"""
function parse_reaction(data_string::AbstractString)
- function for reading the reactions present in the REACTIONS... END block
"""
function parse_reaction(data_string::T, gasphase_species::Array{T}, rxn_id) where T <: AbstractString
    data_string = strip(uppercase(data_string))        
    if data_string != "REACTIONS" && data_string != "REAC" 
        bang = findfirst("!", data_string)    
        if bang != nothing            # this means there is a comment at the end of the reaction specification
            data_string =  data_string[1:bang[1]-1]
        end
        
        reversible = false         

        if reaction_string(data_string)  
            
            r_species_ids = Array{Int64,1}()
            p_species_ids = Array{Int64,1}()
                     
            reversible, r_species, p_species = parse_rxn_string(join(split(data_string)[1:end-3],""))                 
            
            k,β,E = parse_rxn_params(join(split(data_string)[end-2:end],"  "))
            
            
            # search for fall off reactions
            falloff = fall_off(data_string)
            third_body = false 
            tbs = ""                
            if falloff 
                # find the fall off species (+M) or (+N2) etc...            
                tbs = fall_off_species(data_string)                
                third_body = true
                # In the case of fall off reactions after parsing the reactant and product 
                # species will contain "(" as part of name. Clean that up! 
                clean_up!(r_species, gasphase_species)
                clean_up!(p_species, gasphase_species)
                delete_tbs!(tbs, r_species, p_species)
            end                                    
            # check for non unity stoichiometric coefficient
            modify_coefficients!(r_species, p_species)  
            # The following is a check for third body reactions 
            if length(tbs) == 0 && in("M", r_species)
                third_body = true             
            end
            # delete the third body species from the reactant and product list            
            delete_tbs!("M", r_species, p_species)          
            r_species_ids = map(x->get_index(x,gasphase_species),r_species)
            p_species_ids = map(x->get_index(x,gasphase_species),p_species)  
            
            return GasphaseReaction(rxn_id, falloff, third_body, tbs, Arrhenius(k,β,E), GasphaseRxnStoichiometry(reversible,r_species_ids,p_species_ids,data_string))
        end        
    end
end



upstrip(data_string::AbstractString) = uppercase(strip(data_string))

reaction_string(data_string::AbstractString)  = occursin("<=>",data_string) || occursin("=",data_string) || occursin("=>",data_string)
low(data_string::AbstractString)  = occursin("LOW", uppercase(data_string))   
high(data_string::AbstractString)  = occursin("HIGH", uppercase(data_string))   
troe(data_string::AbstractString)  = occursin("TROE", uppercase(data_string))   
sri(data_string::AbstractString)  = occursin("SRI", uppercase(data_string))   
rev(data_string::AbstractString)  = occursin("REV", uppercase(data_string))   
dup(data_string::AbstractString)  = occursin("DUP", uppercase(data_string))   
ford(data_string::AbstractString)  = occursin("FORD", uppercase(data_string))   
rord(data_string::AbstractString)  = occursin("RORD", uppercase(data_string))   
fall_off(data_string::AbstractString) = occursin("(+", uppercase(data_string))
landau_teller(data_string::AbstractString) = occursin("LT", uppercase(data_string))
    

"""
clean_up!(species::Array{T}, all_species::Array{T}) where T <: AbstractString
- In the case of fall-off reactions, the reactant and product species may contain 
"(" as part of the species name and this occurs at the end of the name. This
needs to be cleaned up 
"""
function clean_up!(species::Array{T}, all_species::Array{T}) where T <: AbstractString
    
    for i in eachindex(species)
        if !in(species[i], all_species)
            if species[i][end:end] == "(" || species[i][end:end] == ")"
                species[i] = species[i][1:end-1]                
            end
        end
    end    
end


mutable struct NonAllocatingArray{T <: Real}
    g_all::Array{T}    
    Kp::Array{T}
end
export NonAllocatingArray

"""
calculate_molar_production_rates!(ms::MixtureState, gd::GasmechDefinition, thermo_obj)    
-   ms: State object
-   gd: Struct of the type GasMechDefinition 
-   thermo_obj: ThermoObj structure
return value is in mol/m3
"""
function calculate_molar_production_rates!(ms, gd::GasMechDefinition, thermo_obj)
    # calculate the bulk gasphase concentration in mol/cm3 from the mixture state 
    Cb = ms.p*1e-6/R/ms.T  # 1e-6 to convert to mol/cm3
    # concentration of individual species 
    for k in 1:length(ms.mole_frac)
        ms.conc[k] = Cb * ms.mole_frac[k]
    end
    # Get the Gibbs free energy values for the calculation of Equilibrium constant     
    ms.g_all = IdealGas.H_all(thermo_obj, ms.T) - ms.T * IdealGas.S_all(thermo_obj, ms.T)        
    # Calculate Del G for individual reactions
    for rxn in gd.gm.reactions
        delG = sum(ms.g_all[rxn.rxn_stoic.product_ids]) - sum(ms.g_all[rxn.rxn_stoic.reactant_ids])               
        # Calculate equilibrium constant for individual reactions 
        ms.Kp[rxn.id] = exp(-delG/R/ms.T)
       
    end
    # calculation of forward reaction rate constant 
    for rxn in gd.gm.reactions                
        k_forward = rate_constant(ms.T,rxn.params)   
        third_body_conc = 1.0                     
        if rxn.third_body            
            # In the case of fall off reactions an individual species can act as third body 
            if in(rxn.third_body_species, gd.gm.species)                
                third_body_conc = ms.conc[get_index(rxn.third_body_species,gd.gm.species)]
            else
                third_body_conc = third_body_collision(rxn.id, ms.conc, gd.gm.aux_data.third_body_reaction)                 
            end
        end
        
        # check for fall off reactions
        if rxn.fall_off                                 
            if length(rxn.reactant_ids) == 1 || length(rxn.product_ids) == 1
                chk = 1
            end
            Pr, k_chk = reduced_pressure(rxn.id, k_forward, ms.T, third_body_conc, gd.gm.aux_data.low_reactions, gd.gm.aux_data.high_reactions, chk)
            # By default Lindemann is assumed i.e. F=1
            F = 1.0            
            # Troe and Sri models are supported for fall off reactions 
            if in(rxn.id, collect(keys(gd.gm.aux_data.troe_reactions)))            
                fcent = Fcent(ms.T, gd.gm.aux_data.troe_reactions[rxn.id])                
                log10_f = log10f(Pr, fcent)
                F = 10^log10_f
            elseif in(rxn.id, collect(keys(gd.gm.aux_data.sri_reactions)))
                F = sri_f(gd.gm.aux_data.sri_reactions[rxn.id], Pr, ms.T)                
            end
            
            rPr = 1/(1+Pr)
            if chk == 1
                rPr = rPr*Pr
            end
            
            k_forward = k_chk * rPr * F  
        end
                
        # Calcaulate the product of reactant concentrations                
        if in(rxn.id, collect(keys(gd.gm.aux_data.f_order)))      
            prdt_rct_conc = prod_conc(ms.conc, rxn.rxn_stoic.reactant_ids, gd.gm.aux_data.f_order[rxn.id])
        else
            prdt_rct_conc = prod_conc(ms.conc, rxn.rxn_stoic.reactant_ids)            
        end
        # Rate of the forward reaction 
        ms.g_rxn_rate[rxn.id] = k_forward*prdt_rct_conc
        # If reversible get the net rate                         
        k_reverse = 0.0
        if rxn.rxn_stoic.reversible
            # Calculate the product of product concentrations
            if in(rxn.id, collect(keys(gd.gm.aux_data.r_order)))          
                prod_prdt_conc = prod_conc(ms.conc, rxn.rxn_stoic.product_ids , gd.gm.aux_data.r_order[rxn.id])
            else
                prod_prdt_conc = prod_conc(ms.conc, rxn.rxn_stoic.product_ids)
            end
            # Check if explicit reverse reaction rate is specified   
            if in(rxn.id, collect(keys(gd.gm.aux_data.rev_reactions)))
                
                k_reverse = rate_constant(ms.T, gd.gm.aux_data.rev_reactions[rxn.id])
                
            else
                #Calculate the Equilibrim constant Kc from Kp 
                ν = length(rxn.rxn_stoic.product_ids) - length(rxn.rxn_stoic.reactant_ids)
                Kc = ms.Kp[rxn.id]*(p_std/(R*ms.T))^ν                
                k_reverse = k_forward/Kc
                
            end
            # rate of individual reactions             
            ms.g_rxn_rate[rxn.id] -= k_reverse*prod_prdt_conc
            
        end
        ms.g_rxn_rate[rxn.id] *= third_body_conc
    end
    #calculate the species source terms  in mol/m3
    for srm in gd.rxn_species_array        
        ms.source[srm.species_id] = sum(ms.g_rxn_rate[srm.rxns] .* srm.stoic_coeff)*1e6
    end
end



"""
third_body_collision(rxn_id::Int64,conc::Array{Float64},tbc_all::Dict{Int64, Dict{Int64, Float64}})
- Function for the calculation of third body collision effciencies.
    If third body effciency is not specified for a species then that 
        is assumed as 1. If no collision effciencies are specified then 
        [M] is nothing but the concentration of the mixture 
"""
function third_body_collision(rxn_id::Int64,conc::Array{Float64},tbc_all::Dict{Int64, Dict{Int64, Float64}})
    Cb = sum(conc)
    if in(rxn_id,collect(keys(tbc_all)))
        rxn_tbc = tbc_all[rxn_id]
        species_ids = collect(keys(rxn_tbc))
        collision_eff = collect(values(rxn_tbc))        
        return   Cb + sum(conc[species_ids] .* collision_eff) - sum(conc[species_ids])
    else
        return Cb
    end
end


"""
educed_pressure(id, k, T, conc, low, high)    
-   calculate the pressure dependent rate
-   id: reaction id 
-   k : forward reaction rate constant 
-   T : mixture temperature
-   tb_conc : third body concentration
-   low : struct Arrhenius 
-   high : struct Arrhenius 
"""
function reduced_pressure(id, k, T, conc, low, high, chk)    
    
    if in(id, collect(keys(low)))
        k0 = rate_constant( T, low[id])
        k∞ = k
    end
    if in(id, collect(keys(high)))
        k∞ = rate_constant(T, high[id])
        k0 = k
    end
    if chk == 1
        return k0*conc/k∞, k∞
    else
        return k0*conc/k∞, k0
    end
end

expval(T,TS) = TS == 0 ? 0 : exp(-T/TS)
T2term(T,TS) = TS == 0 ? 0 : exp(-TS/T)
Fcent(T, troe_p) = (1-troe_p.a)*expval(T,troe_p.T3star)+troe_p.a*expval(T,troe_p.Tstar) + T2term(T, troe_p.T2star)
sri_f(sri::Sri, Pr, T) = begin 
    x = 1/(1+log10(Pr)^2)
    f1 = (sri.a*exp(-sri.b/T) + exp(-T/sri.c))^(x)
    sri.d*f1*T^sri.e
end

prod_conc(conc::Array{Float64}, ids::Array{Int}) = prod(conc[ids])

function prod_conc(conc::Array{Float64}, ids::Array{Int}, sp_order::Dict{Int, Float64})
    # Get the species id from the order dict 
    order_change_species_ids = collect(keys(sp_order))
    # Species that are not present in the order dict 
    first_order_species = setdiff(ids, order_change_species_ids)
    p1 = prod(conc[first_order_species])
    p2 = prod(conc[order_change_species_ids].^ collect(values(sp_order)))
    p1*p2
end

function log10f(Pr, fcent) 
    log_fcent = log10(fcent)
    c = -4-0.67log_fcent
    n = 0.75 - 1.27log_fcent
    d = 0.14    
    nr = log10(Pr)+c
    f1 = nr/(n-d*nr)
    return log_fcent/(1+f1^2)
end
"""
Calculate Arrhenius rate constant
"""
rate_constant(T::Real, arr::Arrhenius) = arr.k0 * T^arr.β * exp(-arr.E/R/T)



"""
parse_third_body_collision_data(gasphase_species::Array{T}, data_string::T)
-   Function for parsing the third body collision data 
"""
function parse_third_body_collision_data(gasphase_species::Array{T}, data_string::T) where T <: AbstractString    
    tbc = Dict{Int64,Float64}()
    data_string = (strip(data_string))[end:end] == "/" ? uppercase(data_string[1:end-1]) : uppercase(data_string)    
    data =  split(data_string,"/")              
    i = findall(x->x=="", data) #The split operation sometimes creates and empty string as last entry 
    if length(i) != 0 
        deleteat!(data,i) #remove the empty entry 
    end         
    for i in 1:2:length(data)        
        tbc[get_index(string(strip(data[i])),gasphase_species)] = parse(Float64,string(strip(data[i+1])))
    end
    tbc
end

function parse_order(gasphase_species::Array{T}, data_string::T) where T <: AbstractString
    ord = Dict{Int,Real}()
    data_string = (strip(data_string))[end:end] == "/" ? uppercase(data_string[1:end-1]) : uppercase(data_string)
    data = split(data_string[findfirst('/',data_string)+1:end],"/")
    for it in data
        ord[get_index(string(split(it)[1]),gasphase_species)] = parse(Float64,string(split(it)[2]))        
    end
    ord
end


"""
fall_off_species(rxn_string::AbstractString)
- Function to find the species acting as third body in a fall-off
reaction. The species could be (+M) or (+Sp), where sp is any
species present in the mechanism
"""
function fall_off_species(rxn_string::AbstractString)    
    done = false
    a = findfirst('(', rxn_string)
    while !done
        if rxn_string[a+1:a+1] == "+"
            done = true
            break
        else            
            a = findnext('(',rxn_string, a+1)
            if a == nothing
                break 
            end
        end
    end
    
    if a != nothing 
        b = findnext(')', rxn_string, a)        
        if b != nothing    
            return rxn_string[a+2:b-1]        
        end
    else        
        nothing
    end
    
end





function arrhenius_params(gasphase_reactions::AbstractString)            
    gasphase_reactions = gasphase_reactions[end:end] == "/" ? gasphase_reactions[1:end-1] : gasphase_reactions        
    rxn = String(gasphase_reactions[findfirst('/',gasphase_reactions)+1:end])    
    return parse_rxn_params(rxn)
end



function parse_rev_rxn_params(gasphase_reactions::AbstractString)
    return arrhenius_params(gasphase_reactions)    
end


function parse_lt_larams(lt_params::AbstractString)
    params = lt_params[end:end] == "/" ? lt_params[1:end-1] : lt_params
    params = String(params[findfirst('/',params)+1:end])
    rxn_params = split(params)
    a = parse(Float64, rxn_params[1])
    b = parse(Float64, rxn_params[2])
    return a, b
end

"""
parse_pressure_dependent_model_params(gasphase_reactions::AbstractString)
Function for parsing the pressure dependent model parameters TROE and SRI reactions
"""
function pressure_dependent_params(gasphase_reactions::AbstractString)        
    gasphase_reactions = gasphase_reactions[end:end] == "/" ? gasphase_reactions[1:end-1] : gasphase_reactions    
    if occursin("TROE", gasphase_reactions)        
        troe_params = split(gasphase_reactions[findfirst('/',gasphase_reactions)+1:end])
        return Tuple(map(x->parse(Float64,x),troe_params))
    elseif occursin("SRI", gasphase_reactions)
        sri_str = split(gasphase_reactions[findfirst('/',gasphase_reactions)+1:end])
        sri_params = (map(x->parse(Float64,x),sri_str))
        if length(sri_params)==3
            push!(sri_params,1.)
            push!(sri_params,0.)
        end
        return Tuple(sri_params)
    end

end



"""
modify_coefficients!(args...)
- Chemkin gasphase chemistry allows non unity stoichiometric coefficients.
This function identifies the coefficients of any species that is not 1 and
    modifies the reactant and product species list accordingly
"""
function modify_coefficients!(args...)
    function modify!(species)
        @inbounds for i in eachindex(species)
            if isdigit(species[i][1])
                m = parse(Int, string(species[i][1]))
                species[i] = species[i][2:end]
                for k in 1:m-1
                    push!(species, species[i])
                end
            end
        end  
    end

    for species_list in args 
        modify!(species_list)
    end    
end


"""
delete_tbs(tbs, args...)
function to delete the third body species from the reactant and product list 
"""
function delete_tbs!(tbs, args...)
    for species_list in args 
        i = findall(x->x==tbs, species_list)
        if length(i) != 0
            deleteat!(species_list, i)
        end        
    end
end


# End of module
end
