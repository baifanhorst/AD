###################################################################################
# Right-hand sides of all equations
def RHS(t, y):
    
    # rhs stores all the righ-hand sides to be returned
    rhs = []
    # index of a variable
    index = 0
    
    ## Variables
    
    # ISF
    ISF_AbM = y[index]; index += 1    
    ISF_AbON = y[index]; index += 1
    ISF_AbFN = y[index]; index += 1
    ISF_AbF = y[index]; index += 1
    ISF_AbP = y[index]; index += 1
    ISF_mAb = y[index]; index += 1
    ISF_AbF_mAb = y[index]; index += 1
    ISF_AbM_mAb = y[index]; index += 1
    ISF_AbFN1 = y[index]; index += 1 # one antibody bound to a fibril end
    ISF_AbF1 = y[index]; index += 1 
    ISF_AbFN2 = y[index]; index += 1 # two antibody bound to two fibril ends
    ISF_AbF2 = y[index]; index += 1 
    
    ISF_AbFN_total = ISF_AbFN + ISF_AbFN1 + ISF_AbFN2
    ISF_AbF_total = ISF_AbF + ISF_AbF1 + ISF_AbF2
    
    
    # PL
    PL_AbM = y[index]; index += 1  
    PL_mAb = y[index]; index += 1    
    PL_AbM_mAb = y[index]; index += 1    
    
    
    # CSF
    CSF_AbM = y[index]; index += 1    
    CSF_mAb = y[index]; index += 1  
    CSF_AbM_mAb = y[index]; index += 1  
    
    # PE
    PE_mAb = y[index]; index += 1    
    
    
    
    ## Derivatives
    
    # ISF
    d_ISF_AbM_dt = 0
    d_ISF_AbON_dt = 0 
    d_ISF_AbFN_dt = 0 
    d_ISF_AbF_dt = 0
    d_ISF_AbP_dt = 0 
    d_ISF_mAb_dt = 0
    d_ISF_AbF_mAb_dt = 0
    d_ISF_AbM_mAb_dt = 0
    d_ISF_AbFN1_dt = 0
    d_ISF_AbF1_dt = 0
    d_ISF_AbFN2_dt = 0
    d_ISF_AbF2_dt = 0 
    # PL
    d_PL_AbM_dt = 0   
    d_PL_mAb_dt = 0
    d_PL_AbM_mAb_dt = 0
    # CSF
    d_CSF_AbM_dt = 0   
    d_CSF_mAb_dt = 0
    d_CSF_AbM_mAb_dt = 0
    # PE
    d_PE_mAb_dt = 0
    
    
    
    
    ## Calculating the right-hand sides

    # ISF
    # AbM production
    d_ISF_AbM_dt += k_ISF_AbM_syn
    # Primary nucleation
    d_ISF_AbM_dt += -k_ISF_AbO_1stNuc * (ISF_AbM**n_ISF_AbO_1stNuc) * n_ISF_AbO
    d_ISF_AbON_dt += k_ISF_AbO_1stNuc * (ISF_AbM**n_ISF_AbO_1stNuc)
    # Oligomer dissociation
    d_ISF_AbON_dt += -k_ISF_AbO_diss * ISF_AbON
    d_ISF_AbM_dt += k_ISF_AbO_diss * ISF_AbON * n_ISF_AbO
    # Oligomer conversion into fibril
    d_ISF_AbON_dt += -k_ISF_AbF_conv * ISF_AbON * (ISF_AbM**n_ISF_AbF_conv)
    d_ISF_AbFN_dt += k_ISF_AbF_conv * ISF_AbON * (ISF_AbM**n_ISF_AbF_conv)
    d_ISF_AbF_dt += k_ISF_AbF_conv * ISF_AbON * (ISF_AbM**n_ISF_AbF_conv) * n_ISF_AbO
    # Fibril elongation
    d_ISF_AbF_dt += 2 * k_ISF_AbF_on * ISF_AbM * ISF_AbFN
    d_ISF_AbM_dt += -2 * k_ISF_AbF_on * ISF_AbM * ISF_AbFN
    d_ISF_AbF1_dt += k_ISF_AbF_on * ISF_AbM * ISF_AbFN1
    d_ISF_AbM_dt += -k_ISF_AbF_on * ISF_AbM * ISF_AbFN1
    # Secondary nucleation
    d_ISF_AbM_dt += -k_ISF_AbO_2ndNuc * (ISF_AbM**n_ISF_AbO_2ndNuc) * ISF_AbF_total * n_ISF_AbO
    d_ISF_AbON_dt += k_ISF_AbO_2ndNuc * (ISF_AbM**n_ISF_AbO_2ndNuc) * ISF_AbF_total
    # Fibril conversion into plaque
    d_ISF_AbFN_dt += -k_ISF_AbP_conv * ISF_AbFN
    d_ISF_AbF_dt += -k_ISF_AbP_conv * ISF_AbF
    d_ISF_AbP_dt += k_ISF_AbP_conv * ISF_AbF
    d_ISF_AbFN1_dt += -k_ISF_AbP_conv * ISF_AbFN1 # Fibril end binding might affect plaque formation
    d_ISF_AbF1_dt += -k_ISF_AbP_conv * ISF_AbF1
    d_ISF_AbP_dt += k_ISF_AbP_conv * ISF_AbF1
    d_ISF_AbFN2_dt += -k_ISF_AbP_conv * ISF_AbFN2
    d_ISF_AbF2_dt += -k_ISF_AbP_conv * ISF_AbF2
    d_ISF_AbP_dt += k_ISF_AbP_conv * ISF_AbF2
    # Plaque dissociation into monomer
    d_ISF_AbM_dt += k_ISF_AbP_diss * ISF_AbP
    d_ISF_AbP_dt += -k_ISF_AbP_diss * ISF_AbP
    
    
    
    # Clearance of all species, all compartments
    
    # ISF
    # Monomer clearance
    d_ISF_AbM_dt += -k_ISF_AbM_clear * ISF_AbM
    # Fibril clearance
    d_ISF_AbF_dt += -k_ISF_AbF_clear * ISF_AbF
    d_ISF_AbFN_dt += -k_ISF_AbF_clear * ISF_AbFN
    # Plaque clearance
    d_ISF_AbP_dt += -k_ISF_AbP_clear * ISF_AbP
    # Antibody clearance
    d_ISF_mAb_dt += -k_ISF_mAb_clear * ISF_mAb # No such clearance in both Lin and Madrasi
    # Monomer-antibody clearance
    d_ISF_AbM_mAb_dt += -k_ISF_AbM_mAb_clear * ISF_AbM_mAb
    # Fibril surface-antibody clearance
    d_ISF_AbF_mAb_dt += -k_ISF_AbF_mAb_clear * ISF_AbF_mAb
    # Fibril end-entibody clearance
    d_ISF_AbFN1_dt += -k_ISF_AbFN1_clear * ISF_AbFN1
    d_ISF_AbF1_dt += -k_ISF_AbF1_clear * ISF_AbF1
    d_ISF_AbFN2_dt += -k_ISF_AbFN2_clear * ISF_AbFN2
    d_ISF_AbF2_dt += -k_ISF_AbF2_clear * ISF_AbF2
    
    
    # PL
    # Monomer clearance
    d_PL_AbM_dt += -k_PL_AbM_clear * PL_AbM
    # Antibody clearance 
    d_PL_mAb_dt += -k_PL_mAb_clear * PL_mAb
    # Monomer-antibody clearance
    d_PL_AbM_mAb_dt += -k_PL_AbM_mAb_clear * PL_AbM_mAb
    
    
    
    # Abeta monomer transportation
    # ISF -> CSF
    d_ISF_AbM_dt += -k_ISF_CSF_AbM_tran * ISF_AbM
    d_CSF_AbM_dt += k_ISF_CSF_AbM_tran * ISF_AbM
    # ISF -> PL
    d_ISF_AbM_dt += -k_ISF_PL_AbM_tran * ISF_AbM
    d_PL_AbM_dt +=k_ISF_PL_AbM_tran * ISF_AbM
    # PL -> ISF
    d_ISF_AbM_dt += k_PL_ISF_AbM_tran * PL_AbM
    d_PL_AbM_dt += -k_PL_ISF_AbM_tran * PL_AbM
    # CSF -> PL
    d_CSF_AbM_dt += -k_CSF_PL_AbM_tran * CSF_AbM
    d_PL_AbM_dt += k_CSF_PL_AbM_tran * CSF_AbM
    # PL -> CSF
    d_CSF_AbM_dt += k_PL_CSF_AbM_tran * PL_AbM
    d_PL_AbM_dt += -k_PL_CSF_AbM_tran * PL_AbM
    
    
    
    # Antibody transportation
    # ISF -> PL
    d_ISF_mAb_dt += -k_ISF_PL_mAb_tran * ISF_mAb
    d_PL_mAb_dt += k_ISF_PL_mAb_tran * ISF_mAb
    # PL -> ISF
    d_ISF_mAb_dt += k_PL_ISF_mAb_tran * PL_mAb
    d_PL_mAb_dt += -k_PL_ISF_mAb_tran * PL_mAb
    # ISF -> CSF
    d_ISF_mAb_dt += -k_ISF_CSF_mAb_tran * ISF_mAb
    d_CSF_mAb_dt += k_ISF_CSF_mAb_tran * ISF_mAb
    # CSF -> PL
    d_CSF_mAb_dt += -k_CSF_PL_mAb_tran * CSF_mAb
    d_PL_mAb_dt += k_CSF_PL_mAb_tran * CSF_mAb
    # PL -> CSF
    d_CSF_mAb_dt += k_PL_CSF_mAb_tran * PL_mAb
    d_PL_mAb_dt += -k_PL_CSF_mAb_tran * PL_mAb
    # PL -> PE
    d_PL_mAb_dt += -k_PL_PE_mAb_tran * PL_mAb
    d_PE_mAb_dt += k_PL_PE_mAb_tran * PL_mAb
    # PE -> PL
    d_PL_mAb_dt += k_PE_PL_mAb_tran * PE_mAb
    d_PE_mAb_dt += -k_PE_PL_mAb_tran * PE_mAb
    
    
    
    # mAb-AbM transportation
    # ISF->PL
    d_ISF_AbM_mAb_dt += -k_ISF_PL_AbM_mAb_tran * ISF_AbM_mAb
    d_PL_AbM_mAb_dt += k_ISF_PL_AbM_mAb_tran * ISF_AbM_mAb
    # PL->ISF
    d_ISF_AbM_mAb_dt += k_PL_ISF_AbM_mAb_tran * PL_AbM_mAb
    d_PL_AbM_mAb_dt += -k_PL_ISF_AbM_mAb_tran * PL_AbM_mAb
    # ISF->CSF
    d_ISF_AbM_mAb_dt += -k_ISF_CSF_AbM_mAb_tran * ISF_AbM_mAb
    d_CSF_AbM_mAb_dt += k_ISF_CSF_AbM_mAb_tran * ISF_AbM_mAb
    # CSF->PL
    d_CSF_AbM_mAb_dt += -k_CSF_PL_AbM_mAb_tran * CSF_AbM_mAb
    d_PL_AbM_mAb_dt += k_CSF_PL_AbM_mAb_tran * CSF_AbM_mAb
    # PL->CSF
    d_CSF_AbM_mAb_dt += k_PL_CSF_AbM_mAb_tran * PL_AbM_mAb
    d_PL_AbM_mAb_dt += -k_PL_CSF_AbM_mAb_tran * PL_AbM_mAb
    
    
    # Antibody AbF binding unbinding
    # ISF 
    # Fibril binding
    d_ISF_AbF_dt += -k_ISF_AbF_mAb_bind * ISF_AbF * ISF_mAb
    d_ISF_mAb_dt += -k_ISF_AbF_mAb_bind * ISF_AbF * ISF_mAb
    d_ISF_AbF_mAb_dt += k_ISF_AbF_mAb_bind * ISF_AbF * ISF_mAb
    # AbF-mAb unbinding
    d_ISF_AbF_dt += k_ISF_AbF_mAb_diss * ISF_AbF_mAb
    d_ISF_mAb_dt += k_ISF_AbF_mAb_diss * ISF_AbF_mAb
    d_ISF_AbF_mAb_dt += -k_ISF_AbF_mAb_diss * ISF_AbF_mAb
    
    
    # Antibody AbM binding unbinding
    
    # ISF
    # Binding
    d_ISF_AbM_dt += -k_ISF_AbM_mAb_bind * ISF_AbM * ISF_mAb
    d_ISF_mAb_dt += -k_ISF_AbM_mAb_bind * ISF_AbM * ISF_mAb
    d_ISF_AbM_mAb_dt += k_ISF_AbM_mAb_bind * ISF_AbM * ISF_mAb
    # Unbinding
    d_ISF_AbM_dt += k_ISF_AbM_mAb_diss * ISF_AbM_mAb
    d_ISF_mAb_dt += k_ISF_AbM_mAb_diss * ISF_AbM_mAb
    d_ISF_AbM_mAb_dt += -k_ISF_AbM_mAb_diss * ISF_AbM_mAb
    
    # CSF
    # binding
    d_CSF_AbM_dt += -k_CSF_AbM_mAb_bind * CSF_AbM * CSF_mAb
    d_CSF_mAb_dt += -k_CSF_AbM_mAb_bind * CSF_AbM * CSF_mAb
    d_CSF_AbM_mAb_dt += k_CSF_AbM_mAb_bind * CSF_AbM * CSF_mAb
    # unbinding
    d_CSF_AbM_dt += k_CSF_AbM_mAb_diss * CSF_AbM_mAb
    d_CSF_mAb_dt += k_CSF_AbM_mAb_diss * CSF_AbM_mAb
    d_CSF_AbM_mAb_dt += -k_CSF_AbM_mAb_diss * CSF_AbM_mAb
    
    # PL
    # binding
    d_PL_AbM_dt += -k_PL_AbM_mAb_bind * PL_AbM * PL_mAb
    d_PL_mAb_dt += -k_PL_AbM_mAb_bind * PL_AbM * PL_mAb
    d_PL_AbM_mAb_dt += k_PL_AbM_mAb_bind * PL_AbM * PL_mAb
    # unbinding
    d_PL_AbM_dt += k_PL_AbM_mAb_diss * PL_AbM_mAb
    d_PL_mAb_dt += k_PL_AbM_mAb_diss * PL_AbM_mAb
    d_PL_AbM_mAb_dt += -k_PL_AbM_mAb_diss * PL_AbM_mAb
    
    
    
    # AbFN mAb binding (fibril end binding), clearance
    # ISF
    # FN->FN1
    d_ISF_AbFN_dt += -2 * k_ISF_AbFN_mAb_bind * ISF_AbFN * ISF_mAb
    d_ISF_AbFN1_dt += 2 * k_ISF_AbFN_mAb_bind * ISF_AbFN * ISF_mAb
    d_ISF_mAb_dt += -2 * k_ISF_AbFN_mAb_bind * ISF_AbFN * ISF_mAb
    d_ISF_AbF_dt += -2 * k_ISF_AbFN_mAb_bind * ISF_AbF * ISF_mAb
    d_ISF_AbF1_dt += 2 * k_ISF_AbFN_mAb_bind * ISF_AbF * ISF_mAb
    # FN1->FN
    d_ISF_AbFN_dt += k_ISF_AbFN_mAb_diss * ISF_AbFN1
    d_ISF_AbFN1_dt += -k_ISF_AbFN_mAb_diss * ISF_AbFN1
    d_ISF_mAb_dt += k_ISF_AbFN_mAb_diss * ISF_AbFN1
    d_ISF_AbF_dt += k_ISF_AbFN_mAb_diss * ISF_AbF1
    d_ISF_AbF1_dt += -k_ISF_AbFN_mAb_diss * ISF_AbF1
    # FN1->FN2
    d_ISF_AbFN1_dt += -k_ISF_AbFN_mAb_bind * ISF_AbFN1 * ISF_mAb
    d_ISF_AbFN2_dt += k_ISF_AbFN_mAb_bind * ISF_AbFN1 * ISF_mAb
    d_ISF_mAb_dt += -k_ISF_AbFN_mAb_bind * ISF_AbFN1 * ISF_mAb
    d_ISF_AbF1_dt += -k_ISF_AbFN_mAb_bind * ISF_AbF1 * ISF_mAb
    d_ISF_AbF2_dt += k_ISF_AbFN_mAb_bind * ISF_AbF1 * ISF_mAb
    # FN2->FN1
    d_ISF_AbFN1_dt += 2 * k_ISF_AbFN_mAb_diss * ISF_AbFN2
    d_ISF_AbFN2_dt += - 2 * k_ISF_AbFN_mAb_diss * ISF_AbFN2
    d_ISF_mAb_dt += 2 * k_ISF_AbFN_mAb_diss * ISF_AbFN2
    d_ISF_AbF1_dt += 2 * k_ISF_AbFN_mAb_diss * ISF_AbF2
    d_ISF_AbF2_dt += - 2 * k_ISF_AbFN_mAb_diss * ISF_AbF2
    
    
    # Antibody 'production' in ISF, for testing
    d_ISF_mAb_dt += k_ISF_mAb_syn
    
    
    
    rhs.append(d_ISF_AbM_dt)
    rhs.append(d_ISF_AbON_dt)
    rhs.append(d_ISF_AbFN_dt)
    rhs.append(d_ISF_AbF_dt)
    rhs.append(d_ISF_AbP_dt)
    rhs.append(d_ISF_mAb_dt)
    rhs.append(d_ISF_AbF_mAb_dt)
    rhs.append(d_ISF_AbM_mAb_dt)
    rhs.append(d_ISF_AbFN1_dt)
    rhs.append(d_ISF_AbF1_dt)
    rhs.append(d_ISF_AbFN2_dt)
    rhs.append(d_ISF_AbF2_dt)
    
    rhs.append(d_PL_AbM_dt)
    rhs.append(d_PL_mAb_dt)
    rhs.append(d_PL_AbM_mAb_dt)
    
    rhs.append(d_CSF_AbM_dt)
    rhs.append(d_CSF_mAb_dt)
    rhs.append(d_CSF_AbM_mAb_dt)
    
    rhs.append(d_PE_mAb_dt)
    
    return rhs


###################################################################################
def set_basal_parameter_values():
    # This function sets all parameter values from different papers
    # Although the parameters are declared to be global, they are only 'global' within this module when imported.
    
    # Meanings of the suffixes:
    # syn: production or synchronization
    # clear: clearance
    # tran: transportation
    # bind: binding
    # diss: dissociation or breakage
    
    # Abbreviations of species names
    # AbM: Amyloid beta monomer
    # AbO: Amyloid beta oligomer
    # AbF: Amyloid beta fibril mass (the total number of Abeta monomers in fibrils)
    # AbFN: Amyloid beta fibril count (the number of fibrils)
    # AbP: Amyloid beta plaque
    # mAb: antibody
    # AbM_mAb: Abeta monomer bound by antibody
    # AbF_mAb: Abeta fibril bound by antibody (strictly speaking, the refers to Abeta monomers in fibrils that bound by antibodies, used when simulating the effect of fibril surface binding)
    # AbFN1: Abeta fibril with one end bound by an antibody
    # AbFN2: Abeta fibril with both ends bound by antibodies
    
    
    # Commonly used constants
    global m0, t0, week_to_hour, year_to_week
    global V_PL, MW_aducanumab, W_normal
    
    # Abeta production, aggregation
    global k_ISF_AbM_syn # Abeta monomer production rate
    global n_ISF_AbO_1stNuc, k_ISF_AbO_1stNuc # Primary nucleation order and rate
    global n_ISF_AbO # Number of monomers in an oligomer
    global k_ISF_AbO_diss # Oligmer dissociation rate
    global n_ISF_AbF_conv, k_ISF_AbF_conv # Conversion order and rate from oligomer to fibril
    global k_ISF_AbF_on # Monomer attachment rate to fibril ends
    global n_ISF_AbO_2ndNuc, k_ISF_AbO_2ndNuc # Secondary nucleation order and rate
    global k_ISF_AbP_conv, k_ISF_AbP_diss # Conversion rate from fibril to plaque, monomer dissociation rate from plaques
    
    # Clearance
    # ISF
    global k_ISF_AbM_clear, k_ISF_AbF_clear, k_ISF_AbP_clear # naming pattern: k_compartment_species_clear
    global k_ISF_mAb_clear
    global k_ISF_AbF_mAb_clear
    global k_ISF_AbM_mAb_clear 
    global k_ISF_AbFN1_clear 
    global k_ISF_AbF1_clear 
    global k_ISF_AbFN2_clear
    global k_ISF_AbF2_clear
    # PL
    global k_PL_AbM_clear
    global k_PL_mAb_clear
    global k_PL_AbM_mAb_clear
    
    
    # Transportation
    # Monomer transportation
    global k_ISF_CSF_AbM_tran 
    global k_ISF_PL_AbM_tran 
    global k_PL_ISF_AbM_tran
    global k_CSF_PL_AbM_tran
    global k_PL_CSF_AbM_tran
   
    # Antibody transportation
    global k_ISF_PL_mAb_tran
    global k_PL_ISF_mAb_tran
    global k_ISF_CSF_mAb_tran
    global k_CSF_PL_mAb_tran
    global k_PL_CSF_mAb_tran
    global k_PL_PE_mAb_tran
    global k_PE_PL_mAb_tran
    
    # AbM_mAb transportation rates
    global k_ISF_PL_AbM_mAb_tran
    global k_PL_ISF_AbM_mAb_tran
    global k_ISF_CSF_AbM_mAb_tran
    global k_CSF_PL_AbM_mAb_tran
    global k_PL_CSF_AbM_mAb_tran
    
    
    # Fibril surface-antibody binding unbinding
    global k_ISF_AbF_mAb_bind
    global k_ISF_AbF_mAb_diss
    # Monomer-antibody binding unbinding
    global k_ISF_AbM_mAb_bind
    global k_ISF_AbM_mAb_diss
    global k_CSF_AbM_mAb_bind
    global k_CSF_AbM_mAb_diss
    global k_PL_AbM_mAb_bind
    global k_PL_AbM_mAb_diss
    # Fibril end-antibody binding unbinding
    global k_ISF_AbFN_mAb_bind
    global k_ISF_AbFN_mAb_diss

        
    # antibody 'production' rate
    global k_ISF_mAb_syn
    
    
    
    # Commonly used constants
    m0 = 1e-9
    t0 = 3600
    week_to_hour = 7*24
    year_to_week = 52

    # Plasma volume
    V_PL = 3
    # Aducanumab molecular weight (g/mol)
    MW_aducanumab = 1.5e5
    # A normal person's weight (kg)
    W_normal = 65

    # Parameter values (adjusted)
    # All concentrations' unit is nM
    # Time's unit is hour
    
    # Abeta production, aggregation
    #k_ISF_AbM_syn = 1.2 # Raskatov2019
    k_ISF_AbM_syn = 7 # Estimated
    #k_ISF_AbM_syn = 8 # Testing

    n_ISF_AbO_1stNuc = 0.8
    k_ISF_AbO_1stNuc = 6.7e-8
    k_ISF_AbO_1stNuc = t0 * k_ISF_AbO_1stNuc * (m0**(n_ISF_AbO_1stNuc - 1))
    n_ISF_AbO = 2
    k_ISF_AbO_diss = 9.7e-5 * t0
    n_ISF_AbF_conv = 2.7
    k_ISF_AbF_conv = 1.9e9
    k_ISF_AbF_conv = t0 * k_ISF_AbF_conv * (m0**n_ISF_AbF_conv)
    k_ISF_AbF_on = 3e6
    k_ISF_AbF_on = t0 * k_ISF_AbF_on * m0
    n_ISF_AbO_2ndNuc = 0.9
    k_ISF_AbO_2ndNuc = 2
    k_ISF_AbO_2ndNuc = t0 * k_ISF_AbO_2ndNuc * (m0**n_ISF_AbO_2ndNuc)
    k_ISF_AbP_conv = 7e-8 * t0 # Lin2022 AbO->AbP conversion rate
    k_ISF_AbP_diss = 7e-11 * t0 # Lin2022 AbP->AbO conversion rate
    
    
    # Clearace
    # ISF
    #k_ISF_AbM_clear = 1.93e-5 * t0 # Lin2022
    #k_ISF_AbF_clear = 2.2e-8 * t0 # Lin2022
    #k_ISF_AbP_clear = 4.41e-9 * t0 # Lin2022
    k_ISF_AbM_clear = 1e-1 # Adjusted
    k_ISF_AbF_clear = 1e-3 # Adjusted
    #k_ISF_AbP_clear = 1e-4 # Adjusted
    k_ISF_AbP_clear = 1.5e-4 # Adjusted again
    k_ISF_mAb_clear = 0 # No such clearance in Lin or Madrasi, set to zero
    #k_ISF_AbF_mAb_clear = 2.2e-8 * t0 # Assume to be the same as AbF degradation rate as in Lin2022
    #k_ISF_AbF_mAb_clear = 4e-2 # Assumed to be higher than AbF clearance rate in ISF
    #k_ISF_AbM_mAb_clear = 4e-2 # Assumed to be higher than AbF clearance rate in ISF
    k_ISF_AbF_mAb_clear = 1e-3 # To match aducanumab data, this should be 4e-2
    k_ISF_AbM_mAb_clear = 1e-3
    k_ISF_AbFN1_clear = 1e-3
    k_ISF_AbF1_clear = 1e-3
    k_ISF_AbFN2_clear = 1e-3
    k_ISF_AbF2_clear = 1e-3

    # PL
    k_PL_AbM_clear = 9.63e-5 * t0 # Lin2022
    k_PL_mAb_clear = 1.38099820e-05 * t0 # Fitted
    k_PL_AbM_mAb_clear = 1.38099820e-05 * t0 # Fitted; higher than Lin2022 (1.46e-6)
    
    
    
    # Transportation
    # Monomer
    k_ISF_CSF_AbM_tran = 1.55e-5 * t0 # Lin2022
    k_ISF_PL_AbM_tran = 1.48e-5 * t0 # Lin2022     ????????
    k_PL_ISF_AbM_tran = 1.48e-4 * t0 # Lin2022    ???????
    k_CSF_PL_AbM_tran = 4.17e-5 * t0 # Lin2022   
    k_PL_CSF_AbM_tran = 1.72e-9 * t0 # Lin2022
    
    # Antibody
    #k_ISF_PL_mAb_tran = 3e-3 * t0  # Lin2022
    k_ISF_PL_mAb_tran = 3e-5 * t0  # Testing
    #k_PL_ISF_mAb_tran = 1.6e-6 * t0 # Lin2022
    k_PL_ISF_mAb_tran = 2e-6 * t0 # Testing
    k_ISF_CSF_mAb_tran = 1.55e-5 * t0 # Lin2022
    k_CSF_PL_mAb_tran = 4.17e-5 * t0 # Lin2022
    k_PL_CSF_mAb_tran = 1.72e-9 * t0 # Lin2022
    #k_PL_PE_mAb_tran = 2.5e-6 * t0 # Lin2022
    #k_PE_PL_mAb_tran = 1e-6 * t0 # Lin2022
    #k_PL_mAb_clear = 1.46e-6 * t0 # Lin2022
    k_PL_PE_mAb_tran = 1.81946225e-04 * t0 # Fitted
    k_PE_PL_mAb_tran = 8.56369151e-06 * t0 # Fitted
    
    # Monomer-antibody
    k_ISF_PL_AbM_mAb_tran = 3e-5 * t0  # Testing; Lin2022's value is 3e-3 * t0, too high
    k_PL_ISF_AbM_mAb_tran = 2e-6 * t0 # Testing; Lin2022's value is 1.6e-6 * t0, slightly less
    k_ISF_CSF_AbM_mAb_tran = 1.55e-5 * t0 # Lin2022
    k_CSF_PL_AbM_mAb_tran = 4.17e-5 * t0 # Lin2022  
    k_PL_CSF_AbM_mAb_tran = 1.72e-9 * t0 # Lin2022  
    
    
    # Binding and unbinding
    
    # Fibril surface-antibody
    # Estimation: KD = 1nM
    # Estimation binding rate: 1e-3/(nM*s)
    #k_ISF_AbF_mAb_bind = 1e-3 * t0
    k_ISF_AbF_mAb_bind = 2e-2 * t0 # This fits adu data better than 1e-3 * t0
    k_ISF_AbF_mAb_diss = 1e-3 * t0

    # Monomer-antibody 
    k_ISF_AbM_mAb_bind = 1e-3 * t0
    k_ISF_AbM_mAb_diss = 1e-3 * t0
    k_CSF_AbM_mAb_bind = 1e-3 * t0
    k_CSF_AbM_mAb_diss = 1e-3 * t0
    k_PL_AbM_mAb_bind = 1e-3 * t0
    k_PL_AbM_mAb_diss = 1e-3 * t0

    # Fibril end-antibody 
    k_ISF_AbFN_mAb_bind = 1e-3 * t0
    k_ISF_AbFN_mAb_diss = 1e-3 * t0
    
    # ISF antibody 'production' rate
    k_ISF_mAb_syn = 0
    
    
###################################################################################
def set_zero_initial_conditions():
    # This function set zero for all species
    y0 = []

    # ISF
    ISF_AbM = 0; y0.append(ISF_AbM)
    ISF_AbON = 0; y0.append(ISF_AbON);
    ISF_AbFN = 0; y0.append(ISF_AbFN)
    ISF_AbF = 0; y0.append(ISF_AbF)
    ISF_AbP = 0; y0.append(ISF_AbP)
    ISF_mAb = 0; y0.append(ISF_mAb)
    ISF_AbF_mAb = 0; y0.append(ISF_AbF_mAb)
    ISF_AbM_mAb = 0; y0.append(ISF_AbM_mAb)
    ISF_AbFN1 = 0; y0.append(ISF_AbFN1)
    ISF_AbF1 = 0; y0.append(ISF_AbF1)
    ISF_AbFN2 = 0; y0.append(ISF_AbFN2)
    ISF_AbF2 = 0; y0.append(ISF_AbF2)


    # PL
    PL_AbM = 0; y0.append(PL_AbM)
    PL_mAb = 0; y0.append(PL_mAb)
    PL_AbM_mAb = 0; y0.append(PL_AbM_mAb)


    # CSF
    CSF_AbM = 0; y0.append(CSF_AbM)
    CSF_mAb = 0; y0.append(CSF_mAb)
    CSF_AbM_mAb = 0; y0.append(CSF_AbM_mAb)

    # PE
    PE_mAb = 0; y0.append(PE_mAb)
    
    return y0

###################################################################################
def set_equilibrium_initial_conditions(sol):
    # This function extracts the values from the solution object 'sol'
    # 'sol' is obtained by simulating abeta evolution without antibody
    # Thus, the simulation must be performed before this function is called.
    global ISF_AbM_init 
    global ISF_AbON_init 
    global ISF_AbFN_init 
    global ISF_AbF_init 
    global ISF_AbP_init 
    global ISF_mAb_init 
    global ISF_AbF_mAb_init 
    global ISF_AbM_mAb_init 
    global ISF_AbFN1_init
    global ISF_AbF1_init 
    global ISF_AbFN2_init
    global ISF_AbF2_init 


    global PL_AbM_init
    global PL_mAb_init
    global PL_AbM_mAb_init


    global CSF_AbM_init
    global CSF_mAb_init
    global CSF_AbM_mAb_init

    global PE_mAb_init
    
    index = 0

    ISF_AbM_init = sol.y[index][-1]; index += 1
    ISF_AbON_init = sol.y[index][-1]; index += 1
    ISF_AbFN_init = sol.y[index][-1]; index += 1
    ISF_AbF_init = sol.y[index][-1]; index += 1
    ISF_AbP_init = sol.y[index][-1]; index += 1
    ISF_mAb_init = sol.y[index][-1]; index += 1
    ISF_AbF_mAb_init = sol.y[index][-1]; index += 1
    ISF_AbM_mAb_init = sol.y[index][-1]; index += 1
    ISF_AbFN1_init = sol.y[index][-1]; index += 1
    ISF_AbF1_init = sol.y[index][-1]; index += 1
    ISF_AbFN2_init = sol.y[index][-1]; index += 1
    ISF_AbF2_init = sol.y[index][-1]; index += 1


    PL_AbM_init = sol.y[index][-1]; index += 1
    PL_mAb_init = sol.y[index][-1]; index += 1
    PL_AbM_mAb_init = sol.y[index][-1]; index += 1


    CSF_AbM_init = sol.y[index][-1]; index += 1
    CSF_mAb_init = sol.y[index][-1]; index += 1
    CSF_AbM_mAb_init = sol.y[index][-1]; index += 1

    PE_mAb_init = sol.y[index][-1]; index += 1
    
    
    y0 = []

    # ISF
    y0.append(ISF_AbM_init)
    y0.append(ISF_AbON_init);
    y0.append(ISF_AbFN_init)
    y0.append(ISF_AbF_init)
    y0.append(ISF_AbP_init)
    y0.append(ISF_mAb_init)
    y0.append(ISF_AbF_mAb_init)
    y0.append(ISF_AbM_mAb_init)
    y0.append(ISF_AbFN1_init)
    y0.append(ISF_AbF1_init)
    y0.append(ISF_AbFN2_init)
    y0.append(ISF_AbF2_init)


    # PL
    y0.append(PL_AbM_init)
    y0.append(PL_mAb_init)
    y0.append(PL_AbM_mAb_init)


    # CSF
    y0.append(CSF_AbM_init)
    y0.append(CSF_mAb_init)
    y0.append(CSF_AbM_mAb_init)

    # PE
    y0.append(PE_mAb_init)
    
    return y0
    
############################################################################################
def get_hist(sol):
    # Get the time evolution for each species as well as all the time points.
    global t_hist
    

    global ISF_AbM_hist
    global ISF_AbON_hist
    global ISF_AbFN_hist
    global ISF_AbF_hist
    global ISF_AbP_hist
    global ISF_mAb_hist
    global ISF_AbF_mAb_hist
    global ISF_AbM_mAb_hist
    global ISF_AbFN1_hist
    global ISF_AbF1_hist
    global ISF_AbFN2_hist
    global ISF_AbF2_hist

    global PL_AbM_hist
    global PL_mAb_hist
    global PL_AbM_mAb_hist

    global CSF_AbM_hist
    global CSF_mAb_hist
    global CSF_AbM_mAb_hist

    global PE_mAb_hist
    
    t_hist = sol.t
    
    index = 0

    ISF_AbM_hist = sol.y[index]; index += 1
    ISF_AbON_hist = sol.y[index]; index += 1
    ISF_AbFN_hist = sol.y[index]; index += 1
    ISF_AbF_hist = sol.y[index]; index += 1
    ISF_AbP_hist = sol.y[index]; index += 1
    ISF_mAb_hist = sol.y[index]; index += 1
    ISF_AbF_mAb_hist = sol.y[index]; index += 1
    ISF_AbM_mAb_hist = sol.y[index]; index += 1
    ISF_AbFN1_hist = sol.y[index]; index += 1
    ISF_AbF1_hist = sol.y[index]; index += 1
    ISF_AbFN2_hist = sol.y[index]; index += 1
    ISF_AbF2_hist = sol.y[index]; index += 1

    PL_AbM_hist = sol.y[index]; index += 1
    PL_mAb_hist = sol.y[index]; index += 1
    PL_AbM_mAb_hist = sol.y[index]; index += 1

    CSF_AbM_hist = sol.y[index]; index += 1
    CSF_mAb_hist = sol.y[index]; index += 1
    CSF_AbM_mAb_hist = sol.y[index]; index += 1

    PE_mAb_hist = sol.y[index]; index += 1
    
    
#####################################################################################
def set_basal_parameter_values_Michaels(basal_concentration):
    # This function sets parameter values only for abeta aggregation in ISF
    # There are no AbM production, mAb entrance, clearance, binding, or transportation.
    # Thus, the result reflects the evolution of abeta species in ISF with a given initial
    # ...AbM, which is exactly Michaels2020's experiment. 
    
    # Commonly used constants
    global m0, t0, week_to_hour, year_to_week
    global V_PL, MW_aducanumab, W_normal
    
    # Abeta production, aggregation
    global k_ISF_AbM_syn
    global n_ISF_AbO_1stNuc, k_ISF_AbO_1stNuc
    global n_ISF_AbO
    global k_ISF_AbO_diss
    global n_ISF_AbF_conv, k_ISF_AbF_conv
    global k_ISF_AbF_on
    global n_ISF_AbO_2ndNuc, k_ISF_AbO_2ndNuc
    global k_ISF_AbP_conv, k_ISF_AbP_diss
    
    # Clearance
    # ISF
    global k_ISF_AbM_clear, k_ISF_AbF_clear, k_ISF_AbP_clear 
    global k_ISF_mAb_clear
    global k_ISF_AbF_mAb_clear
    global k_ISF_AbM_mAb_clear 
    global k_ISF_AbFN1_clear 
    global k_ISF_AbF1_clear 
    global k_ISF_AbFN2_clear
    global k_ISF_AbF2_clear
    # PL
    global k_PL_AbM_clear
    global k_PL_mAb_clear
    global k_PL_AbM_mAb_clear
    
    
    # Transportation
    # Monomer transportation
    global k_ISF_CSF_AbM_tran 
    global k_ISF_PL_AbM_tran 
    global k_PL_ISF_AbM_tran
    global k_CSF_PL_AbM_tran
    global k_PL_CSF_AbM_tran
   
    # Antibody transportation
    global k_ISF_PL_mAb_tran
    global k_PL_ISF_mAb_tran
    global k_ISF_CSF_mAb_tran
    global k_CSF_PL_mAb_tran
    global k_PL_CSF_mAb_tran
    global k_PL_PE_mAb_tran
    global k_PE_PL_mAb_tran
    
    # AbM_mAb transportation rates
    global k_ISF_PL_AbM_mAb_tran
    global k_PL_ISF_AbM_mAb_tran
    global k_ISF_CSF_AbM_mAb_tran
    global k_CSF_PL_AbM_mAb_tran
    global k_PL_CSF_AbM_mAb_tran
    
    
    # Fibril surface-antibody binding unbinding
    global k_ISF_AbF_mAb_bind
    global k_ISF_AbF_mAb_diss
    # Monomer-antibody binding unbinding
    global k_ISF_AbM_mAb_bind
    global k_ISF_AbM_mAb_diss
    global k_CSF_AbM_mAb_bind
    global k_CSF_AbM_mAb_diss
    global k_PL_AbM_mAb_bind
    global k_PL_AbM_mAb_diss
    # Fibril end-antibody binding unbinding
    global k_ISF_AbFN_mAb_bind
    global k_ISF_AbFN_mAb_diss

        
    # antibody 'production' rate
    global k_ISF_mAb_syn
    
    
    
    # Commonly used constants
    m0 = basal_concentration # Note that the rescaling concentration is different from ours (1e-9)
    t0 = 3600
    week_to_hour = 7*24
    year_to_week = 52

    # Plasma volume
    V_PL = 3
    # Aducanumab molecular weight (g/mol)
    MW_aducanumab = 1.5e5
    # A normal person's weight (kg)
    W_normal = 65

    # Parameter values (adjusted)
    # All concentrations' unit is nM
    # Time's unit is hour
    
    # Abeta production, aggregation
    #k_ISF_AbM_syn = 1.2 # Raskatov2019
    k_ISF_AbM_syn = 0

    n_ISF_AbO_1stNuc = 0.8
    k_ISF_AbO_1stNuc = 6.7e-8
    k_ISF_AbO_1stNuc = t0 * k_ISF_AbO_1stNuc * (m0**(n_ISF_AbO_1stNuc - 1))
    n_ISF_AbO = 2
    k_ISF_AbO_diss = 9.7e-5 * t0
    n_ISF_AbF_conv = 2.7
    k_ISF_AbF_conv = 1.9e9
    k_ISF_AbF_conv = t0 * k_ISF_AbF_conv * (m0**n_ISF_AbF_conv)
    k_ISF_AbF_on = 3e6
    k_ISF_AbF_on = t0 * k_ISF_AbF_on * m0
    n_ISF_AbO_2ndNuc = 0.9
    k_ISF_AbO_2ndNuc = 2
    k_ISF_AbO_2ndNuc = t0 * k_ISF_AbO_2ndNuc * (m0**n_ISF_AbO_2ndNuc)
    k_ISF_AbP_conv = 0
    k_ISF_AbP_diss = 0
    
    
    # Clearace
    # ISF
    k_ISF_AbM_clear = 0
    k_ISF_AbF_clear = 0
    k_ISF_AbP_clear = 0
    k_ISF_mAb_clear = 0 
    k_ISF_AbF_mAb_clear = 0
    k_ISF_AbM_mAb_clear = 0
    k_ISF_AbFN1_clear = 0
    k_ISF_AbF1_clear = 0
    k_ISF_AbFN2_clear = 0
    k_ISF_AbF2_clear = 0

    # PL
    k_PL_AbM_clear = 0
    k_PL_mAb_clear = 0
    k_PL_AbM_mAb_clear = 0
    
    
    
    # Transportation
    # Monomer
    k_ISF_CSF_AbM_tran = 0
    k_ISF_PL_AbM_tran = 0
    k_PL_ISF_AbM_tran = 0
    k_CSF_PL_AbM_tran = 0
    k_PL_CSF_AbM_tran = 0
    
    # Antibody
    k_ISF_PL_mAb_tran = 0
    k_PL_ISF_mAb_tran = 0
    k_ISF_CSF_mAb_tran = 0
    k_CSF_PL_mAb_tran = 0
    k_PL_CSF_mAb_tran = 0
    k_PL_PE_mAb_tran = 0
    k_PE_PL_mAb_tran = 0
    
    # Monomer-antibody
    k_ISF_PL_AbM_mAb_tran = 0
    k_PL_ISF_AbM_mAb_tran = 0
    k_ISF_CSF_AbM_mAb_tran = 0
    k_CSF_PL_AbM_mAb_tran = 0
    k_PL_CSF_AbM_mAb_tran = 0
    
    
    # Binding and unbinding
    k_ISF_AbF_mAb_bind = 0
    k_ISF_AbF_mAb_diss = 0

    # Monomer-antibody 
    k_ISF_AbM_mAb_bind = 0
    k_ISF_AbM_mAb_diss = 0
    k_CSF_AbM_mAb_bind = 0
    k_CSF_AbM_mAb_diss = 0
    k_PL_AbM_mAb_bind = 0
    k_PL_AbM_mAb_diss = 0

    # Fibril end-antibody 
    k_ISF_AbFN_mAb_bind = 0
    k_ISF_AbFN_mAb_diss = 0
    
    # ISF antibody 'production' rate
    k_ISF_mAb_syn = 0
    
    
###################################################################################
def show_parameters():
    # This function displays all parameters for diagnostic purpose.
    # Commonly used constants
    global m0, t0, week_to_hour, year_to_week
    global V_PL, MW_aducanumab, W_normal
    
    # Abeta production, aggregation
    global k_ISF_AbM_syn
    global n_ISF_AbO_1stNuc, k_ISF_AbO_1stNuc
    global n_ISF_AbO
    global k_ISF_AbO_diss
    global n_ISF_AbF_conv, k_ISF_AbF_conv
    global k_ISF_AbF_on
    global n_ISF_AbO_2ndNuc, k_ISF_AbO_2ndNuc
    global k_ISF_AbP_conv, k_ISF_AbP_diss
    
    # Clearance
    # ISF
    global k_ISF_AbM_clear, k_ISF_AbF_clear, k_ISF_AbP_clear 
    global k_ISF_mAb_clear
    global k_ISF_AbF_mAb_clear
    global k_ISF_AbM_mAb_clear 
    global k_ISF_AbFN1_clear 
    global k_ISF_AbF1_clear 
    global k_ISF_AbFN2_clear
    global k_ISF_AbF2_clear
    # PL
    global k_PL_AbM_clear
    global k_PL_mAb_clear
    global k_PL_AbM_mAb_clear
    
    
    # Transportation
    # Monomer transportation
    global k_ISF_CSF_AbM_tran 
    global k_ISF_PL_AbM_tran 
    global k_PL_ISF_AbM_tran
    global k_CSF_PL_AbM_tran
    global k_PL_CSF_AbM_tran
   
    # Antibody transportation
    global k_ISF_PL_mAb_tran
    global k_PL_ISF_mAb_tran
    global k_ISF_CSF_mAb_tran
    global k_CSF_PL_mAb_tran
    global k_PL_CSF_mAb_tran
    global k_PL_PE_mAb_tran
    global k_PE_PL_mAb_tran
    
    # AbM_mAb transportation rates
    global k_ISF_PL_AbM_mAb_tran
    global k_PL_ISF_AbM_mAb_tran
    global k_ISF_CSF_AbM_mAb_tran
    global k_CSF_PL_AbM_mAb_tran
    global k_PL_CSF_AbM_mAb_tran
    
    
    # Fibril surface-antibody binding unbinding
    global k_ISF_AbF_mAb_bind
    global k_ISF_AbF_mAb_diss
    # Monomer-antibody binding unbinding
    global k_ISF_AbM_mAb_bind
    global k_ISF_AbM_mAb_diss
    global k_CSF_AbM_mAb_bind
    global k_CSF_AbM_mAb_diss
    global k_PL_AbM_mAb_bind
    global k_PL_AbM_mAb_diss
    # Fibril end-antibody binding unbinding
    global k_ISF_AbFN_mAb_bind
    global k_ISF_AbFN_mAb_diss

        
    # antibody 'production' rate
    global k_ISF_mAb_syn
    
    
    
    # Commonly used constants
    print('m0', m0)
    print('t0', t0)
    print('week_to_hour', week_to_hour)
    print('year_to_week', year_to_week)

    # Plasma volume
    print('V_PL', V_PL)
    # Aducanumab molecular weight (g/mol)
    print('MW_aducanumab', MW_aducanumab)
    # A normal person's weight (kg)
    print('W_normal', W_normal)

    # Parameter values (adjusted)
    # All concentrations' unit is nM
    # Time's unit is hour
    
    # Abeta production, aggregation
    print("AbM production")
    print('k_ISF_AbM_syn', k_ISF_AbM_syn)
   
    print("Ab aggregation")
    print('n_ISF_AbO_1stNuc', n_ISF_AbO_1stNuc)
    print('k_ISF_AbO_1stNuc', k_ISF_AbO_1stNuc)
    print('n_ISF_AbO', n_ISF_AbO)
    print('k_ISF_AbO_diss', k_ISF_AbO_diss)
    print('n_ISF_AbF_conv', n_ISF_AbF_conv)
    print('k_ISF_AbF_conv', k_ISF_AbF_conv)
    print('k_ISF_AbF_on', k_ISF_AbF_on)
    print('n_ISF_AbO_2ndNuc', n_ISF_AbO_2ndNuc)
    print('k_ISF_AbO_2ndNuc', k_ISF_AbO_2ndNuc)
    print('k_ISF_AbP_conv', k_ISF_AbP_conv)
    print('k_ISF_AbP_diss', k_ISF_AbP_diss)
    
    
    
    # Clearace
    print("Clearance")
    # ISF
    print("ISF")
    print('k_ISF_AbM_clear', k_ISF_AbM_clear)
    print('k_ISF_AbF_clear', k_ISF_AbF_clear)
    print('k_ISF_AbP_clear', k_ISF_AbP_clear)
    print('k_ISF_mAb_clear', k_ISF_mAb_clear)
    print('k_ISF_AbF_mAb_clear', k_ISF_AbF_mAb_clear)
    print('k_ISF_AbM_mAb_clear', k_ISF_AbM_mAb_clear)
    print('k_ISF_AbFN1_clear', k_ISF_AbFN1_clear)
    print('k_ISF_AbF1_clear', k_ISF_AbF1_clear)
    print('k_ISF_AbFN2_clear', k_ISF_AbFN2_clear)
    print('k_ISF_AbF2_clear', k_ISF_AbF2_clear)

    # PL
    print("PL")
    print('k_PL_AbM_clear', k_PL_AbM_clear)
    print('k_PL_mAb_clear', k_PL_mAb_clear)
    print('k_PL_AbM_mAb_clear', k_PL_AbM_mAb_clear)
    
    
    
    # Transportation
    print("Transportation")
    # Monomer
    print("Monomer")
    print('k_ISF_CSF_AbM_tran', k_ISF_CSF_AbM_tran)
    print('k_ISF_PL_AbM_tran', k_ISF_PL_AbM_tran)
    print('k_PL_ISF_AbM_tran', k_PL_ISF_AbM_tran)
    print('k_CSF_PL_AbM_tran', k_CSF_PL_AbM_tran)
    print('k_PL_CSF_AbM_tran', k_PL_CSF_AbM_tran)
    
    # Antibody
    print("Antibody")
    print('k_ISF_PL_mAb_tran', k_ISF_PL_mAb_tran)
    print('k_PL_ISF_mAb_tran', k_PL_ISF_mAb_tran)
    print('k_ISF_CSF_mAb_tran', k_ISF_CSF_mAb_tran)
    print('k_CSF_PL_mAb_tran', k_CSF_PL_mAb_tran)
    print('k_PL_CSF_mAb_tran', k_PL_CSF_mAb_tran)
    print('k_PL_PE_mAb_tran', k_PL_PE_mAb_tran)
    print('k_PE_PL_mAb_tran', k_PE_PL_mAb_tran)
    
    # Monomer-antibody
    print("AbM-mAb")
    print('k_ISF_PL_AbM_mAb_tran', k_ISF_PL_AbM_mAb_tran)
    print('k_PL_ISF_AbM_mAb_tran', k_PL_ISF_AbM_mAb_tran)
    print('k_ISF_CSF_AbM_mAb_tran', k_ISF_CSF_AbM_mAb_tran)
    print('k_CSF_PL_AbM_mAb_tran', k_CSF_PL_AbM_mAb_tran)
    print('k_PL_CSF_AbM_mAb_tran', k_PL_CSF_AbM_mAb_tran)
    
    
    # Binding and unbinding
    print("Binding")
    
    # Fibril surface-antibody
    print("Fibril surface")
    print('k_ISF_AbF_mAb_bind', k_ISF_AbF_mAb_bind)
    print('k_ISF_AbF_mAb_diss', k_ISF_AbF_mAb_diss)

    # Monomer-antibody 
    print("Monomer")
    print('k_ISF_AbM_mAb_bind', k_ISF_AbM_mAb_bind)
    print('k_ISF_AbM_mAb_diss', k_ISF_AbM_mAb_diss)
    print('k_CSF_AbM_mAb_bind', k_CSF_AbM_mAb_bind)
    print('k_CSF_AbM_mAb_diss', k_CSF_AbM_mAb_diss)
    print('k_PL_AbM_mAb_bind', k_PL_AbM_mAb_bind)
    print('k_PL_AbM_mAb_diss', k_PL_AbM_mAb_diss)

    # Fibril end-antibody 
    print("Fibril end")
    print('k_ISF_AbFN_mAb_bind', k_ISF_AbFN_mAb_bind)
    print('k_ISF_AbFN_mAb_diss', k_ISF_AbFN_mAb_diss)
    
    # ISF antibody 'production' rate
    print("Antibody production for testing")
    print('k_ISF_mAb_syn', k_ISF_mAb_syn)
    
###############################################################################

def show_currents():
    # Show currents
    print("ISF Ab Aggregation")
    print("AbM production:", k_ISF_AbM_syn)
    print("Primary nucleation (AbON):", k_ISF_AbO_1stNuc * (ISF_AbM_hist[-1]**n_ISF_AbO_1stNuc))
    print("AbON dissociation:", k_ISF_AbO_diss * ISF_AbON_hist[-1])
    print("Oligomer conversion into fibril:", 
          k_ISF_AbF_conv * ISF_AbON_hist[-1] * (ISF_AbM_hist[-1]**n_ISF_AbF_conv))
    print("Fibril elongation:", 
          2 * k_ISF_AbF_on * ISF_AbM_hist[-1] * ISF_AbFN_hist[-1])
    ISF_AbF_total = ISF_AbF_hist[-1] + ISF_AbF1_hist[-1] + ISF_AbF2_hist[-1]
    print("Secondary nucleation:", 
          k_ISF_AbO_2ndNuc * (ISF_AbM_hist[-1]**n_ISF_AbO_2ndNuc) * ISF_AbF_total)

    print("ISF Clearance")
    print("Monomer clearance:", k_ISF_AbM_clear * ISF_AbM_hist[-1])
    print("Fibril clearance:", k_ISF_AbF_clear * ISF_AbF_hist[-1])
    print("Antibody clearance:", k_ISF_mAb_clear * ISF_mAb_hist[-1])
    print("AbM-mAb clearance:", k_ISF_AbM_mAb_clear * ISF_AbM_mAb_hist[-1])
    print("AbF-mAb clearance:", k_ISF_AbF_mAb_clear * ISF_AbF_mAb_hist[-1])
    print("Fibril end binding")
    print("AbFN1-mAb clearance:", k_ISF_AbFN1_clear * ISF_AbFN1_hist[-1])
    print("AbF1-mAb clearance:", k_ISF_AbF1_clear * ISF_AbF1_hist[-1])
    print("AbFN2-mAb clearance:", k_ISF_AbFN2_clear * ISF_AbFN2_hist[-1])
    print("AbF2-mAb clearance:", k_ISF_AbF2_clear * ISF_AbF2_hist[-1])
    
##################################################################################################
def set_equilibrium_initial_conditions(sol):
    # This should be run after simulating abeta evolution without antibody
    global ISF_AbM_init 
    global ISF_AbON_init 
    global ISF_AbFN_init 
    global ISF_AbF_init 
    global ISF_AbP_init 
    global ISF_mAb_init 
    global ISF_AbF_mAb_init 
    global ISF_AbM_mAb_init 
    global ISF_AbFN1_init
    global ISF_AbF1_init 
    global ISF_AbFN2_init
    global ISF_AbF2_init 


    global PL_AbM_init
    global PL_mAb_init
    global PL_AbM_mAb_init


    global CSF_AbM_init
    global CSF_mAb_init
    global CSF_AbM_mAb_init

    global PE_mAb_init
    
    index = 0

    ISF_AbM_init = sol.y[index][-1]; index += 1
    ISF_AbON_init = sol.y[index][-1]; index += 1
    ISF_AbFN_init = sol.y[index][-1]; index += 1
    ISF_AbF_init = sol.y[index][-1]; index += 1
    ISF_AbP_init = sol.y[index][-1]; index += 1
    ISF_mAb_init = sol.y[index][-1]; index += 1
    ISF_AbF_mAb_init = sol.y[index][-1]; index += 1
    ISF_AbM_mAb_init = sol.y[index][-1]; index += 1
    ISF_AbFN1_init = sol.y[index][-1]; index += 1
    ISF_AbF1_init = sol.y[index][-1]; index += 1
    ISF_AbFN2_init = sol.y[index][-1]; index += 1
    ISF_AbF2_init = sol.y[index][-1]; index += 1


    PL_AbM_init = sol.y[index][-1]; index += 1
    PL_mAb_init = sol.y[index][-1]; index += 1
    PL_AbM_mAb_init = sol.y[index][-1]; index += 1


    CSF_AbM_init = sol.y[index][-1]; index += 1
    CSF_mAb_init = sol.y[index][-1]; index += 1
    CSF_AbM_mAb_init = sol.y[index][-1]; index += 1

    PE_mAb_init = sol.y[index][-1]; index += 1
    
    
    y0 = []

    # ISF
    y0.append(ISF_AbM_init)
    y0.append(ISF_AbON_init);
    y0.append(ISF_AbFN_init)
    y0.append(ISF_AbF_init)
    y0.append(ISF_AbP_init)
    y0.append(ISF_mAb_init)
    y0.append(ISF_AbF_mAb_init)
    y0.append(ISF_AbM_mAb_init)
    y0.append(ISF_AbFN1_init)
    y0.append(ISF_AbF1_init)
    y0.append(ISF_AbFN2_init)
    y0.append(ISF_AbF2_init)


    # PL
    y0.append(PL_AbM_init)
    y0.append(PL_mAb_init)
    y0.append(PL_AbM_mAb_init)


    # CSF
    y0.append(CSF_AbM_init)
    y0.append(CSF_mAb_init)
    y0.append(CSF_AbM_mAb_init)

    # PE
    y0.append(PE_mAb_init)
    
    return y0
    
    
    
    
    
    