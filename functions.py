import numpy as np
from scipy.integrate import solve_ivp


# Create indices
def set_indices_species():
    indices = {}
    index = 0
    # Create name:index pair for all relevant species
    indices['ISF_AbM'] = index; index += 1
    indices['ISF_AbON'] = index; index += 1
    indices['ISF_AbFN'] = index; index += 1
    indices['ISF_AbF'] = index; index += 1
    indices['ISF_AbP'] = index; index += 1
    indices['ISF_mAb'] = index; index += 1
    indices['ISF_AbF_mAb'] = index; index += 1
    indices['ISF_AbM_mAb'] = index; index += 1
    indices['ISF_AbFN1'] = index; index += 1
    indices['ISF_AbF1'] = index; index += 1
    indices['ISF_AbFN2'] = index; index += 1
    indices['ISF_AbF2'] = index; index += 1
    
    # PL
    indices['PL_AbM'] = index; index += 1
    indices['PL_mAb'] = index; index += 1
    indices['PL_AbM_mAb'] = index; index += 1
   
    # CSF
    indices['CSF_AbM'] = index; index += 1
    indices['CSF_mAb'] = index; index += 1
    indices['CSF_AbM_mAb'] = index; index += 1
    
    # PE
    indices['PE_mAb'] = index; index += 1

    return [indices, index]

indices_species, num_species = set_indices_species()

#############################################################
def set_indices_parameters():
    indices = {}
    index = 0
    
    # Commonly used constants
    indices['m0'] = index; index += 1
    indices['t0'] = index; index += 1
    indices['week_to_hour'] = index; index += 1
    indices['year_to_week'] = index; index += 1
    
    # Plasma volume
    indices['V_PL'] = index; index += 1
    # Aducanumab molecular weight 
    indices['MW_aducanumab'] = index; index += 1
    # A normal person's weight 
    indices['W_normal'] = index; index += 1

    # Parameter values 
    
    # Abeta production, aggregation
    indices['k_ISF_AbM_syn'] = index; index += 1

    indices['n_ISF_AbO_1stNuc'] = index; index += 1
    indices['k_ISF_AbO_1stNuc'] = index; index += 1
    indices['n_ISF_AbO = 2'] = index; index += 1
    indices['k_ISF_AbO_diss'] = index; index += 1
    indices['n_ISF_AbF_conv'] = index; index += 1
    indices['k_ISF_AbF_conv'] = index; index += 1
    indices['k_ISF_AbF_on'] = index; index += 1
    indices['n_ISF_AbO_2ndNuc'] = index; index += 1
    indices['k_ISF_AbO_2ndNuc'] = index; index += 1  
    indices['k_ISF_AbP_conv'] = index; index += 1
    indices['k_ISF_AbP_diss'] = index; index += 1
    
    # Clearace
    # ISF
    indices['k_ISF_AbM_clear'] = index; index += 1
    indices['k_ISF_AbF_clear'] = index; index += 1
    indices['k_ISF_AbP_clear'] = index; index += 1
    indices['k_ISF_mAb_clear'] = index; index += 1
    indices['k_ISF_AbF_mAb_clear'] = index; index += 1
    indices['k_ISF_AbM_mAb_clear'] = index; index += 1   
    indices['k_ISF_AbFN1_clear'] = index; index += 1    
    indices['k_ISF_AbF1_clear'] = index; index += 1    
    indices['k_ISF_AbFN2_clear'] = index; index += 1    
    indices['k_ISF_AbF2_clear'] = index; index += 1
    
    # PL
    indices['k_PL_AbM_clear'] = index; index += 1    
    indices['k_PL_mAb_clear'] = index; index += 1    
    indices['k_PL_AbM_mAb_clear'] = index; index += 1
    
    # Transportation
    # Monomer
    indices['k_ISF_CSF_AbM_tran'] = index; index += 1
    indices['k_ISF_PL_AbM_tran'] = index; index += 1
    indices['k_PL_ISF_AbM_tran'] = index; index += 1    
    indices['k_CSF_PL_AbM_tran'] = index; index += 1    
    indices['k_PL_CSF_AbM_tran'] = index; index += 1
    
    
    # Antibody
    indices['k_ISF_PL_mAb_tran'] = index; index += 1    
    indices['k_PL_ISF_mAb_tran'] = index; index += 1    
    indices['k_ISF_CSF_mAb_tran'] = index; index += 1    
    indices['k_CSF_PL_mAb_tran'] = index; index += 1    
    indices['k_PL_CSF_mAb_tran'] = index; index += 1    
    indices['k_PL_PE_mAb_tran'] = index; index += 1    
    indices['k_PE_PL_mAb_tran'] = index; index += 1
    
    
    # Monomer-antibody
    indices['k_ISF_PL_AbM_mAb_tran'] = index; index += 1
    indices['k_PL_ISF_AbM_mAb_tran'] = index; index += 1
    indices['k_ISF_CSF_AbM_mAb_tran'] = index; index += 1
    indices['k_CSF_PL_AbM_mAb_tran'] = index; index += 1
    indices['k_PL_CSF_AbM_mAb_tran'] = index; index += 1
    
    # Binding and unbinding
    indices['k_ISF_AbF_mAb_bind'] = index; index += 1
    indices['k_ISF_AbF_mAb_diss'] = index; index += 1
    

    # Monomer-antibody 
    indices['k_ISF_AbM_mAb_bind'] = index; index += 1
    indices['k_ISF_AbM_mAb_diss'] = index; index += 1        
    indices['k_CSF_AbM_mAb_bind'] = index; index += 1    
    indices['k_CSF_AbM_mAb_diss'] = index; index += 1    
    indices['k_PL_AbM_mAb_bind'] = index; index += 1    
    indices['k_PL_AbM_mAb_diss'] = index; index += 1
    
    # Fibril end-antibody
    indices['k_ISF_AbFN_mAb_bind'] = index; index += 1
    indices['k_ISF_AbFN_mAb_diss'] = index; index += 1
    
    # ISF antibody 'production' rate
    indices['k_ISF_mAb_syn'] = index; index += 1
        
    return indices


###################################################################################
def RHS(t, y):
    
    # rhs stores all the righ-hand sides to be returned
    rhs = np.zeros(num_species)
    # index of a variable
    index = 0
    
    ## Variables
    
    idx = indices_species
    
    # ISF
    ISF_AbM = y[idx['ISF_AbM']]
    ISF_AbON = y[idx['ISF_AbON']]; 
    ISF_AbFN = y[idx['ISF_AbFN']]; 
    ISF_AbF = y[idx['ISF_AbF']]; 
    ISF_AbP = y[idx['ISF_AbP']]; 
    ISF_mAb = y[idx['ISF_mAb']]; 
    ISF_AbF_mAb = y[idx['ISF_AbF_mAb']]; 
    ISF_AbM_mAb = y[idx['ISF_AbM_mAb']]; 
    ISF_AbFN1 = y[idx['ISF_AbFN1']];  # one antibody bound to a fibril end
    ISF_AbF1 = y[idx['ISF_AbF1']]; 
    ISF_AbFN2 = y[idx['ISF_AbFN2']];  # two antibody bound to two fibril ends
    ISF_AbF2 = y[idx['ISF_AbF2']]; 
    
    ISF_AbFN_total = ISF_AbFN + ISF_AbFN1 + ISF_AbFN2
    ISF_AbF_total = ISF_AbF + ISF_AbF1 + ISF_AbF2
    
    
    # PL
    PL_AbM = y[idx['PL_AbM']]; 
    PL_mAb = y[idx['PL_mAb']]; 
    PL_AbM_mAb = y[idx['PL_AbM_mAb']];  
    
    
    # CSF
    CSF_AbM = y[idx['CSF_AbM']];     
    CSF_mAb = y[idx['CSF_mAb']]; 
    CSF_AbM_mAb = y[idx['CSF_AbM_mAb']]; 
    
    # PE
    PE_mAb = y[idx['PE_mAb']]; 
    
    
    
    p = parameters
    
    
    ## Calculating the right-hand sides

    # ISF
    # AbM production
    rhs[idx['ISF_AbM']] += p['k_ISF_AbM_syn']
    # Primary nucleation
    rhs[idx['ISF_AbM']] += - p['k_ISF_AbO_1stNuc'] * (ISF_AbM**p['n_ISF_AbO_1stNuc']) * p['n_ISF_AbO']
    rhs[idx['ISF_AbON']] += p['k_ISF_AbO_1stNuc'] * (ISF_AbM**p['n_ISF_AbO_1stNuc'])
    # Oligomer dissociation
    rhs[idx['ISF_AbON']] += - p['k_ISF_AbO_diss'] * ISF_AbON
    rhs[idx['ISF_AbM']] += p['k_ISF_AbO_diss'] * ISF_AbON * p['n_ISF_AbO']
    # Oligomer conversion into fibril
    rhs[idx['ISF_AbON']] += - p['k_ISF_AbF_conv'] * ISF_AbON * (ISF_AbM**p['n_ISF_AbF_conv'])
    rhs[idx['ISF_AbFN']] += p['k_ISF_AbF_conv'] * ISF_AbON * (ISF_AbM**p['n_ISF_AbF_conv'])
    rhs[idx['ISF_AbF']] += p['k_ISF_AbF_conv'] * ISF_AbON * (ISF_AbM**p['n_ISF_AbF_conv']) * p['n_ISF_AbO']
    # Fibril elongation
    rhs[idx['ISF_AbF']] += 2 * p['k_ISF_AbF_on'] * ISF_AbM * ISF_AbFN
    rhs[idx['ISF_AbM']] += - 2 * p['k_ISF_AbF_on'] * ISF_AbM * ISF_AbFN
    rhs[idx['ISF_AbF1']] += p['k_ISF_AbF_on'] * ISF_AbM * ISF_AbFN1
    rhs[idx['ISF_AbM']] += - p['k_ISF_AbF_on'] * ISF_AbM * ISF_AbFN1
    # Secondary nucleation
    rhs[idx['ISF_AbM']] += -p['k_ISF_AbO_2ndNuc'] * (ISF_AbM**p['n_ISF_AbO_2ndNuc']) * ISF_AbF_total * p['n_ISF_AbO']
    rhs[idx['ISF_AbON']] += p['k_ISF_AbO_2ndNuc'] * (ISF_AbM**p['n_ISF_AbO_2ndNuc']) * ISF_AbF_total
    # Fibril conversion into plaque
    rhs[idx['ISF_AbFN']] += - p['k_ISF_AbP_conv'] * ISF_AbFN   
    rhs[idx['ISF_AbF']] += - p['k_ISF_AbP_conv'] * ISF_AbF   
    rhs[idx['ISF_AbP']] += p['k_ISF_AbP_conv'] * ISF_AbF    
    rhs[idx['ISF_AbFN1']] += - p['k_ISF_AbP_conv'] * ISF_AbFN1 # Fibril end binding might affect plaque formation    
    rhs[idx['ISF_AbF1']] += - p['k_ISF_AbP_conv'] * ISF_AbF1    
    rhs[idx['ISF_AbP']] += p['k_ISF_AbP_conv'] * ISF_AbF1    
    rhs[idx['ISF_AbFN2']] += - p['k_ISF_AbP_conv'] * ISF_AbFN2   
    rhs[idx['ISF_AbF2']] += - p['k_ISF_AbP_conv'] * ISF_AbF2   
    rhs[idx['ISF_AbP']] += p['k_ISF_AbP_conv'] * ISF_AbF2
    # Plaque dissociation into monomer    
    rhs[idx['ISF_AbM']] += p['k_ISF_AbP_diss'] * ISF_AbP    
    rhs[idx['ISF_AbP']] += - p['k_ISF_AbP_diss'] * ISF_AbP
    
    
    
    # Clearance of all species, all compartments
    
    # ISF
    # Monomer clearance
    rhs[idx['ISF_AbM']] += - p['k_ISF_AbM_clear'] * ISF_AbM
    # Fibril clearance    
    rhs[idx['ISF_AbF']] += - p['k_ISF_AbF_clear'] * ISF_AbF    
    rhs[idx['ISF_AbFN']] += - p['k_ISF_AbF_clear'] * ISF_AbFN
    # Plaque clearance    
    rhs[idx['ISF_AbP']] += - p['k_ISF_AbP_clear'] * ISF_AbP
    # Antibody clearance    
    rhs[idx['ISF_mAb']] += - p['k_ISF_mAb_clear'] * ISF_mAb # No such clearance in both Lin and Madrasi
    # Monomer-antibody clearance    
    rhs[idx['ISF_AbM_mAb']] += -p['k_ISF_AbM_mAb_clear'] * ISF_AbM_mAb
    # Fibril surface-antibody clearance    
    rhs[idx['ISF_AbF_mAb']] += -p['k_ISF_AbF_mAb_clear'] * ISF_AbF_mAb
    # Fibril end-antibody clearance   
    rhs[idx['ISF_AbFN1']] += -p['k_ISF_AbFN1_clear'] * ISF_AbFN1    
    rhs[idx['ISF_AbF1']] += -p['k_ISF_AbF1_clear'] * ISF_AbF1
    rhs[idx['ISF_AbFN2']] += -p['k_ISF_AbFN2_clear'] * ISF_AbFN2    
    rhs[idx['ISF_AbF2']] += -p['k_ISF_AbF2_clear'] * ISF_AbF2
    
    # PL
    # Monomer clearance
    rhs[idx['PL_AbM']] += -p['k_PL_AbM_clear'] * PL_AbM
    # Antibody clearance    
    rhs[idx['PL_mAb']] += -p['k_PL_mAb_clear'] * PL_mAb
    # Monomer-antibody clearance    
    rhs[idx['PL_AbM_mAb']] += -p['k_PL_AbM_mAb_clear'] * PL_AbM_mAb
    
    
    # Abeta monomer transportation
    # ISF -> CSF   
    rhs[idx['ISF_AbM']] += -p['k_ISF_CSF_AbM_tran'] * ISF_AbM   
    rhs[idx['CSF_AbM']] += p['k_ISF_CSF_AbM_tran'] * ISF_AbM
    # ISF -> PL
    rhs[idx['ISF_AbM']] += -p['k_ISF_PL_AbM_tran'] * ISF_AbM
    rhs[idx['PL_AbM']] +=p['k_ISF_PL_AbM_tran'] * ISF_AbM
    # PL -> ISF    
    rhs[idx['ISF_AbM']] += p['k_PL_ISF_AbM_tran'] * PL_AbM    
    rhs[idx['PL_AbM']] += -p['k_PL_ISF_AbM_tran'] * PL_AbM
    # CSF -> PL
    rhs[idx['CSF_AbM']] += -p['k_CSF_PL_AbM_tran'] * CSF_AbM
    rhs[idx['PL_AbM']] += p['k_CSF_PL_AbM_tran'] * CSF_AbM
    # PL -> CSF    
    rhs[idx['CSF_AbM']] += p['k_PL_CSF_AbM_tran'] * PL_AbM    
    rhs[idx['PL_AbM']] += -p['k_PL_CSF_AbM_tran'] * PL_AbM
    
    
    
    # Antibody transportation
    # ISF -> PL    
    rhs[idx['ISF_mAb']] += -p['k_ISF_PL_mAb_tran'] * ISF_mAb   
    rhs[idx['PL_mAb']] += p['k_ISF_PL_mAb_tran'] * ISF_mAb
    # PL -> ISF
    rhs[idx['ISF_mAb']] += p['k_PL_ISF_mAb_tran'] * PL_mAb   
    rhs[idx['PL_mAb']] += -p['k_PL_ISF_mAb_tran'] * PL_mAb
    # ISF -> CSF
    rhs[idx['ISF_mAb']] += -p['k_ISF_CSF_mAb_tran'] * ISF_mAb
    rhs[idx['CSF_mAb']] += p['k_ISF_CSF_mAb_tran'] * ISF_mAb
    # CSF -> PL   
    rhs[idx['CSF_mAb']] += -p['k_CSF_PL_mAb_tran'] * CSF_mAb  
    rhs[idx['PL_mAb']] += p['k_CSF_PL_mAb_tran'] * CSF_mAb
    # PL -> CSF
    rhs[idx['CSF_mAb']] += p['k_PL_CSF_mAb_tran'] * PL_mAb
    rhs[idx['PL_mAb']] += -p['k_PL_CSF_mAb_tran'] * PL_mAb
    # PL -> PE
    rhs[idx['PL_mAb']] += -p['k_PL_PE_mAb_tran'] * PL_mAb
    rhs[idx['PE_mAb']] += p['k_PL_PE_mAb_tran'] * PL_mAb
    # PE -> PL
    rhs[idx['PL_mAb']] += p['k_PE_PL_mAb_tran'] * PE_mAb
    rhs[idx['PE_mAb']] += -p['k_PE_PL_mAb_tran'] * PE_mAb
    
    
    # mAb-AbM transportation
    # ISF->PL 
    rhs[idx['ISF_AbM_mAb']] += -p['k_ISF_PL_AbM_mAb_tran'] * ISF_AbM_mAb
    rhs[idx['PL_AbM_mAb']] += p['k_ISF_PL_AbM_mAb_tran'] * ISF_AbM_mAb
    # PL->ISF
    rhs[idx['ISF_AbM_mAb']] += p['k_PL_ISF_AbM_mAb_tran'] * PL_AbM_mAb
    rhs[idx['PL_AbM_mAb']] += -p['k_PL_ISF_AbM_mAb_tran'] * PL_AbM_mAb
    # ISF->CSF
    rhs[idx['ISF_AbM_mAb']] += -p['k_ISF_CSF_AbM_mAb_tran'] * ISF_AbM_mAb
    rhs[idx['CSF_AbM_mAb']] += p['k_ISF_CSF_AbM_mAb_tran'] * ISF_AbM_mAb
    # CSF->PL
    rhs[idx['CSF_AbM_mAb']] += -p['k_CSF_PL_AbM_mAb_tran'] * CSF_AbM_mAb
    rhs[idx['PL_AbM_mAb']] += p['k_CSF_PL_AbM_mAb_tran'] * CSF_AbM_mAb
    # PL->CSF
    rhs[idx['CSF_AbM_mAb']] += p['k_PL_CSF_AbM_mAb_tran'] * PL_AbM_mAb
    rhs[idx['PL_AbM_mAb']] += -p['k_PL_CSF_AbM_mAb_tran'] * PL_AbM_mAb
    
    
    # Antibody AbF binding unbinding
    # ISF 
    # Fibril binding    
    rhs[idx['ISF_AbF']] += -p['k_ISF_AbF_mAb_bind'] * ISF_AbF * ISF_mAb    
    rhs[idx['ISF_mAb']] += -p['k_ISF_AbF_mAb_bind'] * ISF_AbF * ISF_mAb    
    rhs[idx['ISF_AbF_mAb']] += p['k_ISF_AbF_mAb_bind'] * ISF_AbF * ISF_mAb
    # AbF-mAb unbinding   
    rhs[idx['ISF_AbF']] += p['k_ISF_AbF_mAb_diss'] * ISF_AbF_mAb   
    rhs[idx['ISF_mAb']] += p['k_ISF_AbF_mAb_diss'] * ISF_AbF_mAb   
    rhs[idx['ISF_AbF_mAb']] += -p['k_ISF_AbF_mAb_diss'] * ISF_AbF_mAb
    
    
    # Antibody AbM binding unbinding
    
    # ISF
    # Binding
    rhs[idx['ISF_AbM']] += -p['k_ISF_AbM_mAb_bind'] * ISF_AbM * ISF_mAb
    rhs[idx['ISF_mAb']] += -p['k_ISF_AbM_mAb_bind'] * ISF_AbM * ISF_mAb
    rhs[idx['ISF_AbM_mAb']] += p['k_ISF_AbM_mAb_bind'] * ISF_AbM * ISF_mAb
    # Unbinding
    rhs[idx['ISF_AbM']] += p['k_ISF_AbM_mAb_diss'] * ISF_AbM_mAb
    rhs[idx['ISF_mAb']] += p['k_ISF_AbM_mAb_diss'] * ISF_AbM_mAb
    rhs[idx['ISF_AbM_mAb']] += -p['k_ISF_AbM_mAb_diss'] * ISF_AbM_mAb
    
    # CSF
    # binding
    rhs[idx['CSF_AbM']] += -p['k_CSF_AbM_mAb_bind'] * CSF_AbM * CSF_mAb
    rhs[idx['CSF_mAb']] += -p['k_CSF_AbM_mAb_bind'] * CSF_AbM * CSF_mAb
    rhs[idx['CSF_AbM_mAb']] += p['k_CSF_AbM_mAb_bind'] * CSF_AbM * CSF_mAb
    # unbinding
    rhs[idx['CSF_AbM']] += p['k_CSF_AbM_mAb_diss'] * CSF_AbM_mAb
    rhs[idx['CSF_mAb']] += p['k_CSF_AbM_mAb_diss'] * CSF_AbM_mAb
    rhs[idx['CSF_AbM_mAb']] += -p['k_CSF_AbM_mAb_diss'] * CSF_AbM_mAb
    
    # PL
    # binding
    rhs[idx['PL_AbM']] += -p['k_PL_AbM_mAb_bind'] * PL_AbM * PL_mAb
    rhs[idx['PL_mAb']] += -p['k_PL_AbM_mAb_bind'] * PL_AbM * PL_mAb
    rhs[idx['PL_AbM_mAb']] += p['k_PL_AbM_mAb_bind'] * PL_AbM * PL_mAb
    # unbinding
    rhs[idx['PL_AbM']] += p['k_PL_AbM_mAb_diss'] * PL_AbM_mAb
    rhs[idx['PL_mAb']] += p['k_PL_AbM_mAb_diss'] * PL_AbM_mAb
    rhs[idx['PL_AbM_mAb']] += -p['k_PL_AbM_mAb_diss'] * PL_AbM_mAb
    
    
    
    # AbFN mAb binding (fibril end binding), clearance
    # ISF
    # FN->FN1
    rhs[idx['ISF_AbFN']] += -2 * p['k_ISF_AbFN_mAb_bind'] * ISF_AbFN * ISF_mAb
    rhs[idx['ISF_AbFN1']] += 2 * p['k_ISF_AbFN_mAb_bind'] * ISF_AbFN * ISF_mAb
    rhs[idx['ISF_mAb']] += -2 * p['k_ISF_AbFN_mAb_bind'] * ISF_AbFN * ISF_mAb
    rhs[idx['ISF_AbF']] += -2 * p['k_ISF_AbFN_mAb_bind'] * ISF_AbF * ISF_mAb
    rhs[idx['ISF_AbF1']] += 2 * p['k_ISF_AbFN_mAb_bind'] * ISF_AbF * ISF_mAb
    # FN1->FN
    rhs[idx['ISF_AbFN']] += p['k_ISF_AbFN_mAb_diss'] * ISF_AbFN1
    rhs[idx['ISF_AbFN1']] += - p['k_ISF_AbFN_mAb_diss'] * ISF_AbFN1
    rhs[idx['ISF_mAb']] += p['k_ISF_AbFN_mAb_diss'] * ISF_AbFN1
    rhs[idx['ISF_AbF']] += p['k_ISF_AbFN_mAb_diss'] * ISF_AbF1
    rhs[idx['ISF_AbF1']] += - p['k_ISF_AbFN_mAb_diss'] * ISF_AbF1
    # FN1->FN2
    rhs[idx['ISF_AbFN1']] += - p['k_ISF_AbFN_mAb_bind'] * ISF_AbFN1 * ISF_mAb
    rhs[idx['ISF_AbFN2']] += p['k_ISF_AbFN_mAb_bind'] * ISF_AbFN1 * ISF_mAb
    rhs[idx['ISF_mAb']] += - p['k_ISF_AbFN_mAb_bind'] * ISF_AbFN1 * ISF_mAb
    rhs[idx['ISF_AbF1']] += - p['k_ISF_AbFN_mAb_bind'] * ISF_AbF1 * ISF_mAb
    rhs[idx['ISF_AbF2']] += p['k_ISF_AbFN_mAb_bind'] * ISF_AbF1 * ISF_mAb
    # FN2->FN1
    rhs[idx['ISF_AbFN1']] += 2 * p['k_ISF_AbFN_mAb_diss'] * ISF_AbFN2
    rhs[idx['ISF_AbFN2']] += - 2 * p['k_ISF_AbFN_mAb_diss'] * ISF_AbFN2
    rhs[idx['ISF_mAb']] += 2 * p['k_ISF_AbFN_mAb_diss'] * ISF_AbFN2
    rhs[idx['ISF_AbF1']] += 2 * p['k_ISF_AbFN_mAb_diss'] * ISF_AbF2
    rhs[idx['ISF_AbF2']] += - 2 * p['k_ISF_AbFN_mAb_diss'] * ISF_AbF2
    
    
    # Antibody 'production' in ISF, for testing
    rhs[idx['ISF_mAb']] += p['k_ISF_mAb_syn']
    
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

    
    pars = {}
    
    # Commonly used constants
    pars['m0'] = 1e-9; m0 = pars['m0'] # set the variable m0 for notational convenience
    pars['t0'] = 3600; t0 = pars['t0'] # set the variable t0 for notational convenience
    pars['week_to_hour'] = 7*24
    pars['year_to_week'] = 52

    # Plasma volume
    pars['V_PL'] = 3
    # Aducanumab molecular weight (g/mol)
    pars['MW_aducanumab'] = 1.5e5
    # A normal person's weight (kg)
    pars['W_normal'] = 65

    # Parameter values (adjusted)
    # All concentrations' unit is nM
    # Time's unit is hour
    
    # Abeta production, aggregation
    #k_ISF_AbM_syn = 1.2 # Raskatov2019
    pars['k_ISF_AbM_syn'] = 7
    pars['n_ISF_AbO_1stNuc'] = 0.8
    
    k_ISF_AbO_1stNuc = 6.7e-8 
    k_ISF_AbO_1stNuc = t0 * k_ISF_AbO_1stNuc * (m0**(pars['n_ISF_AbO_1stNuc'] - 1))
    pars['k_ISF_AbO_1stNuc'] = k_ISF_AbO_1stNuc
    
    pars['n_ISF_AbO'] = 2
    pars['k_ISF_AbO_diss'] = 9.7e-5 * t0
    pars['n_ISF_AbF_conv'] = 2.7
    
    k_ISF_AbF_conv = 1.9e9
    k_ISF_AbF_conv = t0 * k_ISF_AbF_conv * (m0**pars['n_ISF_AbF_conv'])
    pars['k_ISF_AbF_conv'] = k_ISF_AbF_conv
    
    k_ISF_AbF_on = 3e6
    k_ISF_AbF_on = t0 * k_ISF_AbF_on * m0
    pars['k_ISF_AbF_on'] = k_ISF_AbF_on
    
    pars['n_ISF_AbO_2ndNuc'] = 0.9
    
    k_ISF_AbO_2ndNuc = 2
    k_ISF_AbO_2ndNuc = t0 * k_ISF_AbO_2ndNuc * (m0**pars['n_ISF_AbO_2ndNuc'])
    pars['k_ISF_AbO_2ndNuc'] = k_ISF_AbO_2ndNuc
    
    pars['k_ISF_AbP_conv'] = 7e-8 * t0 # Lin2022 AbO->AbP conversion rate
    pars['k_ISF_AbP_diss'] = 7e-11 * t0 # Lin2022 AbP->AbO conversion rate
   
    
    
    # Clearace
    # ISF
    #k_ISF_AbM_clear = 1.93e-5 * t0 # Lin2022
    #k_ISF_AbF_clear = 2.2e-8 * t0 # Lin2022
    #k_ISF_AbP_clear = 4.41e-9 * t0 # Lin2022
    pars['k_ISF_AbM_clear'] = 1e-1 # Adjusted
    pars['k_ISF_AbF_clear'] = 1e-3 # Adjusted
    pars['k_ISF_AbP_clear'] =1.5e-4
    pars['k_ISF_mAb_clear'] =0
    
    #k_ISF_AbF_mAb_clear = 2.2e-8 * t0 # Assume to be the same as AbF degradation rate as in Lin2022
    #k_ISF_AbF_mAb_clear = 4e-2 # Assumed to be higher than AbF clearance rate in ISF
    #k_ISF_AbM_mAb_clear = 4e-2 # Assumed to be higher than AbF clearance rate in ISF
    pars['k_ISF_AbF_mAb_clear'] = 1e-3 # To match aducanumab data, this should be 4e-2
    pars['k_ISF_AbM_mAb_clear'] = 1e-3
    pars['k_ISF_AbFN1_clear'] = 1e-3
    pars['k_ISF_AbF1_clear'] = 1e-3
    pars['k_ISF_AbFN2_clear'] = 1e-3
    pars['k_ISF_AbF2_clear'] = 1e-3

    # PL
    pars['k_PL_AbM_clear'] = 9.63e-5 * t0 # Lin2022
    pars['k_PL_mAb_clear'] = 1.38099820e-05 * t0 # Fitted
    pars['k_PL_AbM_mAb_clear'] = 1.38099820e-05 * t0 # Fitted; higher than Lin2022 (1.46e-6)
    
    
    
    # Transportation
    # Monomer
    pars['k_ISF_CSF_AbM_tran'] = 1.55e-5 * t0 # Lin2022
    pars['k_ISF_PL_AbM_tran'] = 1.48e-5 * t0 # Lin2022     ????????
    pars['k_PL_ISF_AbM_tran'] = 1.48e-4 * t0 # Lin2022    ???????
    pars['k_CSF_PL_AbM_tran'] = 4.17e-5 * t0 # Lin2022   
    pars['k_PL_CSF_AbM_tran'] = 1.72e-9 * t0 # Lin2022
   
    
    # Antibody
    #k_ISF_PL_mAb_tran = 3e-3 * t0  # Lin2022
    pars['k_ISF_PL_mAb_tran'] = 3e-5 * t0  # Testing   
    #k_PL_ISF_mAb_tran = 1.6e-6 * t0 # Lin2022
    pars['k_PL_ISF_mAb_tran'] = 2e-6 * t0 # Testing    
    pars['k_ISF_CSF_mAb_tran'] = 1.55e-5 * t0 # Lin2022   
    pars['k_CSF_PL_mAb_tran'] = 4.17e-5 * t0 # Lin2022
    pars['k_PL_CSF_mAb_tran'] = 1.72e-9 * t0 # Lin2022   
    #k_PL_PE_mAb_tran = 2.5e-6 * t0 # Lin2022
    #k_PE_PL_mAb_tran = 1e-6 * t0 # Lin2022
    #k_PL_mAb_clear = 1.46e-6 * t0 # Lin2022
    pars['k_PL_PE_mAb_tran'] = 1.81946225e-04 * t0 # Fitted
    pars['k_PE_PL_mAb_tran'] = 8.56369151e-06 * t0 # Fitted
    
    
    # Monomer-antibody
    pars['k_ISF_PL_AbM_mAb_tran'] = 3e-5 * t0  # Testing; Lin2022's value is 3e-3 * t0, too high
    pars['k_PL_ISF_AbM_mAb_tran'] = 2e-6 * t0 # Testing; Lin2022's value is 1.6e-6 * t0, slightly less
    pars['k_ISF_CSF_AbM_mAb_tran'] = 1.55e-5 * t0 # Lin2022
    pars['k_CSF_PL_AbM_mAb_tran'] = 4.17e-5 * t0 # Lin2022  
    pars['k_PL_CSF_AbM_mAb_tran'] = 1.72e-9 * t0 # Lin2022  
   
    
    
    # Binding and unbinding
    
    # Fibril surface-antibody
    # Estimation: KD = 1nM
    # Estimation binding rate: 1e-3/(nM*s)
    #k_ISF_AbF_mAb_bind = 2e-2 * t0 # This fits adu data better than 1e-3 * t0
    pars['k_ISF_AbF_mAb_bind'] = 1e-3 * t0
    pars['k_ISF_AbF_mAb_diss'] = 1e-3 * t0

    # Monomer-antibody 
    pars['k_ISF_AbM_mAb_bind'] = 1e-3 * t0
    pars['k_ISF_AbM_mAb_diss'] = 1e-3 * t0
    pars['k_CSF_AbM_mAb_bind'] = 1e-3 * t0
    pars['k_CSF_AbM_mAb_diss'] = 1e-3 * t0
    pars['k_PL_AbM_mAb_bind'] = 1e-3 * t0
    pars['k_PL_AbM_mAb_diss'] = 1e-3 * t0
    

    # Fibril end-antibody 
    pars['k_ISF_AbFN_mAb_bind'] = 1e-3 * t0
    pars['k_ISF_AbFN_mAb_diss'] = 1e-3 * t0
    
    # ISF antibody 'production' rate
    pars['k_ISF_mAb_syn'] = 0
    
    return pars

    
###################################################################################
def set_zero_initial_conditions():
    return np.zeros(num_species)

###################################################################################
def set_equilibrium_initial_conditions(sol):
    return sol.y[:, -1]    
    
#####################################################################################
def set_basal_parameter_values_Michaels(basal_concentration):
    # This function sets parameter values only for abeta aggregation in ISF
    # There are no AbM production, mAb entrance, clearance, binding, or transportation.
    # Thus, the result reflects the evolution of abeta species in ISF with a given initial
    # ...AbM, which is exactly Michaels2020's experiment. 
    
    set_basal_parameter_values()
    parameters['m0'] = basal_concentration
    parameters['k_ISF_AbP_conv'] = 0
    parameters['k_ISF_AbP_diss'] = 0
    
    
    # Clearace
    # ISF
    parameters['k_ISF_AbM_clear'] = 0
    parameters['k_ISF_AbF_clear'] = 0
    parameters['k_ISF_AbP_clear'] = 0
    parameters['k_ISF_mAb_clear'] = 0
    parameters['k_ISF_AbF_mAb_clear'] = 0
    parameters['k_ISF_AbM_mAb_clear'] = 0
    parameters['k_ISF_AbFN1_clear'] = 0
    parameters['k_ISF_AbF1_clear'] = 0
    parameters['k_ISF_AbFN2_clear'] = 0
    parameters['k_ISF_AbF2_clear'] = 0

    # PL
    parameters['k_PL_AbM_clear'] = 0
    parameters['k_PL_mAb_clear'] = 0
    parameters['k_PL_AbM_mAb_clear'] = 0

    
    # Transportation
    # Monomer
    parameters['k_ISF_CSF_AbM_tran'] = 0
    parameters['k_ISF_PL_AbM_tran'] = 0 
    parameters['k_PL_ISF_AbM_tran'] = 0  
    parameters['k_CSF_PL_AbM_tran'] = 0  
    parameters['k_PL_CSF_AbM_tran'] = 0
  
    
    # Antibody
    parameters['k_ISF_PL_mAb_tran'] = 0  
    parameters['k_PL_ISF_mAb_tran'] = 0
    parameters['k_ISF_CSF_mAb_tran'] = 0 
    parameters['k_CSF_PL_mAb_tran'] = 0 
    parameters['k_PL_CSF_mAb_tran'] = 0 
    parameters['k_PL_PE_mAb_tran'] = 0  
    parameters['k_PE_PL_mAb_tran'] = 0
 
    
    # Monomer-antibody
    parameters['k_ISF_PL_AbM_mAb_tran'] = 0 
    parameters['k_PL_ISF_AbM_mAb_tran'] = 0 
    parameters['k_ISF_CSF_AbM_mAb_tran'] = 0 
    parameters['k_CSF_PL_AbM_mAb_tran'] = 0  
    parameters['k_PL_CSF_AbM_mAb_tran'] = 0
    
    # Binding and unbinding
    parameters['k_ISF_AbF_mAb_bind'] = 0 
    parameters['k_ISF_AbF_mAb_diss'] = 0
  

    # Monomer-antibody 
    parameters['k_ISF_AbM_mAb_bind'] = 0 
    parameters['k_ISF_AbM_mAb_diss'] = 0   
    parameters['k_CSF_AbM_mAb_bind'] = 0    
    parameters['k_CSF_AbM_mAb_diss'] = 0    
    parameters['k_PL_AbM_mAb_bind'] = 0    
    parameters['k_PL_AbM_mAb_diss'] = 0

    # Fibril end-antibody 
    parameters['k_ISF_AbFN_mAb_bind'] = 0    
    parameters['k_ISF_AbFN_mAb_diss'] = 0
    
    # ISF antibody 'production' rate
    parameters['k_ISF_mAb_syn'] = 0
    
    
    
###################################################################################

    



################################################################
################################################################
################################################################
################################################################
################################################################

# Functions for Aducanumab PK fitting

def AduPKFitting_func_obj(pars, data):
    # Object function for fitting aducanumab PK
    parameters['k_PL_PE_mAb_tran'] = pars[0]
    parameters['k_PE_PL_mAb_tran'] = pars[1]
    parameters['k_PL_mAb_clear'] = pars[2]
    
    
    # data: Aducanumab PK data in Lin2022Fig2a
    time60, concentration60, time30, concentration30, time20, concentration20, time10, concentration10, time3, concentration3, time1, concentration1, time0p3, concentration0p3 = data
    
    
    
    
    # Initialize the loss
    L = 0
    
    dose = 60 * 1e-3
    for t, c in zip(time60, concentration60):
        # Prediction
        c_pred = AduPKFitting_cal_c_pred(t, dose)
        # Normalized squared difference
        L += (c_pred-c)**2 / (c**2 + 1e-8)
        
    
    dose = 30 * 1e-3 # dose unit: g/kg, the original unit is mg/kg
    for t, c in zip(time30, concentration30):
        # Prediction
        c_pred = AduPKFitting_cal_c_pred(t, dose)
        # Normalized squared difference
        L += (c_pred-c)**2 / (c**2 + 1e-8)
        
    dose = 20 * 1e-3
    for t, c in zip(time20, concentration20):
        # Prediction
        c_pred = AduPKFitting_cal_c_pred(t, dose)
        # Normalized squared difference
        L += (c_pred-c)**2 / (c**2 + 1e-8)
        
    dose = 10 * 1e-3
    for t, c in zip(time10, concentration10):
        # Prediction
        c_pred = AduPKFitting_cal_c_pred(t, dose)
        # Normalized squared difference
        L += (c_pred-c)**2 / (c**2 + 1e-8)
        
    dose = 3 * 1e-3
    for t, c in zip(time3, concentration3):
        # Prediction
        c_pred = AduPKFitting_cal_c_pred(t, dose)
        # Normalized squared difference
        L += (c_pred-c)**2 / (c**2 + 1e-8)
    
    dose = 1 * 1e-3
    for t, c in zip(time1, concentration1):
        # Prediction
        c_pred = AduPKFitting_cal_c_pred(t, dose)
        # Normalized squared difference
        L += (c_pred-c)**2 / (c**2 + 1e-8)
        
    
    dose = 0.3 * 1e-3
    for t, c in zip(time0p3, concentration0p3):
        # Prediction
        c_pred = AduPKFitting_cal_c_pred(t, dose)
        # Normalized squared difference
        L += (c_pred-c)**2 / (c**2 + 1e-8)
    
    return L

###########################################################
def AduPKFitting_cal_c_pred(t, dose):
    # Predict the PL mAb at t with dose
    
    # The initial plasma mAb concentration due to intravenous dosing
    # dose: g/kg
    # W_normal: kg
    # MW_aducanumab: g/mol
    # V_PL: liter
    # m0: 1e-9 mol, rescale the concentration to nm/liter
    
    p = parameters
    
    PL_mAb_init = dose * p['W_normal']/p['MW_aducanumab']/p['V_PL']/p['m0']
    
    # Initial condition
    y0 = set_zero_initial_conditions_except_PL_mAb(PL_mAb_init)
    
    # Solving the system up to time t
    sol = solve_ivp(fun=RHS, t_span=[0, t], y0=y0, 
                    method='LSODA', rtol=1e-3, max_step=1e2)
    # Get the current PL_mAb value
    c_pred = sol.y[indices_species['PL_mAb']][-1]
    
    return c_pred


###################################################################################
def set_zero_initial_conditions_except_PL_mAb(PL_mAb_init):
    # This function set zero for all species, except PL_mAb
    # This function is used for calculating PK given initial PL_mAb
        
    y0 = set_zero_initial_conditions()
    y0[indices_species['PL_mAb']] = PL_mAb_init
    
    return y0


#######################################################
def AduPKFitting_turnoff_binding():
    # Set all binding and dissociation rates to be zeros
    # Used for solving PK 
    
    # Fibril surface-antibody binding unbinding
    parameters['k_ISF_AbF_mAb_bind'] = 0
    parameters['k_ISF_AbF_mAb_diss'] = 0
    # Monomer-antibody binding unbinding
    parameters['k_ISF_AbM_mAb_bind'] = 0
    parameters['k_ISF_AbM_mAb_diss'] = 0
    parameters['k_CSF_AbM_mAb_bind'] = 0
    parameters['k_CSF_AbM_mAb_diss'] = 0
    parameters['k_PL_AbM_mAb_bind'] = 0
    parameters['k_PL_AbM_mAb_diss'] = 0
    # Fibril end-antibody binding unbinding
    parameters['k_ISF_AbFN_mAb_bind'] = 0
    parameters['k_ISF_AbFN_mAb_diss'] = 0
    
    
    
    
    
indices_species, num_species = set_indices_species()
indices_parameters = set_indices_parameters()
parameters = set_basal_parameter_values()
    
    
   

    
    
    
    
    
    
    
