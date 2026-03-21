# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 21:13:17 2019
"""

materialproperty ={
##########################################################################################
##########################################################################################
#          ZINCBLENDE
##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
#          I N S U L A T O R S    A N D    M E T A L S
##########################################################################################
##########################################################################################

######### air #############################################

'Air':{                                                   
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'IV_IV'   
    ,'lattice_consts':{
        'a'           :5.5
        ,'a_expansion' :0.0
    }
    
    ,'dielectric_consts':{
        'static_a'  :1.0                                        # vacuum value
        ,'optical_a' :1.0                                        # vacuum value
    }

    ,'elastic_consts':{
        'c11' :0.01 ,'c12' :0.01 ,'c44' :0.01                       # zero elastic energy
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :0.0                                              # no piezoelectric effect
    }                                                
                                                                      
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.156                        # Si value
          ,'bandgap'               :15.0                         # high barrier
          ,'bandgap_alpha'         :0
          ,'bandgap_beta'          :0
          ,'defpot_absolute'       :0
       }
       ,'L':{ 
          'mass_l'                :1.420                        # Si value
          ,'mass_t'                :0.130                        # Si value
          ,'bandgap'               :15.0                         # high barrier
          ,'bandgap_alpha'         :0
          ,'bandgap_beta'          :0
          ,'defpot_absolute'       :0
          ,'defpot_uniaxial'       :0
       }
       ,'Delta':{ 
          'mass_l'                :0.5                          # SiO2 value
          ,'mass_t'                :0.5                          # SiO2 value
          ,'bandgap'               :15.0                         # high barrier
          ,'bandgap_alpha'         :0
          ,'bandgap_beta'          :0
          ,'defpot_absolute'       :0
          ,'defpot_uniaxial'       :0
          ,'position'              :0.85                         # Si value
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-5.0                               # high barrier
        
        ,'HH':{ 'mass'          :0.537 }                            # Si value
        ,'LH':{ 'mass'          :0.153 }                            # Si value
        ,'SO':{ 'mass'          :0.234 }                            # Si value
        
        ,'defpot_absolute'   :0
        ,'defpot_uniaxial_b' :0   ,'defpot_uniaxial_d' :0
       
        ,'delta_SO'          :0.044                              # Si value
       
    }

    ,'kp_6_bands':{
        'L' :-6.69   ,'M' :-4.62   ,'N' :-8.56   ,'kappa' :-0.42      # Si values
    }
 
    ,'kp_8_bands':{                                                # No useful model for indirect ,'bandgap' materials!
        'S' :6.41                                               # inverse effective mass
        ,'E_P' :0                                                # decoupled bands
        ,'B' :0                                                  # This value is zero in group IV semiconductors (inversion symmetry).
        ,'L' :-6.69   ,'M' :-4.62   ,'N' :-8.56   ,'kappa' :-0.42      # 6-band parameters
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :142    ,'exponent' :2.5  }           # 1/10 constant Si value 
        ,'holes':{      'mumax' :47     ,'exponent' :2.2  }           # 1/10 constant Si value 
    }

    ,'mobility_masetti':{
        'electrons':{  'mumax'  :142   ,'exponent' :2.5              # 1/10 constant Si value 
                    ,'mumin1' :0     ,'mumin2'   :0       ,'mu1' :0  # dummy
                    ,'pc'     :0     ,'cr'       :1e100    ,'cs' :0   # dummy
                    ,'alpha'  :1     ,'beta'     :1                # dummy
        }
        ,'holes':{      'mumax'  :47    ,'exponent' :2.2              # 1/10 constant Si value
                    ,'mumin1' :0     ,'mumin2'   :0       ,'mu1' :0  # dummy
                    ,'pc'     :0     ,'cr'       :1e100   ,'cs' :0   # dummy
                    ,'alpha'  :1     ,'beta'     :1                # dummy
        }
    }

    ,'mobility_arora':{
        'electrons':{  'mumin' :142    ,'alm' :-2.5                  # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'mud' :0        ,'ald' :0                     # dummy
                    ,'n0'  :1e20     ,'aln' :1                     # dummy
                    ,'a'     :1      ,'ala' :1                     # dummy
        }
        ,'holes':{      'mumin' :47     ,'alm' :-2.2                  # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'mud' :0        ,'ald' :0                     # dummy
                    ,'n0'  :1e20     ,'aln' :0                     # dummy
                    ,'a'   :1        ,'ala' :0                     # dummy
        }
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :142       ,'muLexpT'  :-2.5     # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'muLImin300' :0         ,'TSwitch'    :200    # dummy
                    ,'muLIexpTabove' :0      ,'muLIexpTbelow' :0   # dummy
                    ,'Cref300'    :1e100     ,'CrefexpT'      :0   # dummy
                    ,'alpha300'   :1         ,'alphaexpT'     :0   # dummy
        }
        ,'holes':{      'muL300'     :47        ,'muLexpT'  :-2.2     # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'muLImin300' :0         ,'TSwitch'  :200      # dummy
                    ,'muLIexpTabove' :0      ,'muLIexpTbelow' :0   # dummy
                    ,'Cref300'    :1e100     ,'CrefexpT'      :0   # dummy
                    ,'alpha300'   :1         ,'alphaexpT'     :0   # dummy
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :4.26e-4     ,'nref_n' :7.1e15         # Si values
                   ,'tau_p' :3.95e-4     ,'nref_p' :7.1e15         # Si values
        }
        ,'Auger':{     'c_n' :2.8e-31       ,'c_p' :9.9e-31           # Si values
        }
    }
}  
#


######### silicon dioxide #############################################

#######################################################################
# E_gap is 9 eV.
# Conduction band offset SiO2/Si :3.1 eV
# Conduction band offset SiO2/Si :3.2 eV (M. Fischetti, JAP 83, 270 (1998))
# Note: SiO2 is hexagonal (wurtzite) and not cubic (zincblende)!!
#######################################################################

,'SiO2':{                                                   
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'IV_IV'                                            # Si value
   
    ,'lattice_consts':{
        'a'           :5.4304                                   # Si value
        ,'a_expansion' :1.8138e-5                                # Si value
    }
    
    ,'dielectric_consts':{
        'static_a'  :3.9                                        # PhD thesis R. Oberhuber
      # 'static_a'  :4.34                                       # www.crystran.co.uk/qutzdata.htm  4.34 (parallel) 4.27 (perpendicular) at 30MHz
        ,'optical_a' :3.5                                        # guess
    }

    ,'elastic_consts':{
      # 'c11' :87.0   ,'c12' :7.0   ,'c13' :13.0   ,'c33' :18.0   ,'c44' :58.0  # wurtzite, www.crystran.co.uk/qutzdata.htm
        'c11' :87.0   ,'c12' :7.0                             ,'c44' :58.0
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :0.000                                            # Si value
    }                                                
                                                                      
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.156                        # Si value
          ,'bandgap'               :11.25                        # guess (0 K)
          ,'bandgap_alpha'         :0
          ,'bandgap_beta'          :0
          ,'defpot_absolute'       :0
       }
       ,'L':{ 
          'mass_l'                :1.420                        # Si value
          ,'mass_t'                :0.130                        # Si value
          ,'bandgap'               :9.5                          # guess (0 K)
          ,'bandgap_alpha'         :0
          ,'bandgap_beta'          :0
          ,'defpot_absolute'       :0
          ,'defpot_uniaxial'       :0
       }
       ,'Delta':{ 
          'mass_l'                :0.5                          # M. Fischetti, JAP 83, 270 (1998)
          ,'mass_t'                :0.5                          # M. Fischetti, JAP 83, 270 (1998)
          ,'bandgap'               :9.0                          # (0 K)
          ,'bandgap_alpha'         :0
          ,'bandgap_beta'          :0
          ,'defpot_absolute'       :0
          ,'defpot_uniaxial'       :0
          ,'position'              :0.85                         # Si value
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-3.66                              # conduction band offset SiO2/Si :3.1 eV
                                                               # conduction band offset SiO2/Si :3.2 eV (M. Fischetti, JAP 83, 270 (1998))
        
        ,'HH':{ 'mass'          :0.537 }                            # Si value
        ,'LH':{ 'mass'          :0.153 }                            # Si value
        ,'SO':{ 'mass'          :0.234 }                            # Si value
        
        ,'defpot_absolute'   :0
        ,'defpot_uniaxial_b' :0   ,'defpot_uniaxial_d' :0
       
        ,'delta_SO'          :0.044                              # Si value
       
    }

    ,'kp_6_bands':{
        'L' :-6.69   ,'M' :-4.62   ,'N' :-8.56  ,'kappa' :-0.42       # Si values
    }
 
    ,'kp_8_bands':{                                                # No useful model for indirect ,'bandgap' materials!
        'S' :6.41                                               # inverse effective mass
        ,'E_P' :0                                                # decoupled bands
        ,'B' :0                                                  # This value is zero in group IV semiconductors (inversion symmetry).
        ,'L' :-6.69   ,'M' :-4.62   ,'N' :-8.56  ,'kappa' :-0.42       # 6-band parameters
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :142    ,'exponent' :2.5  }           # 1/10 constant Si value 
        ,'holes':{      'mumax' :47     ,'exponent' :2.2  }           # 1/10 constant Si value 
    }

    ,'mobility_masetti':{
        'electrons':{  'mumax'  :142   ,'exponent' :2.5              # 1/10 constant Si value 
                    ,'mumin1' :0     ,'mumin2'   :0       ,'mu1' :0  # dummy
                    ,'pc'     :0     ,'cr'       :1e100   ,'cs' :0   # dummy
                    ,'alpha'  :1     ,'beta'     :1                # dummy
        }
        ,'holes':{      'mumax'  :47    ,'exponent' :2.2              # 1/10 constant Si value
                    ,'mumin1' :0     ,'mumin2'   :0       ,'mu1' :0  # dummy
                    ,'pc'     :0     ,'cr'       :1e100   ,'cs' :0   # dummy
                    ,'alpha'  :1     ,'beta'     :1                # dummy
        }
    }

    ,'mobility_arora':{
        'electrons':{  'mumin' :142    ,'alm' :-2.5                  # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'mud' :0        ,'ald' :0                     # dummy
                    ,'n0'  :1e20     ,'aln' :1                     # dummy
                    ,'a'     :1      ,'ala' :1                     # dummy
        }
        ,'holes':{      'mumin' :47     ,'alm' :-2.2                  # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'mud' :0        ,'ald' :0                     # dummy
                    ,'n0'  :1e20     ,'aln' :0                     # dummy
                    ,'a'   :1        ,'ala' :0                     # dummy
        }
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :142       ,'muLexpT'  :-2.5     # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'muLImin300' :0         ,'TSwitch'    :200    # dummy
                    ,'muLIexpTabove' :0      ,'muLIexpTbelow' :0   # dummy
                    ,'Cref300'    :1e100     ,'CrefexpT'      :0   # dummy
                    ,'alpha300'   :1         ,'alphaexpT'     :0   # dummy
        }
        ,'holes':{      'muL300'     :47        ,'muLexpT'  :-2.2     # 1/10 constant Si value, but opposite ,'exponent' sign
                    ,'muLImin300' :0         ,'TSwitch'  :200      # dummy
                    ,'muLIexpTabove' :0      ,'muLIexpTbelow' :0   # dummy
                    ,'Cref300'    :1e100     ,'CrefexpT'      :0   # dummy
                    ,'alpha300'   :1         ,'alphaexpT'     :0   # dummy
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :4.26e-4     ,'nref_n' :7.1e15         # Si values
                   ,'tau_p' :3.95e-4     ,'nref_p' :7.1e15         # Si values
        }
        ,'Auger':{     'c_n' :2.8e-31       ,'c_p' :9.9e-31           # Si values
        }
    }
}  
#




##########################################################################################
##########################################################################################
#          B I N A R I E S     --    IV - IV       V A L E N C E
##########################################################################################
##########################################################################################




######### diamond ###################################################
,'C':{                                                   
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'IV_IV'   

    ,'lattice_consts':{
        'a'           :3.56683                                  # (300 K)
        ,'a_expansion' :0                                        # ???
    }
      
    ,'dielectric_consts':{
        'static_a'  :5.68                                       # http://www.kobelco.co.jp/showroom/np0802e/np08022e.htm
        ,'optical_a' :4                                          # ???
    }

    ,'elastic_consts':{ 
        'c11' :1079.0   ,'c12' :124.0   ,'c44' :578.0               #
      # 'c11' :1076.4   ,'c12' :125.2   ,'c44' :577.4               # Landolt-Boernstein, 298 K
    }     

    ,'piezoelectric_consts':{
        'e14' :0                                                # Piezoelectricity only occurs in III-V materials but not in group IV.
    }
  
# band gap 5.47 at 300 K (indirect)   5.46-5.6  E(Gamma)=7.3-7.4
    ,'conduction_bands':{
        'Gamma':{ 
          'mass'                  :1.9
          ,'bandgap'               :5.47                         # 300 K
          ,'bandgap_alpha'         :2.33e-4                      # 
          ,'bandgap_beta'          :1000                         # 
          ,'defpot_absolute'       :-10.41                       # ??? Ge value
        }
        ,'L':{ 
          'mass_l'                :1.57                         # ??? Ge value
          ,'mass_t'                :0.0807                       # ??? Ge value
          ,'bandgap'               :7                            # ???
          ,'bandgap_alpha'         :2.33e-4                      # 
          ,'bandgap_beta'          :1000                         # 
          ,'defpot_absolute'       :-4.35                        # ??? Ge value
          ,'defpot_uniaxial'       :15.13                        # ??? Ge value
        }
        ,'Delta':{
          'mass_l'                :1.40                         #  85   K 
          ,'mass_t'                :0.36                         #  85   K 
          ,'bandgap'               :7                            # ???
          ,'bandgap_alpha'         :2.33e-4                      # 
          ,'bandgap_beta'          :1000                         # 
          ,'defpot_absolute'       :0.14                         # ??? Ge value
          ,'defpot_uniaxial'       :9.42                         # ??? Ge value
          ,'position'              :0.85                         # ??? Ge value 0.85 for DELTA instead of 1.0 for X valley
          ,'g_l'                   :0.82                         # ??? Ge value
          ,'g_t'                   :1.93                         # ??? Ge value
        }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.09                               # ??? Si value
                                   
        ,'HH':{ 'mass'          :2.18 }                             # Landolt-Boernstein, cyclotron resonance at 1.2 K, value along [111]
      # ,'HH':{ 'mass'          :2.12 }                             # http://www.ioffe.ru/SVA/NSM/Semicond/Diamond/bandstr.html
        ,'LH':{ 'mass'          :0.70 }                             # Landolt-Boernstein, cyclotron resonance at 1.2 K, value along [111]
      # ,'LH':{ 'mass'          :0.7  }                             # http://www.ioffe.ru/SVA/NSM/Semicond/Diamond/bandstr.html
        ,'SO':{ 'mass'          :1.06 }                             # http://www.ioffe.ru/SVA/NSM/Semicond/Diamond/bandstr.html, Landolt-Boernstein, cyclotron resonance at 1.2 K, value along [111]
                                   
        ,'defpot_absolute'   :-0.35                              # ??? Ge value
        ,'defpot_uniaxial_b' :-2.86   ,'defpot_uniaxial_d' :-5.28  # ??? Ge value

        ,'delta_SO'          :0.006                              #

    }
    ,'kp_6_bands':{
      # gamma1 :2.54   gamma2 :-0.10   gamma3 :0.63         # M. Willatzen, M. Cardona, N.E. Christensen PRB 50, 18054 (1994)
        'L' :-3.140      ,'M' :-3.740       ,'N' :-3.780            # 
        ,'kappa' :-0.63                                          # P. Lawaetz, PRB 4, 3460 (1971)
    }                                   
                                       
    ,'kp_8_bands':{                                                # No useful model for indirect ,'bandgap' materials!
        'S' :1                                                  # ??? inverse effective mass
        ,'E_P' :49.8                                             # P. Lawaetz, PRB 4, 3460 (1971)
        ,'B' :0                                                  # This value is zero in group IV semiconductors (inversion symmetry).
        ,'L' :-3.140      ,'M' :-3.740       ,'N' :-3.780            # 6-band parameters
        ,'kappa' :-0.63                                          # 6-band parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :3800     ,'exponent' :1.66 }         # ??? Ge value
        ,'holes':{      'mumax' :1800     ,'exponent' :2.33 }         # ??? Ge value
    }

}
#

######### silicon #####################################################
,'Si':{                                                   
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'IV_IV'   
   
    ,'lattice_consts':{
        'a'           :5.4304                                   # Landolt-Boernstein  298.15 K
        ,'a_expansion' :1.8138e-5                                # Landolt-Boernstein
    }
    
    ,'dielectric_consts':{
        'static_a'  :11.7                                       # Silvaco
      # 'static_a'  :11.9                                       # K.W. Boer, Survey of Semiconductor Physics, Vol. 2 (1990)
        ,'optical_a' :11.85                                      # Landolt-Boernstein
      # ,'optical_a' :13                                         # Dielectric constant at 10GHz  www.crystran.co.uk/sidata.htm
    }

    ,'elastic_consts':{
        'c11' :165.77   ,'c12' :63.93   ,'c44' :79.62               # 298 K  Landolt-Boernstein
      # 'c11' :167.7    ,'c12' :64.98   ,'c44' :80.35               #  77 K  A. Dargys, J. Kundrotas: Handbook on Physical Properties of Ge,Si,GaAs and InP
      # 'c11' :167      ,'c12' :65      ,'c44' :80                  #        www.crystran.co.uk/sidata.htm
    }
                                             
    ,'piezoelectric_consts':{
        'e14' :0                                                # Piezoelectricity only occurs in III-V materials but not in group IV.
    }                                                
                                                                      
    ,'conduction_bands':{
        'Gamma':{ 
          'mass'                  :0.156
          ,'bandgap'               :3.385                        # 0 K
          ,'bandgap_alpha'         :0.4730e-3                    # This is ,'a' guess! (DELTA valley value was taken.)
          ,'bandgap_beta'          :636                          # This is ,'a' guess! (DELTA valley value was taken.)
          ,'defpot_absolute'       :-10.39                       # A. Zunger: a_c :a_v + a_gap :2.05 - 12.44 :-10.39
        }
        ,'L':{ 
          'mass_l'                :1.420
          ,'mass_t'                :0.130
          ,'bandgap'               :2.01                         # 0 K [J. Weber et al., PRB 40, 5683 (1989)]
          ,'bandgap_alpha'         :0.4730e-3                    # This is ,'a' guess! (DELTA valley value was taken.)
          ,'bandgap_beta'          :636                          # This is ,'a' guess! (DELTA valley value was taken.)
          ,'defpot_absolute'       :-2.02                        # A. Zunger: a_c :a_v + a_gap :2.05 - 4.07 :-2.02
        # ,'defpot_absolute'       :-0.66                        # C. van de Walle et al., PRB 34, 5621 (1986)
          ,'defpot_uniaxial'       :16.14                        # C. van de Walle et al., PRB 34, 5621 (1986) - Xi_u(at minimum), theoretical value
        }
        ,'Delta':{                          
          'mass_l'                :0.916                        # K.W. Boer, Survey of Semiconductor Physics, Vol. 2 (1990)
          ,'mass_t'                :0.190                        # K.W. Boer, Survey of Semiconductor Physics, Vol. 2 (1990)
          ,'bandgap'               :1.17                         #   0 K: 1.17 eV  www.ioffe.rssi.ru/SVA/NSM/Semicond/Si/bandstr.html
        # ,'bandgap'               :1.12                         # 300 K: 1.12 eV  www.ioffe.rssi.ru/SVA/NSM/Semicond/Si/bandstr.html
        # ,'bandgap_alpha'         :0.5367e-3                    # R. Hull: Prop. of Cryst. Si
        # ,'bandgap_beta'          :745.8                        # R. Hull: Prop. of Cryst. Si
          ,'bandgap_alpha'         :0.4730e-3                    # S.M. Sze (1981) and DESSIS
          ,'bandgap_beta'          :636                          # S.M. Sze (1981) and DESSIS
          ,'defpot_absolute'       :3.40                         # A. Zunger: a_c :a_v + a_gap :2.05 + 1.35 :3.40
        # ,'defpot_absolute'       :3.3                          # experimental value of Cargill et al., PRL 61, 1748 (1988)
        # ,'defpot_absolute'       :4.18                         # C. van de Walle et al., PRB 34, 5621 (1986)
          ,'defpot_uniaxial'       :9.16                         # C. van de Walle et al., PRB 34, 5621 (1986) - Xi_u(at minimum), theoretical value
        # ,'defpot_uniaxial'       :8.6                          # 8.6+-0.4 is experimental value, Laude et al., PRB 3, 2623 (1971)
          ,'position'              :0.85                         # 0.85 for DELTA instead of 1.0 for X valley
          ,'g_l'                   :2.00232                      # C. Tahan et al., PRB 71, 075315 (2005)
          ,'g_t'                   :2.00232                      # C. Tahan et al., PRB 71, 075315 (2005)
        }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.090                              # take Qteish value of -6.93 and shift it by 8.02 to align it with Zunger's average ,'valence' band energy (van der Walle model)
        
        ,'HH':{ 'mass'          :0.537 }                            # K.W. Boer, Survey of Semiconductor Physics, Vol. 2 (1990)
        ,'LH':{ 'mass'          :0.153 }                            # K.W. Boer, Survey of Semiconductor Physics, Vol. 2 (1990)
        ,'SO':{ 'mass'          :0.234 }                            # K.W. Boer, Survey of Semiconductor Physics, Vol. 2 (1990) 
        
        ,'defpot_absolute'   : 2.05                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 2.46                              # C. van de Walle, PRB 39, 1871 (1989), theoretical value
      # ,'defpot_absolute'   : 1.80                              # calculated by van de Walle from experimental values of Laude et al. PRB 3, 2623 (1971) and Bardeen et al. PR 80, 72 (1950)
      # ,'defpot_uniaxial_b' :-2.35   ,'defpot_uniaxial_d' :-5.32  # C. van de Walle et al., PRB 34, 5621 (1986), theoretical value
        ,'defpot_uniaxial_b' :-2.10   ,'defpot_uniaxial_d' :-4.85  # Laude et al., PRB 3, 2623 (1971), experimental value (-2.10+-0.10, -4.85+-0.15)

        ,'delta_SO'          :0.044
    }

    ,'kp_6_bands':{
      # gamma1 :   gamma2 :   gamma3 :
        'L' :-6.69   ,'M' :-4.62   ,'N' :-8.56                      # M. Rieger, P. Vogl, PRB 48, 14276 (1993) - theoretical value
      # 'L' :-6.64   ,'M' :-4.60   ,'N' :-8.68                      # M. Rieger, P. Vogl, PRB 48, 14276 (1993) - calculated from experimental Luttinger parameters of O. Madelung (Landolt-Boernstein)
      #  ,'kappa' :-0.26                                          # P. Lawaetz, PRB 4, 3460 (1971)
        ,'kappa' :-0.42                                          # Landolt-Börnstein
    }
 
    ,'kp_8_bands':{                                                # No useful model for indirect ,'bandgap' materials!
        'S' :6.41                                               # inverse effective mass
        ,'E_P' :0                                                # decoupled bands
        ,'B' :0                                                  # This value is zero in group IV semiconductors (inversion symmetry).
        ,'L' :-6.69   ,'M' :-4.62   ,'N' :-8.56                      # 6-band parameters
        ,'kappa' :-0.42                                          # 6-band parameter
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax'  :1417.0   ,'exponent' :2.5  }                 # DESSIS
      # 'electrons':{  'mumax'  :1430     ,'exponent' :2    }                 # PhD thesis V. Palankovski but opposite sign for ,'exponent'
        ,'holes':{      'mumax'  :470.5    ,'exponent' :2.2  }                 # DESSIS
      # ,'holes':{      'mumax'  :460      ,'exponent' :2.18 }                 # PhD thesis V. Palankovski but opposite sign for ,'exponent'
    }
    
    ,'mobility_masetti':{
        'electrons':{  'mumax'  :1417.0   ,'exponent' :2.5                    # DESSIS (same as ,'mobility_constant':{})
                    ,'mumin1' :52.2     ,'mumin2'   :52.2     ,'mu1' :43.4    # DESSIS
                    ,'pc'     :0        ,'cr'       :9.68e16  ,'cs' :3.34e20  # DESSIS
                    ,'alpha'  :0.680    ,'beta'     :2.0                    # DESSIS
        }
        ,'holes':{      'mumax'  :470.5    ,'exponent' :2.2                    # DESSIS (same as ,'mobility_constant':{})
                    ,'mumin1' :44.9     ,'mumin2'   :0        ,'mu1' :29.0    # DESSIS
                    ,'pc'     :9.23e16  ,'cr'       :2.23e17  ,'cs' :6.10e20  # DESSIS
                    ,'alpha'  :0.719    ,'beta'     :2.0                    # DESSIS
        }
    }

    ,'mobility_arora':{
        'electrons':{  'mumin' :88        ,'alm' :-0.57                       # DESSIS
                    ,'mud'   :1252      ,'ald' :-2.33                       # DESSIS
                    ,'n0'    :1.25e17   ,'aln' : 2.4                        # DESSIS
                    ,'a'     :0.88      ,'ala' :-0.146                      # DESSIS
        }
        ,'holes':{      'mumin' :54.3      ,'alm' :-0.57                       # DESSIS
                    ,'mud'   :407       ,'ald' :-2.23                       # DESSIS
                    ,'n0'    :2.35e17   ,'aln' : 2.4                        # DESSIS
                    ,'a'     :0.88      ,'ala' :-0.146                      # DESSIS
        }
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :1430      ,'muLexpT'       :-2           # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')
                    ,'muLImin300' :80        ,'TSwitch'    :200             # PhD thesis V. Palankovski
                    ,'muLIexpTabove' :-0.45  ,'muLIexpTbelow' :-0.15        # PhD thesis V. Palankovski
                    ,'Cref300'    :1.12e17   ,'CrefexpT'      :3.2          # PhD thesis V. Palankovski
                    ,'alpha300'   :0.72      ,'alphaexpT'     :0.065        # PhD thesis V. Palankovski
        }
        ,'holes':{      'muL300'     :460       ,'muLexpT'       :-2.18        # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')
                    ,'muLImin300' :45        ,'TSwitch'    :200             # PhD thesis V. Palankovski
                    ,'muLIexpTabove' :-0.45  ,'muLIexpTbelow' :-0.15        # PhD thesis V. Palankovski
                    ,'Cref300'    :2.23e17   ,'CrefexpT'      :3.2          # PhD thesis V. Palankovski
                    ,'alpha300'   :0.72      ,'alphaexpT'     :0.065        # PhD thesis V. Palankovski
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :4.26e-4     ,'nref_n' :7.1e15                  # SIMBA
                   ,'tau_p' :3.95e-4     ,'nref_p' :7.1e15                  # SIMBA
        }
        ,'Auger':{     'c_n' :2.8e-31       ,'c_p' :9.9e-31                    # SIMBA
        }
    }
}  
#
   
######### germanium ###################################################
,'Ge':{                                                   
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'IV_IV'   

    ,'lattice_consts':{
        'a'           :5.6579                                   # Landolt-Boernstein
        ,'a_expansion' :5.8e-5                                   # [1/K?] ? S.M. Sze, Physics of Semiconductor Devices (1981)
    }
      
    ,'dielectric_consts':{
        'static_a'  :16.2                                       # Landolt-Boernstein
      # 'static_a'  :16.6                                       # at 9.37GHz at 300K   www.crystran.co.uk/germdata.htm
        ,'optical_a' :10.10                                      # reference? GaAs value?
    }

    ,'elastic_consts':{ 
        'c11' :128.53   ,'c12' :48.26   ,'c44' :66.80               # Landolt-Boernstein
      # 'c11' :129      ,'c12' :48.3    ,'c44' :67.1                # www.crystran.co.uk/germdata.htm
    }     

    ,'piezoelectric_consts':{
        'e14' :0                                                # Piezoelectricity only occurs in III-V materials but not in group IV.
    }

    ,'conduction_bands':{
        'Gamma':{ 
          'mass'                  :0.038
        # ,'bandgap'               :0.80                         # 300 K (F. Schaeffler, Semicond. Sci. Technol. 12, 1515 (1997))
          ,'bandgap'               :0.9                          #   0 K (guess for 0 K, any reference available?)
          ,'bandgap_alpha'         :0.4774e-3                    # This is ,'a' guess! (L valley value was taken.)
          ,'bandgap_beta'          :235                          # This is ,'a' guess! (L valley value was taken.)
          ,'defpot_absolute'       :-10.41                       # A. Zunger: a_c :a_v + a_gap :-0.35 - 10.06 :-10.41
        }
        ,'L':{ 
          'mass_l'                :1.57
          ,'mass_t'                :0.0807
          ,'bandgap'               :0.74                         # 0 K [4.2 K, F. Schaeffler, Semicond. Sci. Technol. 12 (1997)]
        # ,'bandgap'               :0.664                        # 300 K, H. Grahn, Semiconductor Physics
          ,'bandgap_alpha'         :0.4774e-3                    # S.M. Sze, Physics of Semiconductor Devices (1981)
          ,'bandgap_beta'          :235                          # S.M. Sze, Physics of Semiconductor Devices (1981)
          ,'defpot_absolute'       :-4.35                        # A. Zunger: a_c :a_v + a_gap :-0.35 - 4.00 :-4.35
        # ,'defpot_absolute'       :-1.54                        # C. van de Walle et al., PRB 34, 5621 (1986)
          ,'defpot_uniaxial'       :15.13                        # C. van de Walle et al., PRB 34, 5621 (1986)
        # ,'defpot_uniaxial'       :10.4                         # M. Fischetti
        }
        ,'Delta':{
          'mass_l'                :1.350
          ,'mass_t'                :0.290 
        # ,'bandgap'               :1.094                        #   0   K
          ,'bandgap'               :0.931                        #   4.2 K, J. Weber et al., PRB 40, 5683 (1989) (X or DELTA?)
        # ,'bandgap'               :1.254                        # 300   K
          ,'bandgap_alpha'         :0.4774e-3                    # This is ,'a' guess! (L valley value was taken.)
          ,'bandgap_beta'          :235                          # This is ,'a' guess! (L valley value was taken.)
          ,'defpot_absolute'       :0.14                         # A. Zunger: a_c :a_v + a_gap :-0.35 + 0.49 :0.14
        # ,'defpot_absolute'       :2.55                         # C. van de Walle et al., PRB 34, 5621 (1986)
          ,'defpot_uniaxial'       :9.42                         # C. van de Walle et al., PRB 34, 5621 (1986)
        # ,'defpot_uniaxial'       :9.75                         # M. Fischetti
          ,'position'              :0.85                         # 0.85 for DELTA instead of 1.0 for X valley
          ,'g_l'                   :0.82                         # F.A. Baron et al., PRB 68, 195306 (2003)
          ,'g_t'                   :1.93                         # F.A. Baron et al., PRB 68, 195306 (2003)
        }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.67                               # see comments at beginning of Si database entry
      # 'bandoffset'        :1.830                              # take Qteish value of -6.19 and shift it by 8.02 to align it with A. Zunger's average ,'valence' band energy (van der Walle model)
                                   
      # ,'HH':{ 'mass'          :0.316  }                           # (Reference?)
      # ,'LH':{ 'mass'          :0.0424 }                           # (Reference?)
      # ,'SO':{ 'mass'          :0.095  }                           # (Reference?)
        ,'HH':{ 'mass'          :0.33   }                           # http://www.ioffe.ru/SVA/NSM/Semicond/Ge/bandstr.html
        ,'LH':{ 'mass'          :0.043  }                           # http://www.ioffe.ru/SVA/NSM/Semicond/Ge/bandstr.html
        ,'SO':{ 'mass'          :0.084  }                           # http://www.ioffe.ru/SVA/NSM/Semicond/Ge/bandstr.html
                                   
        ,'defpot_absolute'   :-0.35                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 1.24                              # C. van de Walle, PRB 39, 1871 (1989), theoretical value
      # ,'defpot_uniaxial_b' :-2.55   ,'defpot_uniaxial_d' :-5.5   # C. van de Walle et al., PRB 34, 5621 (1986), theoretical value
        ,'defpot_uniaxial_b' :-2.86   ,'defpot_uniaxial_d' :-5.28  # M. Chandrasekhar et al., PRB 15, 2127 (1977), experimental value (-2.86+-0.15,-5.28+-0.50)

        ,'delta_SO'          :0.30                               # M. Cardona et al. in Landolt-Boernstein

    }

    ,'kp_6_bands':{
      # gamma1 :13.38   gamma2 :4.24   gamma3 :5.69         # J.C. Hensel, K. Suzuki PRB 9, 4219 (1974)
        'L' :-31.34       ,'M' :-5.90       ,'N' :-34.14            # M. Rieger Diploma thesis - experimental value, calculated from experimental Luttinger parameters of J.C. Hensel, K. Suzuki PRB 9, 4219 (1974)
      # 'L' :-21.65       ,'M' :-5.02       ,'N' :-23.48            # M. Rieger, P. Vogl, PRB 48, 14276 (1993) - theoretical value
        ,'kappa' :3.41                                           # P. Lawaetz, PRB 4, 3460 (1971) + Landolt-Börnstein
    }                                   
                                       
    ,'kp_8_bands':{                                                # No useful model for indirect ,'bandgap' materials!
        'S' :26.32                                              # inverse effective mass
        ,'E_P' :0                                                # decoupled bands
        ,'B' :0                                                  # This value is zero in group IV semiconductors (inversion symmetry).
        ,'L' :-31.34       ,'M' :-5.90       ,'N' :-34.14            # 6-band parameters
        ,'kappa' :3.41                                           # 6-band parameter
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :3800     ,'exponent' :1.66 }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
        ,'holes':{      'mumax' :1800     ,'exponent' :2.33 }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :3800      ,'muLexpT'       :-1.66        # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')
                    ,'muLImin300' :850       ,'TSwitch'    :200             # PhD thesis V. Palankovski
                    ,'muLIexpTabove' :0      ,'muLIexpTbelow' :0            # PhD thesis V. Palankovski
                    ,'Cref300'    :2.6e17    ,'CrefexpT'      :0            # PhD thesis V. Palankovski
                    ,'alpha300'   :0.56      ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
        ,'holes':{      'muL300'     :1800      ,'muLexpT'       :-2.33        # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')
                    ,'muLImin300' :300       ,'TSwitch'    :200             # PhD thesis V. Palankovski
                    ,'muLIexpTabove' :0      ,'muLIexpTbelow' :0            # PhD thesis V. Palankovski
                    ,'Cref300'    :1.0e17    ,'CrefexpT'      :0            # PhD thesis V. Palankovski
                    ,'alpha300'   :1.0       ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :4.26e-4     ,'nref_n' :7.1e15                  # SIMBA
                   ,'tau_p' :3.95e-4     ,'nref_p' :7.1e15                  # SIMBA
        }
        ,'Auger':{     'c_n' :1.0e-31       ,'c_p' :1.0e-31                    # SIMBA
        }
    }
}
#









##########################################################################################
##########################################################################################
#          B I N A R I E S     --    III - V       V A L E N C E
##########################################################################################
##########################################################################################




######### gallium arsenide ############################################
,'GaAs':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :5.65325                                  # Vurgaftman1 (300 K)
        ,'a_expansion' :3.88e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :12.93
        ,'optical_a' :10.10
    }

    ,'elastic_consts':{
        'c11' :122.1   ,'c12' :56.6   ,'c44' :60.0                  # Vurgaftman1
    }                    

    ,'piezoelectric_consts':{
     #  'e14' :-0.175                                           # calculated by      S. Gironcoli et al., PRL 62(24), 2853 (1989)
        'e14' :-0.160                                           # experimental value S. Gironcoli et al., PRL 62(24), 2853 (1989)
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.067                        # Vurgaftman1
          ,'bandgap'               :1.519                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.5405e-3                    # Vurgaftman1
          ,'bandgap_beta'          :204                          # Vurgaftman1
          ,'defpot_absolute'       :-9.36                        # A. Zunger: a_c :a_v + a_gap :-1.21 - 8.15 :-9.36
        # ,'defpot_absolute'       :-7.17                        # Vurgaftman1
        # ,'g'                     :-0.44                        #   4 K, M. Oestreich et al., PRB 53, 7911 (1996)
          ,'g'                     :-0.30                        # 280 K, M. Oestreich et al., PRB 53, 7911 (1996)
       }
       ,'L':{ 
          'mass_l'                :1.9                          # Vurgaftman1
          ,'mass_t'                :0.0754                       # Vurgaftman1
          ,'bandgap'               :1.815                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.605e-3                     # Vurgaftman1
          ,'bandgap_beta'          :204                          # Vurgaftman1
          ,'defpot_absolute'       :-4.91                        # A. Zunger: a_c :a_v + a_gap :-1.21 - 3.70 :-4.91
          ,'defpot_uniaxial'       :14.26                        # C. van de Walle
       }
       ,'X':{ 
          'mass_l'                :1.3                          # Vurgaftman1
          ,'mass_t'                :0.23                         # Vurgaftman1
          ,'bandgap'               :1.981                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.460e-3                     # Vurgaftman1
          ,'bandgap_beta'          :204                          # Vurgaftman1
          ,'defpot_absolute'       :-0.16                        # A. Zunger: a_c :a_v + a_gap :-1.21 + 1.05 :-0.16
          ,'defpot_uniaxial'       :8.61                         # C. van de Walle
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.346                              # A. Zunger
        
        ,'HH':{ 'mass'          :0.51  ,'g' :-7.86 }                  # m_hh :http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/bandstr.html
        ,'LH':{ 'mass'          :0.082 ,'g' :-2.62 }                  # m_lh :http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/bandstr.html
#       ,'SO':{ 'mass'          :0.15  }                            # m_so :http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/bandstr.html
        ,'SO':{ 'mass'          :0.172 }                            # Vurgaftman1

        ,'defpot_absolute'   :-1.21                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 1.16                              # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -1.16
        ,'defpot_uniaxial_b' :-2.0   ,'defpot_uniaxial_d' :-4.8
       
        ,'delta_SO'          :0.341                              # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :6.98   gamma2 :2.06   gamma3 :2.93          # Vurgaftman1
        'L' :-16.22      ,'M' : -3.86      ,'N' :-17.58                
        ,'kappa' :1.2                                            # Kiselev, PRB 64, 125303 (2001)
      # ,'kappa' :1.72                                           # P. Lawaetz, PRB 4, 3460 (1971)
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :1.519
        'S' :-2.88                                              # S :1 + 2F :1 + 2 (-1.94) :-2.88 (Vurgaftman1)
        ,'E_P' :28.8                                             # Vurgaftman1
        ,'B' : 0                     
        ,'L' : 1.420   ,'M' :-3.86    ,'N' : 0.056
        ,'kappa' :-1.74                                          # Kiselev
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :8500     ,'exponent' :2.2  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
        ,'holes':{      'mumax' :800      ,'exponent' :0.9  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :8500      ,'muLexpT'       :-2.2         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :800       ,'muLIexpTabove' :-0.9         # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :-0.9         # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :1.0e17    ,'CrefexpT'      :6.2          # PhD thesis V. Palankovski
                    ,'alpha300'   :0.5       ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
        ,'holes':{      'muL300'     :800       ,'muLexpT'       :-0.9         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :40        ,'muLIexpTabove' :0            # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :0            # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :1.0e17    ,'CrefexpT'      :0.5          # PhD thesis V. Palankovski
                    ,'alpha300'   :1.0       ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19                  # SIMBA
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18                  # SIMBA
        }
      # ,'Auger':{     'c_n' :1.0e-31       ,'c_p' :1.0e-31                    # SIMBA
        ,'Auger':{     'c_n' :1.0e-30       ,'c_p' :1.0e-30                    # 300 K, http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/electric.html#Recombination
        }
      # ,'radiative':{ 'c' :2.0e-10                                          # DESSIS
        ,'radiative':{ 'c' :7.2e-10                                          # Ioffe, 300 K, V. P. Varshni, Phys. Status Solidi 19, 459 (1967); 20, 9 (1967), http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/electric.html#Recombination
        }
    }
}
#
   
######### aluminum arsenide ###########################################
,'AlAs':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :5.6611                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :2.90e-5                                  # Vurgaftman1
    }
   
    ,'dielectric_consts':{
        'static_a'  :10.064                                     # Landolt-Boernstein
        ,'optical_a' :8.162                                      # Landolt-Boernstein
    }

    ,'elastic_consts':{
        'c11' :125.0   ,'c12' :53.4   ,'c44' :54.2                  # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
      # 'e14' :-0.22                                            # guess from ,'a' picture http://nina.ecse.rpi.edu/shur/Tutorial/GaNtutorial1/sld036.htm
        'e14' :-0.015                                           # calculated by S. Gironcoli et al., PRL 62(24), 2853 (1989)
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.15                         # Vurgaftman1
          ,'bandgap'               :3.099                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.885e-3                     # Vurgaftman1
          ,'bandgap_beta'          :530                          # Vurgaftman1
          ,'defpot_absolute'       :-7.40                        # A. Zunger: a_c :a_v + a_gap :1.53 - 8.93 :-7.40
        # ,'defpot_absolute'       :-5.64                        # Vurgaftman1
          ,'g'                     :1.52                         # J.-M. Jancu, R. Scholz, PRB 72, 193201 (2005)
       }
       ,'L':{ 
          'mass_l'                :1.32                         # Vurgaftman1
          ,'mass_t'                :0.15                         # Vurgaftman1
          ,'bandgap'               :2.46                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.605e-3                     # Vurgaftman1
          ,'bandgap_beta'          :204                          # Vurgaftman1
          ,'defpot_absolute'       :-3.07                        # A. Zunger: a_c :a_v + a_gap :1.53 - 4.60 :-3.07
          ,'defpot_uniaxial'       :11.35                        # InAs value !!!
       }
       ,'X':{ 
          'mass_l'                :0.97                         # Vurgaftman1
          ,'mass_t'                :0.22                         # Vurgaftman1
          ,'bandgap'               :2.24                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.70e-3                      # Vurgaftman1
          ,'bandgap_beta'          :530                          # Vurgaftman1
          ,'defpot_absolute'       :2.54                         # A. Zunger: a_c :a_v + a_gap :1.53 + 1.01 :2.54
          ,'defpot_uniaxial'       :6.11                         # Munoz
          ,'g_l'                   :1.9                          # J.D. Caldwell et al., PRB 72, 115339 (2005)
          ,'g_t'                   :1.9                          # J.D. Caldwell et al., PRB 72, 115339 (2005)
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.857                              # A. Zunger
        
        ,'HH':{ 'mass'          :0.5   }                            # Landolt-Boernstein
        ,'LH':{ 'mass'          :0.26  }                            # Landolt-Boernstein
        ,'SO':{ 'mass'          :0.28  }                            # Vurgaftman1
        
        ,'defpot_absolute'   :1.53                               # A. Zunger: a_v
      # ,'defpot_absolute'   :2.47                               # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -2.47
        ,'defpot_uniaxial_b' :-2.3   ,'defpot_uniaxial_d' :-3.4    # Vurgaftman1
       
        ,'delta_SO'          :0.28                               # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :3.76   gamma2 :0.82   gamma3 :1.42          # Vurgaftman1
        'L' :-8.04       ,'M' :-3.12       ,'N' :-8.52               
        ,'kappa' :0.12                                           # P. Lawaetz, PRB 4, 3460 (1971) 
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :3.099
        'S' : 0.04                                              # S :1 + 2F :1 + 2 (-0.48) :0.04 (Vurgaftman1)
        ,'E_P' :21.1                                             # Vurgaftman1
        ,'B' : 0                       
        ,'L' :-1.430   ,'M' :-3.12   ,'N' :-1.910
        ,'kappa' :-0.982                
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :410      ,'exponent' :2.1  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
        ,'holes':{      'mumax' :130      ,'exponent' :2.2  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :410       ,'muLexpT'       :-2.1         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :10        ,'muLIexpTabove' :0            # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :0            # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :1.0e17    ,'CrefexpT'      :0            # PhD thesis V. Palankovski
                    ,'alpha300'   :0.5       ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
        ,'holes':{      'muL300'     :130       ,'muLexpT'       :-2.2         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :5         ,'muLIexpTabove' :0            # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :0            # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :2.9e17    ,'CrefexpT'      :0.5          # PhD thesis V. Palankovski
                    ,'alpha300'   :1.0       ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19                  # SIMBA
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18                  # SIMBA
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                          # SIMBA
        }
        ,'radiative':{ 'c' :0                                                # ?
        }
    }
}
#

######### indium arsenide #############################################
,'InAs':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :6.0583                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :2.74e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :15.15                                      # Landolt-Boernstein epsilon(0)
        ,'optical_a' :12.25                                      # Landolt-Boernstein epsilon(infinity)
    }

    ,'elastic_consts':{
        'c11' :83.29   ,'c12' :45.26   ,'c44' :39.59                # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
      # 'e14' :-0.035                                           # calculated by      S. Gironcoli et al., PRL 62(24), 2853 (1989)
        'e14' :-0.044                                           # experimental value S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.0459                                          # Landolt-Boernstein
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.026                        # Vurgaftman1
          ,'bandgap'               :0.417                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.276e-3                     # Vurgaftman1
          ,'bandgap_beta'          :93                           # Vurgaftman1
        # ,'defpot_absolute'       :-5.08                        # Vurgaftman1
          ,'defpot_absolute'       :-6.66                        # A. Zunger: a_c :a_v + a_gap :-1.00 - 5.66 :-6.66
        # ,'g'                     :-15.6                        # PhD thesis S. Hackenbuchner, p. 104
          ,'g'                     :-14.9                        # J.-M. Jancu, R. Scholz, PRB 72, 193201 (2005)
       }
       ,'L':{ 
          'mass_l'                :0.64                         # Vurgaftman1
          ,'mass_t'                :0.05                         # Vurgaftman1
          ,'bandgap'               :1.133                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.276e-3                     # Vurgaftman1
          ,'bandgap_beta'          :93                           # Vurgaftman1
          ,'defpot_absolute'       :-3.89                        # A. Zunger: a_c :a_v + a_gap :-1.00 - 2.89 :-3.89
          ,'defpot_uniaxial'       :11.35                        # C. van de Walle
       }
       ,'X':{ 
          'mass_l'                :1.13                         # Vurgaftman1
          ,'mass_t'                :0.16                         # Vurgaftman1
          ,'bandgap'               :1.433                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.276e-3                     # Vurgaftman1
          ,'bandgap_beta'          :93                           # Vurgaftman1
          ,'defpot_absolute'       :-0.08                        # A. Zunger: a_c :a_v + a_gap :-1.00 + 0.92 :-0.08
          ,'defpot_uniaxial'       :3.7                          # Munoz
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.390                              # A. Zunger
        
        ,'HH':{ 'mass'          :0.41  ,'g' :-45.2 }                  # http://www.ioffe.ru/SVA/NSM/Semicond/InAs/bandstr.html
        ,'LH':{ 'mass'          :0.026 ,'g' :-15.1 }                  # Yu, Cardona, Fundamentals of Semiconductors, p. 70, http://www.ioffe.ru/SVA/NSM/Semicond/InAs/bandstr.html
      # ,'LH':{ 'mass'          :0.025 }                            # Landolt-Boernstein
      # ,'SO':{ 'mass'          :0.016 }                            # http://www.ioffe.ru/SVA/NSM/Semicond/InAs/bandstr.html
        ,'SO':{ 'mass'          :0.014 }                            # Vurgaftman1
        
        ,'defpot_absolute'   :-1.00                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 1.00                              # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -1.00
        ,'defpot_uniaxial_b' :-1.8   ,'defpot_uniaxial_d' :-3.6    # Vurgaftman1 - N.E. Christensen et al., PRB 36 (2), 1032 (1987) suggest to revise experimental value of d=-3.6.
       
        ,'delta_SO'          :0.39                               # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :20.0   gamma2 :8.5   gamma3 :9.2            # Vurgaftman1
        'L' :-55.0       ,'M' :-4.0       ,'N' :-55.2        
        ,'kappa' :7.68                                            # Kiselev, PRB 64, 125303 (2001) + P. Lawaetz, PRB 4, 3460 (1971)       
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :0.417
        'S' : -4.80                                             # S :1 + 2F :1 + 2 (-2.90) :-4.80 (Vurgaftman1)
        ,'E_P' :21.5                                             # Vurgaftman1
        ,'B' :  0                   
        ,'L' :-15.695   ,'M' :-4.0   ,'N' :-15.895
        ,'kappa' :1.129               
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :32500    ,'exponent' :1.7  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
        ,'holes':{      'mumax' :510      ,'exponent' :2.3  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :32500     ,'muLexpT'       :-1.7         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :11700     ,'muLIexpTabove' :-0.33        # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :-0.33        # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :4.4e16    ,'CrefexpT'      : 3.6         # PhD thesis V. Palankovski
                    ,'alpha300'   :0.5       ,'alphaexpT'     : 0           # PhD thesis V. Palankovski
        }
        ,'holes':{      'muL300'     :510       ,'muLexpT'       :-2.3         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :48        ,'muLIexpTabove' :0            # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :0            # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :2.55e17   ,'CrefexpT'      :0.5          # PhD thesis V. Palankovski
                    ,'alpha300'   :1.0       ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19                  # SIMBA
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18                  # SIMBA
        }
        ,'Auger':{     'c_n' :1.0e-31       ,'c_p' :1.0e-31                    # SIMBA
        }
        ,'radiative':{ 'c' :0                                                # ?
        }
    }
}
#

######### gallium phosphide ###########################################
,'GaP':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :5.4505                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :2.92e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :11.1                                       # Landolt-Boernstein epsilon(0)
        ,'optical_a' : 9.075                                     # Landolt-Boernstein epsilon(infinity)
    }

    ,'elastic_consts':{
        'c11' :140.5   ,'c12' :62.03   ,'c44' :70.33                # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
      # 'e14' :-0.097                                           # calculated by      S. Gironcoli et al., PRL 62(24), 2853 (1989)
        'e14' :-0.097                                           # experimental value S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.1                                             # Landolt-Boernstein - was already negative in Landolt-Boernstein
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.13                         # Vurgaftman1
          ,'bandgap'               :2.886                        # Vurgaftman1 (0 K) 2.886 + 0.1081 [1 - coth(164/T)]
          ,'bandgap_alpha'         :0.5771e-3                    # L and X valley ??? => see formula above
          ,'bandgap_beta'          :372                          # L and X valley ??? => see formula above
        # ,'defpot_absolute'       :-8.2                         # Vurgaftman1
          ,'defpot_absolute'       :-9.41                        # A. Zunger: a_c :a_v + a_gap :-0.58 - 8.83 :-9.41
          ,'g'                     :1.9                          # J.-M. Jancu, R. Scholz, PRB 72, 193201 (2005)
       }
       ,'L':{ 
          'mass_l'                :1.2                          # Vurgaftman1
          ,'mass_t'                :0.15                         # Vurgaftman1
          ,'bandgap'               :2.72                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.5771e-3                    # Vurgaftman1
          ,'bandgap_beta'          :372                          # Vurgaftman1
          ,'defpot_absolute'       :-4.41                        # A. Zunger: a_c :a_v + a_gap :-0.58 - 3.83 :-4.41
          ,'defpot_uniaxial'       :11.35                        # InAs value
       }
     # ,'Delta':{                                                  # GaP has DELTA vallye (and not X valley)
       ,'X':{
          'mass_l'                :2.0                          # Vurgaftman1
          ,'mass_t'                :0.253                        # Vurgaftman1
          ,'bandgap'               :2.35                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.5771e-3                    # Vurgaftman1
          ,'bandgap_beta'          :372                          # Vurgaftman1
          ,'defpot_absolute'       :0.69                         # A. Zunger: a_c :a_v + a_gap :-0.58 + 1.27 :0.69
          ,'defpot_uniaxial'       :5.65                         # Munoz
        # ,'position'              :0.95                         # 0.95 for DELTA instead of 1.0 for X valley
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.963                              # A. Zunger
        
      # ,'HH':{ 'mass'          :0.57  }                            # Yu, Cardona, Fundamentals of Semiconductors, p. 70
        ,'HH':{ 'mass'          :0.79  }                            # Landolt-Boernstein (m_p,h), http://www.ioffe.ru/SVA/NSM/Semicond/GaP/bandstr.html
        ,'LH':{ 'mass'          :0.14  }                            # Landolt-Boernstein (m_p,l), http://www.ioffe.ru/SVA/NSM/Semicond/GaP/bandstr.html
        ,'SO':{ 'mass'          :0.25  }                            # Vurgaftman1
        
        ,'defpot_absolute'   :-0.58                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 1.7                               # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -1.7
        ,'defpot_uniaxial_b' :-1.6   ,'defpot_uniaxial_d' :-4.6    # Vurgaftman1
       
        ,'delta_SO'          :0.08                               # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.05   gamma2 :0.49   gamma3 :1.25          # (ERROR in Vurgaftman1 paper!!!)
        'L' :-7.01       ,'M' :-4.07       ,'N' :-7.50   
        ,'kappa' :0.34                                           # Kiselev, PRB 64, 125303 (2001)+ P. Lawaetz, PRB 4, 3460 (1971)            
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :2.886
        'S' :-3.08                                              # S :1 + 2F :1 + 2 (-2.04) :-3.08 (Vurgaftman1)
        ,'E_P' :31.4                                             # Vurgaftman1
        ,'B' : 0                   
        ,'L' : 5.7705   ,'M' :-4.07   ,'N' :3.280  
        ,'kappa' :-1.457             
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :210      ,'exponent' :1.94 }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
        ,'holes':{      'mumax' :160      ,'exponent' :2.0  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :210       ,'muLexpT'       :-1.94        # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :76        ,'muLIexpTabove' :-1.07        # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :-1.07        # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :2.85e17   ,'CrefexpT'      : 1.8         # PhD thesis V. Palankovski
                    ,'alpha300'   :0.5       ,'alphaexpT'     : 0           # PhD thesis V. Palankovski
        }
        ,'holes':{      'muL300'     :160       ,'muLexpT'       :-2.0         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :27        ,'muLIexpTabove' :0            # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' :0            # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :2.33e17   ,'CrefexpT'      :0            # PhD thesis V. Palankovski
                    ,'alpha300'   :1.0       ,'alphaexpT'     :0            # PhD thesis V. Palankovski
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19                  # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18                  # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                          # InP value !!!
        }
        ,'radiative':{ 'c' :0                                                # ?
        }
    }
}
#

######### aluminium phosphide #########################################
,'AlP':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :5.4672                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :2.92e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :9.8                                        # Landolt-Boernstein epsilon(0)
        ,'optical_a' :7.5385                                     # Landolt-Boernstein epsilon(infinity) => epsilon(0) / epsilon(infinity) :1.3
    }

    ,'elastic_consts':{
        'c11' :133.0   ,'c12' :63.0   ,'c44' :61.5                  # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :0.059                                            # POSITIVE value - calculated by S. Gironcoli et al., PRL 62(24), 2853 (1989)
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.22                         # Vurgaftman1
          ,'bandgap'               :3.63                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.5771e-3                    # Vurgaftman1
          ,'bandgap_beta'          :372                          # Vurgaftman1
        # ,'defpot_absolute'       :-5.7                         # Vurgaftman1
          ,'defpot_absolute'       :-6.88                        # A. Zunger: a_c :a_v + a_gap :2.64 - 9.52 :-6.88
          ,'g'                     :1.92                         # J.-M. Jancu, R. Scholz, PRB 72, 193201 (2005)
       }
       ,'L':{ 
          'mass_l'                :1                            # Vurgaftman1
          ,'mass_t'                :0.1                          # ???
          ,'bandgap'               :3.57                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.318e-3                     # Vurgaftman1
          ,'bandgap_beta'          :588                          # Vurgaftman1
          ,'defpot_absolute'       :-1.74                        # A. Zunger: a_c :a_v + a_gap :2.64 - 4.38 :-1.74
          ,'defpot_uniaxial'       :11.35                        # InAs value
       }
       ,'X':{
          'mass_l'                :2.68                         # Vurgaftman1
          ,'mass_t'                :0.155                        # Vurgaftman1
          ,'bandgap'               :2.52                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.318e-3                     # Vurgaftman1
          ,'bandgap_beta'          :588                          # Vurgaftman1
          ,'defpot_absolute'       :3.98                         # A. Zunger: a_c :a_v + a_gap :2.64 + 1.34 :3.98
          ,'defpot_uniaxial'       :6.75                         # Munoz
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.427                              # A. Zunger
        
        ,'HH':{ 'mass'          :0.63  }                            # Landolt-Boernstein (m_p,h)
        ,'LH':{ 'mass'          :0.20  }                            # Landolt-Boernstein (m_p,l)
        ,'SO':{ 'mass'          :0.30  }                            # Vurgaftman1
        
        ,'defpot_absolute'   : 2.64                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 3.0                               # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -3.0
        ,'defpot_uniaxial_b' :-1.5   ,'defpot_uniaxial_d' :-4.6    # Vurgaftman1
       
        ,'delta_SO'          :0.07                               # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :3.35   gamma2 :0.71   gamma3 :1.23          # Vurgaftman1
        'L' :-7.19       ,'M' :-2.93       ,'N' :-7.38  
        ,'kappa' :-0.54                                          # P. Lawaetz, PRB 4, 3460 (1971)             
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :2.52
        'S' :-0.30                                              # S :1 + 2F :1 + 2 (-0.65) :-0.30 (Vurgaftman1)
        ,'E_P' :17.7                                             # Vurgaftman1
        ,'B' : 0                   
        ,'L' :-2.345   ,'M' :-2.93   ,'N' :-2.535  
        ,'kappa' :-1.347             
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :210      ,'exponent' :1.94 }         # GaP values
        ,'holes':{      'mumax' :160      ,'exponent' :2.0  }         # GaP values
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :0                                       # ?
        }
    }
}
#

######### indium phosphide ############################################
,'InP':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :5.8697                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :2.79e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :12.61                                      # Landolt-Boernstein epsilon(0)
        ,'optical_a' : 9.61                                      # Landolt-Boernstein epsilon(infinity)
    }

    ,'elastic_consts':{
        'c11' :101.1   ,'c12' :56.1   ,'c44' :45.6                  # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' : 0.056                                           # POSITIVE value - calculated by S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.042                                           #             experimental value S. Gironcoli et al., PRL 62(24), 2853 (1989)
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.0795                       # Vurgaftman1
          ,'bandgap'               :1.4236                       # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.363e-3                     # Vurgaftman1
          ,'bandgap_beta'          :162                          # Vurgaftman1
        # ,'defpot_absolute'       :-6.0                         # Vurgaftman1
          ,'defpot_absolute'       :-6.34                        # A. Zunger: a_c :a_v + a_gap :-0.41 - 5.93 :-6.34
          ,'g'                     :1.20                         # 4 K (up to 160 K), M. Oestreich et al., PRB 53, 7911 (1996) (small temperature dependence)
       }
       ,'L':{ 
          'mass_l'                :0.47                         # Vurgaftman1 - m_DOS(L)
          ,'mass_t'                :0.47                         # Vurgaftman1 - m_DOS(L)
          ,'bandgap'               :2.014                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.363e-3                     # Vurgaftman1
          ,'bandgap_beta'          :162                          # Vurgaftman1
          ,'defpot_absolute'       :-3.41                        # A. Zunger: a_c :a_v + a_gap :-0.41 - 3.00 :-3.41
          ,'defpot_uniaxial'       :11.35                        # InAs value
       }
       ,'X':{
          'mass_l'                :0.88                         # Vurgaftman1 - m_DOS(X)
          ,'mass_t'                :0.88                         # Vurgaftman1 - m_DOS(X)
          ,'bandgap'               :2.384                        # Vurgaftman1 (0 K)        ,'bandgap' : 2.384 - 3.7e-4 * T
          ,'bandgap_alpha'         :0.363e-3                     # GAMMA and L values !!!
          ,'bandgap_beta'          :162                          # GAMMA and L values !!!
          ,'defpot_absolute'       :0.59                         # A. Zunger: a_c :a_v + a_gap :-0.41 + 1.00 :0.59
          ,'defpot_uniaxial'       :3.3                          # Munoz
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.064                              # A. Zunger
        
      # ,'HH':{ 'mass'          :0.58  }                            # Yu, Cardona, Fundamentals of Semiconductors, p. 70
      # ,'HH':{ 'mass'          :0.6   }                            # http://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html
        ,'HH':{ 'mass'          :0.85  }                            # Landolt-Boernstein (m_p,h)
        ,'LH':{ 'mass'          :0.089 }                            # Landolt-Boernstein (m_p,l), http://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html
        ,'SO':{ 'mass'          :0.21  }                            # Vurgaftman1
      # ,'SO':{ 'mass'          :0.17  }                            # http://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html

        ,'defpot_absolute'   :-0.41                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 0.6                               # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -0.6
        ,'defpot_uniaxial_b' :-2.0   ,'defpot_uniaxial_d' :-5.0    # Vurgaftman1
       
        ,'delta_SO'          :0.108                              # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :5.08   gamma2 :1.60   gamma3 :2.10          # Vurgaftman1
        'L' :-12.48      ,'M' :-2.88       ,'N' :-12.60    
        ,'kappa' :0.97                                           # Kiselev, PRB 64, 125303 (2001)
      # ,'kappa' :1.47                                           # P. Lawaetz, PRB 4, 3460 (1971)           
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :1.4236
        'S' :-1.62                                              # S :1 + 2F :1 + 2 (-1.31) :-1.62 (Vurgaftman1)
        ,'E_P' :20.7                                             # Vurgaftman1
        ,'B' : 0                   
        ,'L' : 1.7020   ,'M' :-2.88   ,'N' :1.5820 
        ,'kappa' :-1.394                                         # Kiselev              
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :5300     ,'exponent' :1.9  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
        ,'holes':{      'mumax' :200      ,'exponent' :1.2  }         # PhD thesis V. Palankovski but opposite sign for ,'exponent'
    }

    ,'mobility_minimos':{
        'electrons':{  'muL300'     :5300      ,'muLexpT'       :-1.9         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :1520      ,'muLIexpTabove' : 2.0         # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' : 2.0         # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :6.4e16    ,'CrefexpT'      : 3.7         # PhD thesis V. Palankovski
                    ,'alpha300'   :0.5       ,'alphaexpT'     : 0           # PhD thesis V. Palankovski
        }
        ,'holes':{      'muL300'     :200       ,'muLexpT'       :-1.2         # PhD thesis V. Palankovski (same as ,'mobility_constant':{} but opposite sign for ,'exponent')

                    ,'muLImin300' :24        ,'muLIexpTabove' : 1.2         # PhD thesis V. Palankovski
                                           ,'muLIexpTbelow' : 1.2         # PhD thesis V. Palankovski
                    ,'TSwitch'    :200                                    # PhD thesis V. Palankovski
                    ,'Cref300'    :2.5e17    ,'CrefexpT'      : 0.47        # PhD thesis V. Palankovski
                    ,'alpha300'   :1.0       ,'alphaexpT'     : 0           # PhD thesis V. Palankovski
        }
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19                  # SIMBA
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18                  # SIMBA
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                          # SIMBA
        }
        ,'radiative':{ 'c' :0                                                # ?
        }
    }
}
#

######### gallium antimonide ##########################################
,'GaSb':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :6.0959                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :4.72e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :15.69                                      # Landolt-Boernstein epsilon(0)
        ,'optical_a' :14.44                                      # Landolt-Boernstein epsilon(infinity)
    }

    ,'elastic_consts':{
        'c11' :88.42   ,'c12' :40.26   ,'c44' :43.22                 # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :-0.172                                           # calculated by      S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.168                                           # experimental value S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.126                                           # Landolt-Boernstein
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.039                        # Vurgaftman1
          ,'bandgap'               :0.812                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.417e-3                     # Vurgaftman1
          ,'bandgap_beta'          :140                          # Vurgaftman1
        # ,'defpot_absolute'       :-7.5                         # Vurgaftman1
          ,'defpot_absolute'       :-9.33                        # A. Zunger: a_c :a_v + a_gap :-1.32 - 8.01 :-9.33
          ,'g'                     :-9.2                         # J.-M. Jancu, R. Scholz, PRB 72, 193201 (2005)
       }
       ,'L':{ 
          'mass_l'                :1.3                          # Vurgaftman1
          ,'mass_t'                :0.10                         # Vurgaftman1
          ,'bandgap'               :0.875                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.597e-3                     # Vurgaftman1
          ,'bandgap_beta'          :140                          # Vurgaftman1
          ,'defpot_absolute'       :-4.38                        # A. Zunger: a_c :a_v + a_gap :-1.32 - 3.06 :-4.38
          ,'defpot_uniaxial'       :15.0                         # Landolt-Boernstein
       }
       ,'X':{
          'mass_l'                :1.51                         # Vurgaftman1
          ,'mass_t'                :0.22                         # Vurgaftman1
          ,'bandgap'               :1.141                        # Vurgaftman1 (0 K)        ,'bandgap' : 2.384 - 3.7e-4 * T
          ,'bandgap_alpha'         :0.475e-3                     # Vurgaftman1
          ,'bandgap_beta'          :94                           # Vurgaftman1
          ,'defpot_absolute'       :-0.20                        # A. Zunger: a_c :a_v + a_gap :-1.32 + 1.12 :-0.20
          ,'defpot_uniaxial'       :6.46                         # Munoz
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.777                              # A. Zunger
        
      # ,'HH':{ 'mass'          :0.8    }                           # Yu, Cardona, Fundamentals of Semiconductors, p. 70
        ,'HH':{ 'mass'          :0.34   }                           # Landolt-Boernstein (m_p,h)
        ,'LH':{ 'mass'          :0.0447 }                           # Landolt-Boernstein (m_p,l)
        ,'SO':{ 'mass'          :0.12   }                           # Vurgaftman1

      # ,'HH':{ 'mass'          :0.4    }                           # http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/bandstr.html
      # ,'LH':{ 'mass'          :0.05   }                           # http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/bandstr.html
      # ,'SO':{ 'mass'          :0.14   }                           # http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/bandstr.html
       
        ,'defpot_absolute'   :-1.32                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 0.8                               # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -0.8
        ,'defpot_uniaxial_b' :-2.0   ,'defpot_uniaxial_d' :-4.7    # Vurgaftman1
       
        ,'delta_SO'          :0.76                               # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :13.4   gamma2 :4.7   gamma3 :6.0            # Vurgaftman1
        'L' :-33.2       ,'M' :-5.0       ,'N' :-36.0  
        ,'kappa' :3.18                                           # P. Lawaetz, PRB 4, 3460 (1971)             
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :0.812
        'S' :-2.26                                              # S :1 + 2F :1 + 2 (-1.63) :-2.26 (Vurgaftman1)
        ,'E_P' :27.0                                             # Vurgaftman1
        ,'B' : 0                   
        ,'L' :-7.856   ,'M' :-5.0   ,'N' :-10.656  
        ,'kappa' :-1.044             
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :5300     ,'exponent' :1.9  }         # D.Martin & C.Algora - Temperature-depent GaSb material parameters for reliable thermophotovoltaic cell modelling (,'exponent' from InP)
        ,'holes':{      'mumax' :1000     ,'exponent' :1.2  }        # D.Martin & C.Algora - Temperature-depent GaSb material parameters for reliable thermophotovoltaic cell modelling (,'exponent' from InP)

    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # SIMBA
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # SIMBA
        }
        ,'Auger':{     #'c_n' :1.0e-31       ,'c_p' :1.0e-31          # SIMBA
                   'c_n' :5.0e-30       ,'c_p' :5.0e-30           # Stollwerck et. al. - Characterization an Simulation of GaSb Device-Related Properties
        }
        ,'radiative':{ #'c' :8.0e-11                                # Stollwerck et. al. - Characterization an Simulation of GaSb Device-Related Properties
                   'c' :1.0e-9                                  # ?
        }
    }
}
#

######### aluminium antimonide ########################################
,'AlSb':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :6.1355                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :2.60e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :12.04                                      # Landolt-Boernstein epsilon(0)
        ,'optical_a' :10.24                                      # Landolt-Boernstein epsilon(infinity)
    }

    ,'elastic_consts':{
        'c11' :87.69   ,'c12' :43.41   ,'c44' :40.76                # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :-0.055                                           # calculated by      S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.068                                           # experimental value S. Gironcoli et al., PRL 62(24), 2853 (1989) & Landolt-Boernstein
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.14                         # Vurgaftman1
          ,'bandgap'               :2.386                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.42e-3                      # Vurgaftman1
          ,'bandgap_beta'          :140                          # Vurgaftman1
        # ,'defpot_absolute'       :-4.5                         # Vurgaftman1
          ,'defpot_absolute'       :-8.12                        # A. Zunger: a_c :a_v + a_gap :0.73 - 8.85 :-8.12
          ,'g'                     :0.84                         # J.-M. Jancu, R. Scholz, PRB 72, 193201 (2005)
       }
       ,'L':{ 
          'mass_l'                :1.64                         # Vurgaftman1
          ,'mass_t'                :0.23                         # Vurgaftman1
          ,'bandgap'               :2.329                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.58e-3                      # Vurgaftman1
          ,'bandgap_beta'          :140                          # Vurgaftman1
          ,'defpot_absolute'       :-2.91                        # A. Zunger: a_c :a_v + a_gap :0.73 - 3.64 :-2.91
          ,'defpot_uniaxial'       :15.0                         # GaSb value
       }
       ,'X':{
          'mass_l'                :1.357                        # Vurgaftman1
          ,'mass_t'                :0.123                        # Vurgaftman1
          ,'bandgap'               :1.696                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.39e-3                      # Vurgaftman1
          ,'bandgap_beta'          :140                          # Vurgaftman1
          ,'defpot_absolute'       :1.91                         # A. Zunger: a_c :a_v + a_gap :0.73 + 1.18 :1.91
          ,'defpot_uniaxial'       :6.0                          # Munoz
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.385                              # A. Zunger
        
        ,'HH':{ 'mass'          :0.8   }                            # Landolt-Boernstein (m_p,h)
        ,'LH':{ 'mass'          :0.13  }                            # Landolt-Boernstein (m_p,l)
        ,'SO':{ 'mass'          :0.22  }                            # Vurgaftman1
        
        ,'defpot_absolute'   : 0.73                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 1.4                               # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -1.4
        ,'defpot_uniaxial_b' :-1.35   ,'defpot_uniaxial_d' :-4.3   # Vurgaftman1
       
        ,'delta_SO'          :0.676                              # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :5.18   gamma2 :1.19   gamma3 :1.97          # Vurgaftman1
        'L' :-10.94      ,'M' :-3.80       ,'N' :-11.82     
        ,'kappa' :0.31                                           # P. Lawaetz, PRB 4, 3460 (1971)          
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :2.386
        'S' :-0.12                                              # S :1 + 2F :1 + 2 (-0.56) :-0.12 (Vurgaftman1)
        ,'E_P' :18.7                                             # Vurgaftman1
        ,'B' : 0                   
        ,'L' :-3.779   ,'M' :-3.80   ,'N' :-4.659    
        ,'kappa' :-0.884           
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :250  ,'exponent' :1.9 }               # R.J.Stirn & W.M.Becker - Electron Mobility in Aluminum Antimonide (,'exponent' from InP)
        ,'holes':{      'mumax' :375  ,'exponent' :1.2 }               # D.Shaw & H.D. McKell - Tantalum doping and high resistivity in aluminium antimonide (,'exponent' from InP)
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # SIMBA
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # SIMBA
        }
        ,'Auger':{     'c_n' :1.0e-31       ,'c_p' :1.0e-31           # SIMBA
        }
        ,'radiative':{ 'c' :0                                       # ?
        }
    }
}
#

######### indium antimonide ###########################################
,'InSb':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :6.4794                                   # Vurgaftman1 (300 K)
        ,'a_expansion' :3.48e-5                                  # Vurgaftman1
    }
    
    ,'dielectric_consts':{
        'static_a'  :17.5                                       # Landolt-Boernstein epsilon(0)
        ,'optical_a' :15.68                                      # Landolt-Boernstein epsilon(infinity)
    }

    ,'elastic_consts':{
        'c11' :68.47   ,'c12' :37.35   ,'c44' :31.11                # Vurgaftman1
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :-0.076                                           # calculated by      S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.069                                           # experimental value S. Gironcoli et al., PRL 62(24), 2853 (1989)
      # 'e14' :-0.0717                                          # Landolt-Boernstein
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.0135                       # Vurgaftman1
          ,'bandgap'               :0.235                        # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.32e-3                      # Vurgaftman1
          ,'bandgap_beta'          :170                          # Vurgaftman1
        # ,'defpot_absolute'       :-6.94                        # Vurgaftman1
          ,'defpot_absolute'       :-6.04                        # Qteish / A. Zunger: a_c :a_v + a_gap :0.31 - 6.35 :-6.04
        #  ,'g'                     :-50                          # R.A. Isaacson: Phys. Rev. 169, 312 (1968)
          ,'g'                     :-51.6                        # J.-M. Jancu, R. Scholz, PRB 72, 193201 (2005)
       }
       ,'L':{ 
          'mass_l'                :0.25                         # Vurgaftman1 m_DOS
          ,'mass_t'                :0.25                         # Vurgaftman1 m_DOS
          ,'bandgap'               :0.93                         # Vurgaftman1 (0 K)
          ,'bandgap_alpha'         :0.32e-3                      # GAMMA value
          ,'bandgap_beta'          :170                          # GAMMA value
          ,'defpot_absolute'       :-2.20                        # Qteish / A. Zunger: a_c :a_v + a_gap :0.31 - 2.51 :-2.20
          ,'defpot_uniaxial'       :15.0                         # GaSb value
       }
       ,'X':{
          'mass_l'                :1.51                         # GaSb value
          ,'mass_t'                :0.22                         # GaSb value
          ,'bandgap'               :1.63                         # Vurgaftman1 (0 K): 0.63 eV but private communication I. Vurgaftman (0 K): 1.63 eV
          ,'bandgap_alpha'         :0.32e-3                      # GAMMA value
          ,'bandgap_beta'          :170                          # GAMMA value
          ,'defpot_absolute'       :1.41                         # Qteish / A. Zunger: a_c :a_v + a_gap :0.31 + 1.10 :1.41
          ,'defpot_uniaxial'       :4.53                         # Munoz
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :1.750                              # A. Zunger
        
        ,'HH':{ 'mass'          :0.405  }                           # Landolt-Boernstein (m_p,h)
        ,'LH':{ 'mass'          :0.0162 }                           # Landolt-Boernstein (m_p,l)
        ,'SO':{ 'mass'          :0.11   }                           # Vurgaftman1

      # ,'HH':{ 'mass'          :0.43   }                           # http://www.ioffe.ru/SVA/NSM/Semicond/InSb/bandstr.html
      # ,'LH':{ 'mass'          :0.015  }                           # http://www.ioffe.ru/SVA/NSM/Semicond/InSb/bandstr.html
      # ,'SO':{ 'mass'          :0.19   }                           # http://www.ioffe.ru/SVA/NSM/Semicond/InSb/bandstr.html
        
        ,'defpot_absolute'   : 0.31                              # Qteish value, no Zunger value available: a_v
      # ,'defpot_absolute'   : 0.36                              # Vurgaftman1 - Note that Vurgaftman1 has different sign convention. => -0.36
        ,'defpot_uniaxial_b' :-2.0   ,'defpot_uniaxial_d' :-4.7    # Vurgaftman1
       
        ,'delta_SO'          :0.81                               # Vurgaftman1
       
    }

    ,'kp_6_bands':{
      # gamma1 :34.8   gamma2 :15.5   gamma3 :16.5          # Vurgaftman1
        'L' :-97.8       ,'M' :-4.8        ,'N' :-99.0   
        #,'kappa' :14.76                                          # P. Lawaetz, PRB 4, 3460 (1971) 
        ,'kappa' :15.6                                           # Landolt-Börnstein           
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :0.235
        'S' :  0.54                                             # S :1 + 2F :1 + 2 (-0.23) :-0.54 (Vurgaftman1)
        ,'E_P' :23.3                                             # Vurgaftman1
        ,'B' :  0                   
        ,'L' :-51.661   ,'M' :-4.8   ,'N' :-52.861   
        ,'kappa' :7.91            
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :70400 ,'exponent' :1.9 }             # D.L.Rode - Electron Transport in InSb, InAs and InP (,'exponent' from InP)
        ,'holes':{      'mumax' :750   ,'exponent' :1.2 }             # R.K.Willardson & A.C.Beer - Semiconductors and Semimetals (,'exponent' from InP)
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # GaSb value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # GaSb value !!!
        }
        ,'Auger':{     'c_n' :1.0e-31       ,'c_p' :1.0e-31           # GaSb value !!!
        }
        ,'radiative':{ 'c' :1.0e-9                                  # GaSb value !!!
        }
    }
}
#

######### gallium nitride (zincblende) ################################
,'GaN_zb':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :4.50                                     # Vurgaftman1/Vurgaftman2 (300 K)
        ,'a_expansion' :0                                        # ?
    }
    
    ,'dielectric_consts':{
        'static_a'  :9.7                                        # http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/optic.html  Bougrov et al. (2001)  (300 K)
        ,'optical_a' :5.3                                        # high frequency (300 K)
    }

    ,'elastic_consts':{
        'c11' :293   ,'c12' :159   ,'c44' :155                      # Vurgaftman1/Vurgaftman2
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :0.56                                             # taken from M. Shur http://nina.ecse.rpi.edu/shur/Tutorial/GaNtutorial1/sld038.htm
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.15                         # Vurgaftman1/Vurgaftman2
          ,'bandgap'               :3.299                        # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.593e-3                     # Vurgaftman2
          ,'bandgap_beta'          :600                          # Vurgaftman2
        # ,'defpot_absolute'       :-6.71                        # Vurgaftman2
          ,'defpot_absolute'       :-6.68                        # A. Zunger: a_c :a_v + a_gap :0.69 - 7.37 :-6.68
       }
       ,'L':{ 
          'mass_l'                :0.2                          # ?
          ,'mass_t'                :0.2                          # ?
          ,'bandgap'               :5.59                         # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.593e-3                     # Vurgaftman2
          ,'bandgap_beta'          :600                          # Vurgaftman2
          ,'defpot_absolute'       :-7.46                        # A. Zunger: a_c :a_v + a_gap :0.69 - 8.15 :-7.46
          ,'defpot_uniaxial'       :14.26                        # GaAs value
       }
       ,'X':{
          'mass_l'                :0.5                          # Vurgaftman1/Vurgaftman2
          ,'mass_t'                :0.3                          # Vurgaftman1/Vurgaftman2
          ,'bandgap'               :4.52                         # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.593e-3                     # Vurgaftman2
          ,'bandgap_beta'          :600                          # Vurgaftman2
          ,'defpot_absolute'       :-0.52                        # A. Zunger: a_c :a_v + a_gap :0.69 - 1.21 :-0.52
          ,'defpot_uniaxial'       :6.5                          # GaAs value
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-0.726                             # A. Zunger
        
        ,'HH':{ 'mass'          :1.3   }                            # http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/bandstr.html Leszczynski et al. (1996), Fan et al. (1996)
        ,'LH':{ 'mass'          :0.19  }                            # http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/bandstr.html Leszczynski et al. (1996), Fan et al. (1996)
        ,'SO':{ 'mass'          :0.29  }                            # Vurgaftman1/Vurgaftman2
        
        ,'defpot_absolute'   : 0.69                              # A. Zunger/Vurgaftman2: a_v  -  Note that Vurgaftman1/Vurgaftman2 has different sign convention. => -0.69
        ,'defpot_uniaxial_b' :-2.0   ,'defpot_uniaxial_d' :-3.7    # Vurgaftman2
       
        ,'delta_SO'          :0.017                              # Vurgaftman1/Vurgaftman2
       
    }

    ,'kp_6_bands':{
      # gamma1 :2.70   gamma2 :0.76   gamma3 :1.11          # Vurgaftman2
        'L' :-6.74       ,'M' :-2.18       ,'N' :-6.66              
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :3.299
        'S' :-0.90                                              # S :1 + 2F :1 + 2 (-0.95) :-0.90 (Vurgaftman2)
        ,'E_P' :25.0                                             # Vurgaftman1/Vurgaftman2
        ,'B' : 0                   
        ,'L' : 0.825   ,'M' :-2.18   ,'N' :0.905               
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         #
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         #
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :0                                       # ?
        }
    }
}
#

######### aluminum nitride (zincblende) ###############################
,'AlN_zb':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :4.38                                     # Vurgaftman1/Vurgaftman2 (300 K)
        ,'a_expansion' :0                                        # ?
    }
    
     ,'dielectric_consts':{
        'static_a'  :9.14                                       # ? Landolt-Boernstein epsilon(0) wurtzite, Collins et al. (1967)
        ,'optical_a' :4.84                                       # ? Landolt-Boernstein epsilon(infinity) wurtzite, Collins et al. (1967)
    }

    ,'elastic_consts':{
        'c11' :304   ,'c12' :160   ,'c44' :193                      # Vurgaftman1/Vurgaftman2
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :0.92                                             # ? taken from M. Shur http://nina.ecse.rpi.edu/shur/Tutorial/GaNtutorial1/sld038.htm (wurtzite)
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.25                         # Vurgaftman1/Vurgaftman2
          ,'bandgap'               :5.4                          # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.593e-3                     # Vurgaftman2
          ,'bandgap_beta'          :600                          # Vurgaftman2
        # ,'defpot_absolute'       :-4.5                         # Vurgaftman2
          ,'defpot_absolute'       :-5.22                        # A. Zunger: a_c :a_v + a_gap :4.94 - 10.16 :-5.22
       }
       ,'L':{ 
          'mass_l'                :0.2                          # ?
          ,'mass_t'                :0.2                          # ?
          ,'bandgap'               :9.3                          # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.593e-3                     # Vurgaftman2
          ,'bandgap_beta'          :600                          # Vurgaftman2
          ,'defpot_absolute'       :-4.95                        # A. Zunger: a_c :a_v + a_gap :4.94 - 9.89 :-4.95
          ,'defpot_uniaxial'       :14.26                        # GaAs value
       }
       ,'X':{
          'mass_l'                :0.53                         # Vurgaftman1/Vurgaftman2
          ,'mass_t'                :0.31                         # Vurgaftman1/Vurgaftman2
          ,'bandgap'               :4.9                          # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.593e-3                     # Vurgaftman2
          ,'bandgap_beta'          :600                          # Vurgaftman2
          ,'defpot_absolute'       :3.81                         # A. Zunger: a_c :a_v + a_gap :4.94 - 1.13 :-3.81
          ,'defpot_uniaxial'       :6.5                          # GaAs value
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-1.526                             # A. Zunger
        
        ,'HH':{ 'mass'          :0.3   }                            # ?
        ,'LH':{ 'mass'          :0.3   }                            # ?
        ,'SO':{ 'mass'          :0.47  }                            # Vurgaftman1/Vurgaftman2
        
        ,'defpot_absolute'   : 4.94                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 4.9                               # Vurgaftman2 - Note that Vurgaftman1/Vurgaftman2 has different sign convention. => -4.9
        ,'defpot_uniaxial_b' :-1.7   ,'defpot_uniaxial_d' :-5.5    # Vurgaftman2
       
        ,'delta_SO'          :0.019                              # Vurgaftman1/Vurgaftman2
       
    }

    ,'kp_6_bands':{
      # gamma1 :1.92   gamma2 :0.47   gamma3 :0.85          # Vurgaftman1/Vurgaftman2
        'L' :-4.80       ,'M' :-1.98       ,'N' :-5.10              
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :5.4 (Vurgaftman2)
        'S' :-1.02                                              # S :1 + 2F :1 + 2 (-1.01) :-1.02 (Vurgaftman2)
        ,'E_P' :27.1                                             # Vurgaftman1/Vurgaftman2
        ,'B' : 0                   
        ,'L' : 0.213   ,'M' :-1.98   ,'N' :-0.087               
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         #
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         #
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :0                                       # ?
        }
    }
}
#

######### indium nitride (zincblende) #################################
,'InN_zb':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :4.98                                     # Vurgaftman1/Vurgaftman2 (300 K)
        ,'a_expansion' :0                                        # ?
    }
    
     ,'dielectric_consts':{
        'static_a'  :15.3                                       # ? wurtzite Zubrilov (2001)
        ,'optical_a' : 9.3                                       # ? Landolt-Boernstein epsilon(infinity) wurtzite
    }

    ,'elastic_consts':{
        'c11' :187   ,'c12' :125   ,'c44' :86                       # Vurgaftman1/Vurgaftman2
    }                    
                                             
    ,'piezoelectric_consts':{
        'e14' :0.37                                             # taken from M. Shur http://nina.ecse.rpi.edu/shur/Tutorial/GaNtutorial1/sld038.htm
    }                                                
   
    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.07                         # Vurgaftman2
          ,'bandgap'               :0.78                         # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.245e-3                     # Vurgaftman2
          ,'bandgap_beta'          :624                          # Vurgaftman2
        # ,'defpot_absolute'       :-2.65                        # Vurgaftman2
          ,'defpot_absolute'       :-2.93                        # A. Zunger: a_c :a_v + a_gap :0.73 - 3.66 :-2.93
       }
       ,'L':{ 
          'mass_l'                :0.2                          # ?
          ,'mass_t'                :0.2                          # ?
          ,'bandgap'               :5.82                         # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.245e-3                     # Vurgaftman2
          ,'bandgap_beta'          :624                          # Vurgaftman2
          ,'defpot_absolute'       :-4.50                        # A. Zunger: a_c :a_v + a_gap :0.73 - 5.23 :-4.50
          ,'defpot_uniaxial'       :14.26                        # GaAs value
       }
       ,'X':{
          'mass_l'                :0.48                         # Vurgaftman1/Vurgaftman2
          ,'mass_t'                :0.27                         # Vurgaftman1/Vurgaftman2
          ,'bandgap'               :2.51                          # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.245e-3                     # Vurgaftman2
          ,'bandgap_beta'          :624                          # Vurgaftman2
          ,'defpot_absolute'       :-0.62                        # A. Zunger: a_c :a_v + a_gap :0.73 - 1.35 :-0.62
          ,'defpot_uniaxial'       :6.5                          # GaAs value
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-0.462                             # A. Zunger
        
        ,'HH':{ 'mass'          :0.2   }                            # ?
        ,'LH':{ 'mass'          :0.2   }                            # ?
        ,'SO':{ 'mass'          :0.3   }                            # Vurgaftman1/Vurgaftman2
        
        ,'defpot_absolute'   : 0.73                              # A. Zunger: a_v
      # ,'defpot_absolute'   : 0.7                               # Vurgaftman2 - Note that Vurgaftman2 has different sign convention. => -0.7
        ,'defpot_uniaxial_b' :-1.2   ,'defpot_uniaxial_d' :-9.3    # Vurgaftman1/Vurgaftman2
       
        ,'delta_SO'          :0.005                              # Vurgaftman2
       
    }

    ,'kp_6_bands':{
      # gamma1 :3.72   gamma2 :1.26   gamma3 :1.63          # Vurgaftman1/Vurgaftman2
        'L' :-9.76       ,'M' :-2.20       ,'N' :-9.78              
    }                    
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :0.78 (Vurgaftman2)
        'S' :-7.72                                              # S :1 + 2F :1 + 2 (-4.36) :-7.72 (Vurgaftman2)
        ,'E_P' :17.2                                             # Vurgaftman2
        ,'B' : 0                   
        ,'L' : 12.24   ,'M' :-2.20   ,'N' :12.224               
    }
    
    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         #
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         #
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :0                                       # ?
        }
    }
}
#

######### zinc selenide ############################################
,'ZnSe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :5.6674                                   # [Angstrom]    300 K, H. Karzel et al., PRB 53, 11425 (1996)
        ,'a_expansion' :7.7e-5                                   # [Angstrom/K]  LB
    }
    
    ,'dielectric_consts':{
        'static_a'  : 8.6                                       # S. Adachi et al., PRB 43, 9569 (1995)
        ,'optical_a' : 5.73                                      # S. Adachi et al., PRB 43, 9569 (1995), R.T. Senger et al., phys. stat. sol. (b) 241, 1896 (2004)
    }

    ,'elastic_consts':{
        'c11' :82.6   ,'c12' :49.8   ,'c44' :40.0                   # C. Van de Walle, PRB 39, 1871 (1989)
    }                    

    ,'piezoelectric_consts':{
        'e14' :0.049                                            # LB
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.145                        # Isshiki, Journal of Crystal Growth 86, 615 (1988)
#         ,'bandgap'               :2.71                         # 300 K, J. Piprek
          ,'bandgap'               :2.825                        #   0 K, J. Piprek
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
#         ,'defpot_absolute'       :-6.96                        # A. Zunger: a_c :a_v + a_gap :... + ... :-6.96, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_absolute'       :-5.93                        # S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.89                        # A. Zunger: a_c :a_v + a_gap :... + ... :-4.89, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.61                        # A. Zunger: a_c :a_v + a_gap :... + ... :-0.61, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.39                               # A. Zunger, average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :1.04  ,'g' :0 }                      # Isshiki, Journal of Crystal Growth 86, 615 (1988)
        ,'LH':{ 'mass'          :0.16  ,'g' :0 }                      # from Luttinger parameters in Venghaus 1979
        ,'SO':{ 'mass'          :0.30  }                            # P. Lawaetz, PRB 4, 3460 (1971)

        ,'defpot_absolute'   :-1.97                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-1.2   ,'defpot_uniaxial_d' :-4.3    # b (van de Walle and Cardona), d (Cardona) [eV]
       
        ,'delta_SO'          :0.45                               # Chelikowski
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.3    gamma2 :1.14   gamma3 :1.84          # J. Piprek
      # gamma1 :4.8    gamma2 :0.67   gamma3 :1.53          # Venghaus 1979
        'L' :-9.86       ,'M' :-3.02       ,'N' :-11.04             # calculated from J. Piprek parameters  gamma1 :4.3, gamma2 :1.14, gamma3 :1.84
        ,'kappa' :0.64                                           # P. Lawaetz, PRB 4, 3460 (1971)
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :2.825
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :24.2                                             # P. Lawaetz, PRB 4, 3460 (1971)
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-3.02       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#

######### zinc telluride ############################################
,'ZnTe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :6.07                                     # [Angstrom]    300 K (?), T. Magalingham et al., Cryst. Res. Technol. 37, 329 (2002)
        ,'a_expansion' :7.7e-5                                   # [Angstrom/K]  ZnSe
    }
    
    ,'dielectric_consts':{
        'static_a'  : 9.67                                      # D. L. Rode, PRB 2, 4036 (1970)
        ,'optical_a' : 7.23                                      # D. L. Rode, PRB 2, 4036 (1970)
    }

    ,'elastic_consts':{
        'c11' :72.2   ,'c12' :40.9   ,'c44' :30.8                   # M. Yamada et al., J. Phys. D: Appl. Phys. 10, 1309 (1977)
    }                    

    ,'piezoelectric_consts':{
        'e14' :0                                                # ???
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.21                         # M. L. Cohen et al., PR 141, 789 (1966)
          ,'bandgap'               :2.34                         #   0 K (?), D. L. Rode, PRB 2, 4036 (1970)
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-7.88                        # A. Zunger: a_c :a_v + a_gap :... + ... :-7.88, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # R. Brazis et al., phys. stat. sol. (b) 244 1662 (2007)
          ,'mass_t'                :0.5                          # R. Brazis et al., phys. stat. sol. (b) 244 1662 (2007)
          ,'bandgap'               :3.6                          # calculated from Gamma-L valley separation energy (R. Brazis)
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.68                        # A. Zunger: a_c :a_v + a_gap :... + ... :-4.68, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.56                        # A. Zunger: a_c :a_v + a_gap :... + ... :-0.56, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.96                               # A. Zunger, average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :0.2   ,'g' :0 }                      # Singh, p. 840 (1993)
        ,'LH':{ 'mass'          :0.16  ,'g' :0 }                      # ZnSe
        ,'SO':{ 'mass'          :0.33  }                            # P. Lawaetz, PRB 4, 3460 (1971)

        ,'defpot_absolute'   :-2.28                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-1.26  ,'defpot_uniaxial_d' :-4.6    # b (van de Walle), d (experiment) [eV]
       
        ,'delta_SO'          :0.91                               # Qteish/Needs
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.00   gamma2 :1.15   gamma3 :1.29          # R. L. Hollis, PRB 15, 932 (1977) (calculations)
        'L' :-9.60       ,'M' :-2.70       ,'N' :-7.74              # calculated from Luttinger parameters of R. L. Hollis, PRB 15, 932 (1977)
        ,'kappa' :0.42                                           # P. Lawaetz, PRB 4, 3460 (1971)
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :2.34
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :19.1                                             # P. Lawaetz, PRB 4, 3460 (1971)
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-2.70       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#

######### magnesium selenide ############################################
,'MgSe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :5.91                                     # [Angstrom]    300 K, B. Jobst et al., APL 69, 97 (1996)
        ,'a_expansion' :7.7e-5                                   # [Angstrom/K]  ??? ZnSe
    }
    
    ,'dielectric_consts':{
        'static_a'  : 8.6                                       # ??? ZnSe
        ,'optical_a' : 4.86                                      # S. Saib et al., Eur. Phys. J. B 73, 185 (2010)
    }

    ,'elastic_consts':{
        'c11' :61.8   ,'c12' :45.1   ,'c44' :22.8                   # S. Saib et al., Eur. Phys. J. B 73, 185 (2010)
    }                    

    ,'piezoelectric_consts':{
        'e14' :0.68955                                          # S. Saib et al., Eur. Phys. J. B 73, 185 (2010)
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.145                        # ??? ZnSe
          ,'bandgap'               :4.0                          # 300 K, B. Jobst et al., APL 69, 97 (1996)
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-5.93                        # ??? ZnSe S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.89                        # ??? ZnSe A. Zunger: a_c :a_v + a_gap :... + ... :...
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.61                        # ??? ZnSe A. Zunger: a_c :a_v + a_gap :... + ... :...
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0          # should be adjusted!!!  # ??? A. Zunger, average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :1.04  ,'g' :0 }                      # ??? ZnSe
        ,'LH':{ 'mass'          :0.16  ,'g' :0 }                      # ??? ZnSe
        ,'SO':{ 'mass'          :0.30  }                            # ??? ZnSe

        ,'defpot_absolute'   :-1.97                              # ??? ZnSe
        ,'defpot_uniaxial_b' :-1.2   ,'defpot_uniaxial_d' :-4.3    # ??? ZnSe
       
        ,'delta_SO'          :0.45                               # ??? ZnSe
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.8    gamma2 :0.67   gamma3 :1.53          # ??? ZnSe
        'L' :-9.86       ,'M' :-3.02       ,'N' :-11.04             # ??? ZnSe
        ,'kappa' :0.64                                           # ??? ZnSe
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :4.0
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :24.2                                             # ??? ZnSe
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-3.02       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#

######### cadmium selenide ############################################
,'CdSe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :6.052                                    # [Angstrom]    300 K, Wei and Zunger
        ,'a_expansion' :7.7e-5                                   # [Angstrom/K]  ??? ZnSe
    }
    
    ,'dielectric_consts':{
        'static_a'  : 9.7                                       # LB
        ,'optical_a' : 5.8                                       # LB
    }

    ,'elastic_consts':{
        'c11' :74.6   ,'c12' :46.1   ,'c44' :13.0                   # LB
    }                    

    ,'piezoelectric_consts':{
        'e14' :0.049                                            # ZnSe value ?
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.12                         # Furdyna PRB 49
#         ,'bandgap'               :1.766                        # 300 K, Furdyna PRB 49 (theory)
          ,'bandgap'               :1.846                        #   0 K, J. Piprek
          ,'bandgap_alpha'         :6.96e-4                      # Furdyna PRB 50
          ,'bandgap_beta'          :281                          # Furdyna PRB 50
          ,'defpot_absolute'       :-3.77                        # A. Zunger: a_c :a_v + a_gap :... + ... :-3.77, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :3.7                          # Furdyna PRB 50
          ,'bandgap_alpha'         :6.96e-4                      # Furdyna PRB 50
          ,'bandgap_beta'          :281                          # Furdyna PRB 50
          ,'defpot_absolute'       :-4.89                        # ??? ZnSe A. Zunger: a_c :a_v + a_gap :... + ... :-4.89, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :4.8                          # Furdyna PRB 50
          ,'bandgap_alpha'         :6.96e-4                      # Furdyna PRB 50
          ,'bandgap_beta'          :281                          # Furdyna PRB 50
          ,'defpot_absolute'       :-0.61                        # ??? ZnSe A. Zunger: a_c :a_v + a_gap :... + ... :-0.61, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.7162                             # 30 % VBO to ZnSe, Egap :1.766 eV, average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :0.9   ,'g' :0 }                      # Furdyna PRB 49
        ,'LH':{ 'mass'          :0.18  ,'g' :0 }                      # Furdyna PRB 49
        ,'SO':{ 'mass'          :0.34  }                            # Furdyna PRB 49

        ,'defpot_absolute'   :-1.81                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-1.2   ,'defpot_uniaxial_d' :-4.3    # ??? ZnSe
       
        ,'delta_SO'          :0.42                               # Furdyna PRB 50
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.8    gamma2 :0.67   gamma3 :1.53          # ??? ZnSe
        'L' :-9.86       ,'M' :-3.02       ,'N' :-11.04             # ??? ZnSe
        ,'kappa' :0                                              # ???
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :1.776
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :24.2                                             # ???
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-3.02       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#

######### beryllium selenide ############################################
,'BeSe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :5.139                                    # [Angstrom]    300 K (?) epitaxy dp
        ,'a_expansion' :7.7e-5                                   # [Angstrom/K]  ??? ZnSe
    }
    
    ,'dielectric_consts':{
        'static_a'  : 8.6                                       # ??? ZnSe
        ,'optical_a' : 5.73                                      # ??? ZnSe
    }

    ,'elastic_consts':{
        'c11' :82.6   ,'c12' :49.8   ,'c44' :40.0                   # ??? ZnSe
    }                    

    ,'piezoelectric_consts':{
        'e14' :0.049                                            # ??? ZnSe
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.145                        # ???
#         ,'bandgap'               :...                          # 300 K
          ,'bandgap'               :5.15                         #   0 K
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-5.93                        # ??? ZnSe
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :6                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.89                        # ??? ZnSe A. Zunger: a_c :a_v + a_gap :... + ... :-4.89, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :6                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.61                        # ??? ZnSe A. Zunger: a_c :a_v + a_gap :... + ... :-0.61, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-0.45                              # average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :1.04  ,'g' :0 }                      # ??? ZnSe
        ,'LH':{ 'mass'          :0.16  ,'g' :0 }                      # ??? ZnSe
        ,'SO':{ 'mass'          :0.30  }                            # ??? ZnSe

        ,'defpot_absolute'   :-1.97                              # ??? ZnSe
        ,'defpot_uniaxial_b' :-1.2   ,'defpot_uniaxial_d' :-4.3    # ??? ZnSe
       
        ,'delta_SO'          :0.45                               # ??? ZnSe
    }

    ,'kp_6_bands':{
      # gamma1 :4.8    gamma2 :0.67   gamma3 :1.53          # ??? ZnSe
        'L' :-9.86       ,'M' :-3.02       ,'N' :-11.04             # ??? ZnSe
        ,'kappa' :0.64                                           # ??? ZnSe
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :5.15
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :24.2                                             # ??? ZnSe
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-3.02       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#

######### cadmium telluride ############################################
,'CdTe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :6.486                                    # [Angstrom] 300 K - Book: J. Piprek
        ,'a_expansion' :3.11e-5                                  # [Angstrom/K]       Book: J. Piprek - [see also R.D. Greenough et al., J. Phys. D: Appl. Phys. 6, 587 (1973) for ,'a' discussion of the temperature dependence)
    }
    
    ,'dielectric_consts':{
        'static_a'  :10.6                                       #                  V. Borsari, C. Jacoboni, phys. stat. sol. (b) 54, 649 (1972)
        ,'optical_a' : 7.13                                      # high frequency - V. Borsari, C. Jacoboni, phys. stat. sol. (b) 54, 649 (1972)
    }

    ,'elastic_consts':{
        'c11' :53.8   ,'c12' :37.4   ,'c44' :20.18                  # [GPa] R.D. Greenough et al., J. Phys. D: Appl. Phys. 6, 587 (1973) at 298 K
    }                    

    ,'piezoelectric_consts':{
        'e14' :0                                                # ???
    }                                                

# E_gap :1.528 eV (300 K)
# E_gap :1.606 eV (4.2 K)
# E_gap :1.43  eV (Gamma, direct semiconductor, T :300 K) J. Piprek
# E_gap :1.606 eV (Gamma, direct semiconductor, T :  0 K) J. Piprek
# E_gap :1.425 eV (Gamma, direct semiconductor, T :300 K) J. Reno et al., APL 49, 106 (1986)
# E_gap :1.550 eV (Gamma, direct semiconductor, T : 77 K)
# E_gap :1.600 eV (Gamma, direct semiconductor, T :  4 K) J. Reno et al., APL 49, 106 (1986)
# E_gap :1.54  eV (Gamma, direct semiconductor, T :?   K) D.L. Rode, PRB 2, 4036 (1970)
# E_gap :3.04  eV (L                          , T :?   K) calculated from Gamma-L valley separation energy (R. Brazis)
    ,'conduction_bands':{
       'Gamma':{ 
#         'mass'                  :0.0963                       # V. Borsari, C. Jacoboni, phys. stat. sol. (b) 54, 649 (1972)
          'mass'                  :0.08992                      # calculated from Novik parameters
          ,'bandgap'               :1.606                        # 4.2 K
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-5.84                        # A. Zunger: a_c :a_v + a_gap :-2.14 + ... :-5.84
          ,'g'                     :-0.30                        # 280 K, M. Oestreich et al., PRB 53, 7911 (1996)
       }
       ,'L':{ 
          'mass_l'                :0.5                          # V. Borsari, C. Jacoboni, phys. stat. sol. (b) 54, 649 (1972)
          ,'mass_t'                :0.5                          # V. Borsari, C. Jacoboni, phys. stat. sol. (b) 54, 649 (1972)
          ,'bandgap'               :3.04                         # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.02                        # A. Zunger: a_c :a_v + a_gap :-2.14 + ... :-4.02
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :3.5                          # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.70                        # A. Zunger: a_c :a_v + a_gap :-2.14 + ... :-0.70
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
      # 'bandoffset'        :0.86                               # A. Zunger, average ,'valence' band energy E_v,av [eV]
      # 'bandoffset'        :0.86 - 1.430328571                 # VBO HgTe :570 meV
      # 'bandoffset'        :      -0.570328571                 # VBO HgTe :570 meV
      #  'bandoffset'       :-0.570 + 0.91/3 - 1.08/3
      #                     -0.570 + 0.056666666               # VBO HgTe :570 meV
        'bandoffset'        :-0.513333333                       # VBO HgTe :570 meV
        
        ,'HH':{ 'mass'          :0.72  ,'g' :0 }                      # J. Piprek, p. 21   ,'g' :???
        ,'LH':{ 'mass'          :0.13  ,'g' :0 }                      # J. Piprek, p. 21   ,'g' :???
        ,'SO':{ 'mass'          :0.28  }                            # [LawaetzPRB1971]

        ,'defpot_absolute'   :-2.14                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-1.1   ,'defpot_uniaxial_d' :-2.8    # b (van de Walle), d (experiment) [eV]
       
      # ,'delta_SO'          :0.93                               # J. Piprek
        ,'delta_SO'          :0.91                               # Novik
       
    }

    ,'kp_6_bands':{
      # gamma1 :5.3     gamma2 :1.7     gamma3 :2.0         # J. Piprek
      # gamma1 :5.37203 gamma2 :1.67102 gamma3 :1.98102     # calculated from Novik parameters
        'L' :-13.0561    ,'M' : -3.03      ,'N' :-11.8861           # calculated from Novik parameters
      # ,'kappa' :1.27                                           # P. Lawaetz, PRB 4, 3460 (1971)
        ,'kappa' :0.641017                                       # calculated from Novik parameters
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :1.519
        'S' :0.82                                               # S :1 + 2F   Novik: F :-0.09
      # ,'E_P' :18.0                                             # Ref. E_P: J. Reno et al., APL 49, 106 (1986)
        ,'E_P' :18.8                                             # Novik
        ,'B' : 0                     
      # gamma1 :1.47   gamma2 :-0.28  gamma3 :0.03          # Novik
        ,'L' :-1.35       ,'M' : -3.03      ,'N' :-0.18              # calculated from Novik parameters
      # ,'kappa' :1.27                                           # P. Lawaetz, PRB 4, 3460 (1971)
        ,'kappa' :-1.31                                          # Novik
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }


}
#

######### mercury telluride ############################################
,'HgTe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :6.486                                    # [Angstrom]     CdTe value
        ,'a_expansion' :3.11e-5                                  # [Angstrom/K]   CdTe value
    }
    
    ,'dielectric_consts':{
        'static_a'  :10.6                                       # CdTe value
        ,'optical_a' : 7.13                                      # CdTe value
    }

    ,'elastic_consts':{
        'c11' :53.8   ,'c12' :37.4   ,'c44' :20.18                  # CdTe value
    }                    

    ,'piezoelectric_consts':{
        'e14' :0                                                # ???
    }                                                

# E_gap :-0.122 eV (Gamma, semimetal, T :300 K) J. Reno et al., APL 49, 106 (1986)
# E_gap :-0.261 eV (Gamma, semimetal, T : 77 K)
# E_gap :-0.302 eV (Gamma, semimetal, T :  4 K) J. Reno et al., APL 49, 106 (1986)
    ,'conduction_bands':{
       'Gamma':{ 
#         'mass'                  :-0.031                       # [LawaetzPRB1971]
          'mass'                  :-0.03096                     # calculated from Novik parameters
#         ,'bandgap'               :-0.302                       # 4 K   J. Reno et al., APL 49, 106 (1986)
          ,'bandgap'               :-0.303                       # 4 K   Novik
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-6.64                        # A. Zunger: a_c :a_v + a_gap :-3.45 + ... :-6.64
          ,'g'                     :-0.30                        # 280 K, M. Oestreich et al., PRB 53, 7911 (1996)
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :3                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.19                        # A. Zunger: a_c :a_v + a_gap :-3.45 + ... :-4.19
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :3.5                          # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-1.48                        # A. Zunger: a_c :a_v + a_gap :-3.45 + ... :-1.48
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
      # 'bandoffset'        :1.18                               # A. Zunger, average ,'valence' band energy E_v,av [eV]
        'bandoffset'        :0.0                                #
        
        ,'HH':{ 'mass'          :1.12  ,'g' :0 }                      # [LawaetzPRB1971]
        ,'LH':{ 'mass'          :-0.026 ,'g' :0 }                      # [LawaetzPRB1971]
        ,'SO':{ 'mass'          :0.102 }                            # [LawaetzPRB1971]

        ,'defpot_absolute'   :-3.45                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-1.15  ,'defpot_uniaxial_d' :-2.8    # b (van de Walle), d (CdTe) [eV]
       
      # ,'delta_SO'          :1.05                               # Qteish/Needs
        ,'delta_SO'          :1.08                               # Novik
       
    }

    ,'kp_6_bands':{
      # gamma1 :-15.5  gamma2 :-8.9   gamma3 :-8.9          #   4 K -  J. Reno et al., APL 49, 106 (1986): gamma2 :gamma3
      # gamma1 :-44.8  gamma2 :-23.55 gamma3 :-23.55        # 300 K -  J. Reno et al., APL 49, 106 (1986): gamma2 :gamma3
      # gamma1 =-16.5821 gamma2 =-9.84103 gamma3 =-9.04103     # calculated from Novik parameters
      # 'L' :50.1        ,'M' :-3.3        ,'N' :53.4               #   4 K -  J. Reno et al., APL 49, 106 (1986): gamma2 :gamma3
        'L' :54.9462	,'M' :-4.1        ,'N' :54.2462            # calculated from Novik parameters
#       ,'kappa' :-10.85                                         # P. Lawaetz, PRB 4, 3460 (1971)
        ,'kappa' :-10.741                                        # calculated from Novik parameters
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :1.519
        'S' :1.0                                                # S :1 + 2F   Novik: F :0
      # ,'E_P' :18.0                                             # Ref. E_P: J. Reno et al., APL 49, 106 (1986)
        ,'E_P' :18.8                                             # Novik
        ,'B' : 0                     
      # gamma1 :4.1    gamma2 :0.5    gamma3 :1.3           # Novik
        ,'L' :-7.1        ,'M' :-4.1        ,'N' :-7.8               # calculated from Novik parameters
        ,'kappa' :-0.4                                           # Novik
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }


}
#

######### magnesium telluride ############################################
# MgTe has normally wurtzite crystal structure! ==> This set of parameters is for Cd(1-x)Mg(x)Te alloys.
,'MgTe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :6.435                                    # [Angstrom] 300 K
        ,'a_expansion' :3.11e-5                                  # ??? CdTe [Angstrom/K]
    }
    
    ,'dielectric_consts':{
        'static_a'  :10.6                                       # ??? CdTe
        ,'optical_a' : 7.13                                      # ??? CdTe
    }

    ,'elastic_consts':{
        'c11' :53.8   ,'c12' :37.4   ,'c44' :20.18                  # [GPa] CdTe
    }                    

    ,'piezoelectric_consts':{
        'e14' :0                                                # ??? CdTe
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.08992                      # ??? CdTe
#         ,'bandgap'               :3.49                         # 300 K
          ,'bandgap'               :3.2                          # 4.2 K
          ,'bandgap_alpha'         :0                            # ??? CdTe
          ,'bandgap_beta'          :0                            # ??? CdTe
          ,'defpot_absolute'       :-5.84                        # ??? CdTe A. Zunger: a_c :a_v + a_gap :...
          ,'g'                     :-0.30                        # ??? CdTe
       }
       ,'L':{ 
          'mass_l'                :0.5                          #
          ,'mass_t'                :0.5                          #
          ,'bandgap'               :3.04                         # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.02                        # ??? A. Zunger: a_c :a_v + a_gap :...
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :3.5                          # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.70                        #  ??? A. Zunger: a_c :a_v + a_gap :...
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-0.513333333                       # ???
        
        ,'HH':{ 'mass'          :0.72  ,'g' :0 }                      # ???
        ,'LH':{ 'mass'          :0.13  ,'g' :0 }                      # ???
        ,'SO':{ 'mass'          :0.28  }                            # ???

        ,'defpot_absolute'   :-2.14                              # ???
        ,'defpot_uniaxial_b' :-1.1   ,'defpot_uniaxial_d' :-2.8    # ???
       
        ,'delta_SO'          :0.91                               # ???
       
    }

    ,'kp_6_bands':{
      # gamma1 :5.3    gamma2 :1.7    gamma3 :2.0           # ???
        'L' :-13.05      ,'M' : -3.03      ,'N' :-11.88             # ???
        ,'kappa' :1.27                                           # ???
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :3.49
        'S' :0.82                                               # ??? S :1 + 2F
        ,'E_P' :18.8                                             # ???
        ,'B' : 0                     
      # gamma1 :1.47   gamma2 :-0.28  gamma3 :0.03          # ???
        ,'L' :-13.05      ,'M' : -3.03      ,'N' :-11.88             # ???
        ,'kappa' :-1.31                                          # ???
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }


}
#

######### manganese selenide ############################################
# So far, everything is taken from ZnSe. Please add appropriate values and add proper references!!!
,'MnSe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :5.6674                                   # [Angstrom]    300 K, H. Karzel et al., PRB 53, 11425 (1996)
        ,'a_expansion' :7.7e-5                                   # [Angstrom/K]  LB
    }
    
    ,'dielectric_consts':{
        'static_a'  : 8.6                                       # S. Adachi et al., PRB 43, 9569 (1995)
        ,'optical_a' : 5.73                                      # S. Adachi et al., PRB 43, 9569 (1995), R.T. Senger et al., phys. stat. sol. (b) 241, 1896 (2004)
    }

    ,'elastic_consts':{
        'c11' :82.6   ,'c12' :49.8   ,'c44' :40.0                   # C. Van de Walle, PRB 39, 1871 (1989)
    }                    

    ,'piezoelectric_consts':{
        'e14' :0.049                                            # LB
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.145                        # Isshiki, Journal of Crystal Growth 86, 615 (1988)
#         ,'bandgap'               :2.71                         # 300 K, J. Piprek
          ,'bandgap'               :2.825                        #   0 K, J. Piprek
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
#         ,'defpot_absolute'       :-6.96                        # A. Zunger: a_c :a_v + a_gap :... + ... :-6.96, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_absolute'       :-5.93                        # S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.89                        # A. Zunger: a_c :a_v + a_gap :... + ... :-4.89, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.61                        # A. Zunger: a_c :a_v + a_gap :... + ... :-0.61, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.39                               # A. Zunger, average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :1.04  ,'g' :0 }                      # Isshiki, Journal of Crystal Growth 86, 615 (1988)
        ,'LH':{ 'mass'          :0.16  ,'g' :0 }                      # from Luttinger parameters in Venghaus 1979
        ,'SO':{ 'mass'          :0.30  }                            # P. Lawaetz, PRB 4, 3460 (1971)

        ,'defpot_absolute'   :-1.97                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-1.2   ,'defpot_uniaxial_d' :-4.3    # b (van de Walle and Cardona), d (Cardona) [eV]
       
        ,'delta_SO'          :0.45                               # Chelikowski
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.3    gamma2 :1.14   gamma3 :1.84          # J. Piprek
      # gamma1 :4.8    gamma2 :0.67   gamma3 :1.53          # Venghaus 1979
        'L' :-9.86       ,'M' :-3.02       ,'N' :-11.04             # calculated from J. Piprek parameters  gamma1 :4.3, gamma2 :1.14, gamma3 :1.84
        ,'kappa' :0.64                                           # P. Lawaetz, PRB 4, 3460 (1971)
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :2.825
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :24.2                                             # P. Lawaetz, PRB 4, 3460 (1971)
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-3.02       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#

######### manganese telluride ############################################
# MnTe has NiAs crystal structure! ==> This set of parameters is for Cd(1-x)Mn(x)Te alloys.
,'MnTe':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :6.33                                     # [Angstrom] 300 K (?), MnTe zincblende, J. Kossut, Acta Phys. Pol. A 100, 111 (2001)
        ,'a_expansion' :3.11e-5                                  # ??? CdTe [Angstrom/K]
    }
    
    ,'dielectric_consts':{
        'static_a'  :10.6                                       # ??? CdTe
        ,'optical_a' : 7.13                                      # ??? CdTe
    }

    ,'elastic_consts':{
        'c11' :53.8   ,'c12' :37.4   ,'c44' :20.18                  # [GPa] CdTe
    }                    

    ,'piezoelectric_consts':{
        'e14' :0                                                # ??? CdTe
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.08992                      # ??? CdTe
#         ,'bandgap'               :2.9                          # 300 K
          ,'bandgap'               :3.2                          # 4.2 K
          ,'bandgap_alpha'         :0                            # ??? CdTe
          ,'bandgap_beta'          :0                            # ??? CdTe
          ,'defpot_absolute'       :-5.84                        # ??? CdTe A. Zunger: a_c :a_v + a_gap :...
          ,'g'                     :-0.30                        # ??? CdTe
       }
       ,'L':{ 
          'mass_l'                :0.5                          #
          ,'mass_t'                :0.5                          #
          ,'bandgap'               :3.04                         # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.02                        # ??? A. Zunger: a_c :a_v + a_gap :...
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :3.5                          # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.70                        #  ??? A. Zunger: a_c :a_v + a_gap :...
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-0.513333333                       # ???
        
        ,'HH':{ 'mass'          :0.72  ,'g' :0 }                      # ???
        ,'LH':{ 'mass'          :0.13  ,'g' :0 }                      # ???
        ,'SO':{ 'mass'          :0.28  }                            # ???

        ,'defpot_absolute'   :-2.14                              # ???
        ,'defpot_uniaxial_b' :-1.1   ,'defpot_uniaxial_d' :-2.8    # ???
       
        ,'delta_SO'          :0.91                               # ???
       
    }

    ,'kp_6_bands':{
      # gamma1 :5.3    gamma2 :1.7    gamma3 :2.0           # ???
        'L' :-13.05      ,'M' : -3.03      ,'N' :-11.88             # ???
        ,'kappa' :1.27                                           # ???
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :3.49
        'S' :0.82                                               # ??? S :1 + 2F
        ,'E_P' :18.8                                             # ???
        ,'B' : 0                     
      # gamma1 :1.47   gamma2 :-0.28  gamma3 :0.03          # ???
        ,'L' :-13.05      ,'M' : -3.03      ,'N' :-11.88             # ???
        ,'kappa' :-1.31                                          # ???
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }


}
#

######### zinc sulfide ############################################
,'ZnS':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :5.4102                                   # [Angstrom]    300 K, J. Piprek
        ,'a_expansion' :3.68e-5                                  # [Angstrom/K]         J. Piprek
    }
    
    ,'dielectric_consts':{
        'static_a'  : 8.9                                       # R. T. Senger et al., phys. stat. sol. (b) 241, 1896 (2004)
        ,'optical_a' : 5.10                                      # R. T. Senger et al., phys. stat. sol. (b) 241, 1896 (2004)
    }

    ,'elastic_consts':{
        'c11' :106.7  ,'c12' :66.6   ,'c44' :45.6                   # G. Martino et al., phys. stat. sol. (,'a') 152, 249 (1995)
    }                    

    ,'piezoelectric_consts':{
        'e14' :0.0                                              # ???
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.20                         # R. T. Senger et al., phys. stat. sol. (b) 241, 1896 (2004)
#         ,'bandgap'               :3.68                         # 300 K
          ,'bandgap'               :3.841                        #   0 K
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-6.9                         # A. Zunger: a_c :a_v + a_gap :-1.74 + (-5.16) :-6.9, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.83                        # A. Zunger: a_c :a_v + a_gap :-1.74 + (-3.09) :-4.83, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.65                        # A. Zunger: a_c :a_v + a_gap :-1.74 + 1.09 :-0.65, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :-0.02                              # A. Zunger, average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :1.76  ,'g' :0 }                      # J. Piprek
        ,'LH':{ 'mass'          :0.17  ,'g' :0 }                      # J. Piprek
        ,'SO':{ 'mass'          :0.40  }                            # P. Lawaetz, PRB 4, 3460 (1971)

        ,'defpot_absolute'   :-1.74                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-0.8   ,'defpot_uniaxial_d' :-3.7    # b (J. Piprek), d (experiment, A. Blacha et al., phys. stat. sol. (b) 126, 11 (1984)) [eV]
       
        ,'delta_SO'          :0.07                               # Qteish/Needs
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.8    gamma2 :0.67   gamma3 :1.53          # ???
        'L' :-9.86       ,'M' :-3.02       ,'N' :-11.04             # ???
        ,'kappa' :0.17                                           # P. Lawaetz, PRB 4, 3460 (1971)
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :3.841
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :20.4                                             # P. Lawaetz, PRB 4, 3460 (1971)
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-3.02       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#

######### cadmium sulfide ############################################
,'CdS':{
    'mat_crys_strc'    :'Zincblende'
    ,'valence' :'II_VI'
   
    ,'lattice_consts':{
        'a'           :5.818                                    # [Angstrom]    300 K, J. Piprek
        ,'a_expansion' :3.68e-5                                  # [Angstrom/K] ??? ZnS
    }
    
    ,'dielectric_consts':{
        'static_a'  : 8.9                                       # ???
        ,'optical_a' : 5.10                                      # ???
    }

    ,'elastic_consts':{
        'c11' :83.1  ,'c12' :50.4   ,'c44' :45.6                   # J. Piprek, (,'c44' ??? ZnS value)
    }                    

    ,'piezoelectric_consts':{
        'e14' :0.0                                              # ???
    }                                                

    ,'conduction_bands':{
       'Gamma':{ 
          'mass'                  :0.21                         # J. Piprek
#         ,'bandgap'               :2.48                         # 300 K
          ,'bandgap'               :2.583                        #   0 K
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-4.45                         # A. Zunger: a_c :a_v + a_gap :-1.51 + (-2.94) :-4.45, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'g'                     :0                            # ???
       }
       ,'L':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-3.74                        # A. Zunger: a_c :a_v + a_gap :-1.51 + (-2.23) :-3.74, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }
       ,'X':{ 
          'mass_l'                :0.5                          # ???
          ,'mass_t'                :0.5                          # ???
          ,'bandgap'               :5                            # ???
          ,'bandgap_alpha'         :0                            # ???
          ,'bandgap_beta'          :0                            # ???
          ,'defpot_absolute'       :-0.63                        # A. Zunger: a_c :a_v + a_gap :-1.51 + 0.88 :-0.63, S.-H. Wei and A. Zunger, APL 72, 2011 (1998)
          ,'defpot_uniaxial'       :0                            # ???
       }      
    }

    ,'valence_bands':{
        'bandoffset'        :0.16                               # A. Zunger, average ,'valence' band energy E_v,av [eV]
        
        ,'HH':{ 'mass'          :0.64  ,'g' :0 }                      # J. Piprek
        ,'LH':{ 'mass'          :0.64  ,'g' :0 }                      # J. Piprek
        ,'SO':{ 'mass'          :0.40  }                            # ???

        ,'defpot_absolute'   :-1.51                              # A. Zunger: a_v
        ,'defpot_uniaxial_b' :-1.18   ,'defpot_uniaxial_d' :-3.7   # b (Qteish, +1.6 eV (exp.)?), d (??? ZnS) [eV]
       
        ,'delta_SO'          :0.064                              # J. Piprek
       
    }

    ,'kp_6_bands':{
      # gamma1 :4.8    gamma2 :0.67   gamma3 :1.53          # ???
        'L' :-9.86       ,'M' :-3.02       ,'N' :-11.04             # ???
        ,'kappa' :0                                              # ???
    }                    

    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :2.583
        'S' :1.0                                                # S :1 + 2F
        ,'E_P' :20.4                                             # ???
        ,'B' : 0                     
      # gamma1 :0.0    gamma2 :0.0    gamma3 :0.0           # ??? can be calculated from 6x6 k.p parameters
        ,'L' :0           ,'M' :-3.02       ,'N' :0                  # ??? can be calculated from 6x6 k.p parameters
        ,'kappa' :0                                              # ??? can be calculated from 6x6 k.p parameters
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

}
#
    
,'Air_wz':{                                                   
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'III_V'                                            # III_VI
   
    ,'lattice_consts':{
        'a'           :3.189                                    
        ,'a_expansion' :0                                        
        ,'c'           :5.185                                    
        ,'c_expansion' :0                                        
    }
    
    ,'dielectric_consts':{    
        'static_a'  :1.0     ,'static_c'  :1.0                            ,'optical_a' :1.0    ,'optical_c' :1.0                   
    }

    ,'elastic_consts':{
        'c11' :0.01       ,'c12' :0.01       ,'c13' :0.01              
        ,'c33' :0.01       ,'c44' :0.01 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :0.0     ,'e33' :0.0                             
        ,'e15' :0.0                                            
    }  
    
    ,'pyroelectric_const' :0.0                                   
       
    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.202                         
          ,'mass_t'                :0.206                         
          ,'bandgap'               :8.0                        
          ,'bandgap_alpha'         :0.0  
          ,'bandgap_beta'          :0    
          ,'defpot_absolute_l'     :0                         
          ,'defpot_absolute_t'     :0                        
       }
    }

    ,'valence_bands':{
        'bandoffset'        :-3.0                            
        
        ,'HH':{ 'mass_l' :1.6      ,'mass_t' :1.1 }                
        ,'LH':{ 'mass_l' :0.15     ,'mass_t' :0.11  }               
        ,'SO':{ 'mass_l' :1.1      ,'mass_t' :0.15  }               
        ,'defpotentials' :{'D1':-3.90,
                            'D2':-4.13,
                            'D3':1.15,
                            'D4':-1.22,
                            'D5':-1.53,
                            'D6':2.83}                   
       
        ,'delta':{'delta_': 0.010,'delta_so':0.00567,'delta_cr':0.00567}    
                  
    }

    ,'kp_6_bands':{  
        'A1' :-7.21     ,'A2' :-0.44      ,'A3' : 6.68              
        ,'A4' :-3.46     ,'A5' :-3.40      ,'A6' :-4.90              
    }
 
    ,'kp_8_bands':{                                                
        'S1'   : 1.0   ,'S2'   : 1.0                        
        ,'E_P1' :0.0      ,'E_P2' :0.0                           
        ,'B1'   : 0        ,'B2'   : 0         ,'B3' : 0            
        ,'A1'   :-7.21    ,'A2'   :-0.44      ,'A3' : 6.68         
        ,'A4'   :-3.46    ,'A5'   :-3.40     ,'A6' : -4.90         
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }        
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }        
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 
        }
        ,'radiative':{ 'c' :1.1e-8                                  
        }
    }
} 
#


######### sapphire #####################################################
,'Al2O3':{
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'III_V'                                            # III_VI
   
    ,'lattice_consts':{
        'a'           :3.112                                    # AlN value
        ,'a_expansion' :0                                        # ?
        ,'c'           :4.982                                    # AlN value
        ,'c_expansion' :0                                        # ?
    }
    
    ,'dielectric_consts':{    
        'static_a'  :9.4     ,'static_c'  :11.5                   # www.crystran.co.uk/sappdata.htm  11.5 (para) 9.4 (perp) at 1MHz
        ,'optical_a' :4.68    ,'optical_c' :4.68                   # AlN value
    }

    ,'elastic_consts':{
        'c11' :496       ,'c12' :164       ,'c13' :115              # www.crystran.co.uk/sappdata.htm
        ,'c33' :498       ,'c44' :148 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :-0.50     ,'e33' :1.79                             # AlN value
        ,'e15' :-0.48                                            # AlN value
    }  
    
    ,'pyroelectric_const' :0.0                                   #
       
    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.32                         # AlN value
          ,'mass_t'                :0.30                         # AlN value
          ,'bandgap'               :6.25                         # AlN value
          ,'bandgap_alpha'         :0.0                          # ?
          ,'bandgap_beta'          :0                            # ?
          ,'defpot_absolute_l'     :-3.4                         # AlN value  along c axis
          ,'defpot_absolute_t'     :-11.8                        # AlN value  perpendicular to c axis
       }
    }

    ,'valence_bands':{
        'bandoffset'        :-1.526                             # AlN value
        
        ,'HH':{ 'mass_l' :3.53      ,'mass_t' :10.42 }                # AlN value
        ,'LH':{ 'mass_l' :3.53      ,'mass_t' :0.24  }                # AlN value
        ,'SO':{ 'mass_l' :0.25      ,'mass_t' :3.81  }                # AlN value
        
        ,'defpotentials' :{'D1':-3.90,
                            'D2':-4.13,
                            'D3':1.15,
                            'D4':-1.22,
                            'D5':-1.53,
                            'D6':2.83}                   # ?
       
        ,'delta':{'delta_': -0.169,'delta_so':0.00633,'delta_cr':0.00633}# AlN value    

    }

    ,'kp_6_bands':{  
        'A1' :-3.86     ,'A2' :-0.25      ,'A3' : 3.58              # AlN value
        ,'A4' :-1.32     ,'A5' :-1.47      ,'A6' :-1.64              # AlN value
    }
 
    ,'kp_8_bands':{                                                #
        'S1'   : 0.805    ,'S2'   : 1.013                         # AlN value
        ,'E_P1' :14.5      ,'E_P2' :14.5                           # AlN value
        ,'B1'   : 0        ,'B2'   : 0         ,'B3' : 0             # AlN value
        ,'A1'   :-1.540    ,'A2'   :-0.25      ,'A3' : 1.260         # AlN value
        ,'A4'   :-0.160    ,'A5'   :-0.310     ,'A6' : 4.877         # AlN value
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # AlN value
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # AlN value
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :1.1e-8                                  # GaN value !!!
        }
    }
}
#




##########################################################################################
##########################################################################################
#          BINARIES     --    III - V       V A L E N C E
##########################################################################################
##########################################################################################


######### gallium nitride #############################################
,'GaN':{
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :3.189                                    # Vurgaftman1/Vurgaftman2 (300 K) and O. Ambacher, Review
        ,'a_expansion' :0                                        # 
      # ,'a_expansion' :5.59e-5                                  # This value corresponds to ,'a' different equation! http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/thermal.html Qian W. et al., MRS Symposium Proceedings, Pittsburgh, 475-486 (1996)
        ,'c'           :5.185                                    # Vurgaftman1/Vurgaftman2 (300 K) and O. Ambacher, Review
        ,'c_expansion' :0                                        # 
      # ,'c_expansion' :3.17e-5                                  # This value corresponds to ,'a' different equation! http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/thermal.html Qian W. et al., MRS Symposium Proceedings, Pittsburgh, 475-486 (1996)
    }
    
    ,'dielectric_consts':{    
        'static_a'  :9.28    ,'static_c'  :10.10                  # Tsai et al., JAP 85, 1475 (1999)
        ,'optical_a' :5.29    ,'optical_c' :5.29                   # S.M. Komirenko et al., PRB 59, 5013 (1999) (value  taken from paper of V.A. Fonoberov et al., JAP 94, 7178 (2003)
    }

    ,'elastic_consts':{
        'c11' :390       ,'c12' :145       ,'c13' :106              # Vurgaftman1/Vurgaftman2
        ,'c33' :398       ,'c44' :105 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :-0.35     ,'e33' :1.27                             # Vurgaftman1 (Vurgaftman2 lists d_ij (/:e_ij !) parameters.)
        ,'e15' :-0.30                                            # O. Ambacher
    }  
    
    ,'pyroelectric_const' :-0.034                                # Vurgaftman2 and O. Ambacher

    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.206                        # m_perp=0.202, m_perp=0.202, m_par=0.206 (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/bandstr.html
          ,'mass_t'                :0.202                        # m_perp=0.202, m_perp=0.202, m_par=0.206 (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/bandstr.html
          ,'bandgap'               :3.510                        # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.909e-3                     # Vurgaftman2
          ,'bandgap_beta'          :830                          # Vurgaftman2
          #-----------------------------------------------------------------------------------------------
          # Note that I. Vurgaftman et al., JAP 94, 3675 (2003) lists a_1 and a_2 parameters.
          # They refer to the interband deformation potentials, i.e. to the deformation of the band gaps.
          # Thus we have to add the deformation potentials of the ,'valence' bands
          # to get the deformation potentials for the conduction band edge.
          # a_c,,'a' :a_2 + D2 :-11.3 +   4.5  :-6.8    [Vurgaftman2]
          # a_c,c :a_1 + D1 : -4.9 + (-3.7) :-8.6    [Vurgaftman2]
          #-----------------------------------------------------------------------------------------------
          ,'defpot_absolute_l'     :-8.6                         # Vurgaftman2 (a1) along c axis
          ,'defpot_absolute_t'     :-6.8                         # Vurgaftman2 (a2) perpendicular to c axis
       }
    }

    ,'valence_bands':{
        'bandoffset'        :-0.726                             # A. Zunger
        
        ,'HH':{ 'mass_l' :1.1       ,'mass_t' :1.6   }                # m_perp=1.6 , m_perp=1.6 , m_par=1.1  (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/bandstr.html
        ,'LH':{ 'mass_l' :1.1       ,'mass_t' :0.15  }                # m_perp=0.15, m_perp=0.15, m_par=1.1  (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/bandstr.html
        ,'SO':{ 'mass_l' :0.15      ,'mass_t' :1.1   }                # m_perp=1.1 , m_perp=1.1 , m_par=0.15 (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaN/bandstr.html
        
        ,'defpotentials' :{'D1':-3.7,
                            'D2':4.5,
                            'D3':8.2,
                            'D4':-4.1,
                            'D5':-4.0,
                            'D6':-5.5
    }
        
        ,'delta':{'delta_': 0.010,'delta_so':0.00567,'delta_cr':0.00567}    
      # ,'delta' :{[ 0.0108, 0.00703, 0.00703 ]}                   # Ren et al.
    }

    # Note: The GaN values 'A1',,'A2',,'A3',,'A4',,'A5',,'A6' are taken from Vurgaftman2. He took the values of Ren et al.
    #       To be consistent with Ren's values one should use Ren's values for the delta splittings.
    #       Delta_1(,'cr')=0.0108d0 and Delta_2=Delta_3=0.0211/3=Delta_so/3=0.00703d0
    ,'kp_6_bands':{  
        'A1' :-7.21     ,'A2' :-0.44      ,'A3' : 6.68              # Vurgaftman2
        ,'A4' :-3.46     ,'A5' :-3.40      ,'A6' :-4.90              # Vurgaftman2
    }
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :3.510 eV  (Vurgaftman2)
        'S1'   : 0.866    ,'S2'   : 0.962
      # 'S1'   : 1        ,'S2'   : 1                             # S :1 + 2F :1 + 2 (0) :1 (Vurgaftman1)
        ,'E_P1' : 14.0     ,'E_P2' : 14.0                          # Vurgaftman1
        ,'B1'   : 0.0937   ,'B2'   : 0.0937    ,'B3' : 0.0937        # Vurgaftman2 A7 :0.0937
        ,'A1'   :-3.221    ,'A2'   :-0.44      ,'A3' : 2.691         # based on Vurgaftman2
        ,'A4'   :-1.466    ,'A5'   :-1.406     ,'A6' :-2.080         # based on Vurgaftman2
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         #
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         #
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :1.1e-8                                  # radiative constant :1.1d-8 [cm^3/s] - J.F.Muth APL 71 (18) 2572 (1997)
        }
    }
}
#

######### aluminum nitride ############################################
,'AlN':{
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'III_V'   
   
    ,'lattice_consts':{
        'a'           :3.112                                    # Vurgaftman1/Vurgaftman2 (300 K) and O. Ambacher, Review
        ,'a_expansion' :0                                        # 
      # ,'a_expansion' :4.15e-5                                  # This value corresponds to ,'a' different equation! http://www.ioffe.rssi.ru/SVA/NSM/Semicond/AlN/basic.html Sirota, N.N., Golodushko, V.Z., Tezisy Dokl., Vses Konf. Khi., Svyazi Poluprovdn. Polumetallakh 5th (1974) 98
        ,'c'           :4.982                                    # Vurgaftman1/Vurgaftman2 (300 K) and O. Ambacher, Review
        ,'c_expansion' :0                                        # 
      # ,'c_expansion' :5.27e-5                                  # This value corresponds to ,'a' different equation! http://www.ioffe.rssi.ru/SVA/NSM/Semicond/AlN/basic.html Sirota, N.N., Golodushko, V.Z., Tezisy Dokl., Vses Konf. Khi., Svyazi Poluprovdn. Polumetallakh 5th (1974) 98
    }
    
    ,'dielectric_consts':{    
        'static_a'  :8.67    ,'static_c'  :8.57                   # (values taken from paper of V.A. Fonoberov et al., JAP 94, 7178 (2003)
        ,'optical_a' :4.68    ,'optical_c' :4.68                   # S.M. Komirenko et al., PRB 59, 5013 (1999) (value  taken from paper of V.A. Fonoberov et al., JAP 94, 7178 (2003)
    }

    ,'elastic_consts':{
        'c11' :396       ,'c12' :137       ,'c13' :108              # Vurgaftman1/Vurgaftman2
        ,'c33' :373       ,'c44' :116 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :-0.50     ,'e33' :1.79                             # Vurgaftman1 (Vurgaftman2 lists d_ij (/:e_ij !) parameters.)
        ,'e15' :-0.48                                            # O. Ambacher
    }  
    
    ,'pyroelectric_const' :-0.090                                # Vurgaftman2 and O. Ambacher
       
    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.32                         # Vurgaftman2  m_perp=0.30, m_perp=0.30, m_par=0.32 (with respect to c-axis)
          ,'mass_t'                :0.30                         # Vurgaftman2  m_perp=0.30, m_perp=0.30, m_par=0.32 (with respect to c-axis)
          ,'bandgap'               :6.25                         # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :1.799e-3                     # Vurgaftman2
          ,'bandgap_beta'          :1462                         # Vurgaftman2
          #-----------------------------------------------------------------------------------------------
          # Note that I. Vurgaftman et al., JAP 94, 3675 (2003) lists a_1 and a_2 parameters.
          # They refer to the interband deformation potentials, i.e. to the deformation of the band gaps.
          # Thus we have to add the deformation potentials of the ,'valence' bands
          # to get the deformation potentials for the conduction band edge.
          # a_c,,'a' :a_2 + D2 :-11.8 +    7.9  :-3.9    [Vurgaftman2]
          # a_c,c :a_1 + D1 : -3.4 + (-17.1) :-20.5   [Vurgaftman2]
          #-----------------------------------------------------------------------------------------------
          ,'defpot_absolute_l'     :-20.5                        # Vurgaftman2 (a1) along c axis
          ,'defpot_absolute_t'     : -3.9                        # Vurgaftman2 (a2) perpendicular to c axis
       }
    }

    ,'valence_bands':{
        'bandoffset'        :-1.526                             # A. Zunger
        
        ,'HH':{ 'mass_l' :3.53      ,'mass_t' :10.42 }                # m_perp=10.42, m_perp=10.42, m_par=3.53 (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/AlN/bandstr.html
        ,'LH':{ 'mass_l' :3.53      ,'mass_t' :0.24  }                # m_perp=0.24 , m_perp=0.24 , m_par=3.53 (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/AlN/bandstr.html
        ,'SO':{ 'mass_l' :0.25      ,'mass_t' :3.81  }                # m_perp=3.81 , m_perp=3.81 , m_par=0.25 (with respect to c-axis) - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/AlN/bandstr.html
        
        # Vurgaftman2
        ,'defpotentials' :{'D1':-17.1,
                            'D2':7.9,
                            'D3':8.8,
                            'D4':-3.9,
                            'D5':-3.4,
                            'D6':-3.4 }        
        #,'delta' :{[-0.169 , 0.00633, 0.00633 ]}                   # Vurgaftman2
        ,'delta':{'delta_': 0.010,'delta_so':0.00567,'delta_cr':0.00567}
    }

    ,'kp_6_bands':{  
        'A1' :-3.86     ,'A2' :-0.25      ,'A3' : 3.58              # Vurgaftman2
        ,'A4' :-1.32     ,'A5' :-1.47      ,'A6' :-1.64              # Vurgaftman2
    }
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :6.25 eV   (Vurgaftman2)
        'S1'   : 0.805    ,'S2'   : 1.013
      # 'S1'   : 1        ,'S2'   : 1                             # S :1 + 2F :1 + 2 (0) :1 (Vurgaftman1)
        ,'E_P1' :14.5      ,'E_P2' :14.5                           # Vurgaftman1
        ,'B1'   : 0        ,'B2'   : 0         ,'B3' : 0             # Vurgaftman2 A7 :0
        ,'A1'   :-1.540    ,'A2'   :-0.25      ,'A3' : 1.260         # based on Vurgaftman2
        ,'A4'   :-0.160    ,'A5'   :-0.310     ,'A6' : 4.877         # based on Vurgaftman2
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         #
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         #
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :1.1e-8                                  # GaN value !!!
        }
    }
}
#

######### indium nitride ##############################################
,'InN':{
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'III_V'    
   
    ,'lattice_consts':{
        'a'           :3.545                                    # Vurgaftman1/Vurgaftman2 (300 K)
        ,'a_expansion' :0                                        # 
      # ,'a_expansion' :3.8e-5                                   # This value corresponds to ,'a' different equation! http://www.ioffe.rssi.ru/SVA/NSM/Semicond/InN/thermal.html
        ,'c'           :5.703                                    # Vurgaftman1/Vurgaftman2 (300 K)
        ,'c_expansion' :0                                        # 
      # ,'c_expansion' :2.9e-5                                   # This value corresponds to ,'a' different equation! http://www.ioffe.rssi.ru/SVA/NSM/Semicond/InN/thermal.html
    }
    
    ,'dielectric_consts':{    
        'static_a'  :13.1    ,'static_c'  :14.4                   # http://www.ioffe.rssi.ru/SVA/NSM/Semicond/InN/optic.html
        ,'optical_a' :6.7     ,'optical_c' :6.7                    # A. Kasic PRB 65, 115206 (2002)
    }

     ,'elastic_consts':{
        'c11' :223       ,'c12' :115       ,'c13' :92               # Vurgaftman1/Vurgaftman2
        ,'c33' :224       ,'c44' :48 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :-0.57     ,'e33' :0.97                             # Vurgaftman1 (Vurgaftman2 lists d_ij (/:e_ij !) parameters.)
        ,'e15' :-0.48                                            # AlN value
    }  
    
    ,'pyroelectric_const' :-0.042                                # Vurgaftman2 and O. Ambacher
       
    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.07                         # J. Wu, PRB 66, 201403, Vurgaftman1/Vurgaftman2 m_perp=0.07, m_perp=0.07, m_par=0.07 (with respect to c-axis)
          ,'mass_t'                :0.07                         # J. Wu, PRB 66, 201403, Vurgaftman1/Vurgaftman2 m_perp=0.07, m_perp=0.07, m_par=0.07 (with respect to c-axis)
          ,'bandgap'               :0.78                         # Vurgaftman2 (0 K)
          ,'bandgap_alpha'         :0.245e-3                     # Vurgaftman2
          ,'bandgap_beta'          :624                          # Vurgaftman2
          #-----------------------------------------------------------------------------------------------
          # Note that I. Vurgaftman et al., JAP 94, 3675 (2003) lists a_1 and a_2 parameters.
          # They refer to the interband deformation potentials, i.e. to the deformation of the band gaps.
          # Thus we have to add the deformation potentials of the ,'valence' bands
          # to get the deformation potentials for the conduction band edge.
          # a_c,,'a' :a_2 + D2 :-3.5 +   4.5  : 1.0    [Vurgaftman2]
          # a_c,c :a_1 + D1 :-3.5 + (-3.7) :-7.2    [Vurgaftman2]
          #-----------------------------------------------------------------------------------------------
          ,'defpot_absolute_l'     :-7.2                         # Vurgaftman2 (a1) along c axis
          ,'defpot_absolute_t'     : 1.0                         # Vurgaftman2 (a2) perpendicular to c axis
       }
    }

    ,'valence_bands':{
        'bandoffset'        :-0.462                             # A. Zunger
        
        ,'HH':{ 'mass_l' :1.63      ,'mass_t' :1.63  }                # http://www.ioffe.rssi.ru/SVA/NSM/Semicond/InN/bandstr.html, Xu & Ching (1993), Yeo et al. (1998), Pugh et al. (1999)
        ,'LH':{ 'mass_l' :0.27      ,'mass_t' :0.27  }                # http://www.ioffe.rssi.ru/SVA/NSM/Semicond/InN/bandstr.html, Xu & Ching (1993), Yeo et al. (1998), Pugh et al. (1999)
        ,'SO':{ 'mass_l' :0.65      ,'mass_t' :0.65  }                # http://www.ioffe.rssi.ru/SVA/NSM/Semicond/InN/bandstr.html, Xu & Ching (1993), Yeo et al. (1998), Pugh et al. (1999)
        
        # Vurgaftman2
        ,'defpotentials' :{'D1':-3.7,
                            'D2':4.5,
                            'D3':8.2,
                            'D4':-4.1,
                            'D5':-4.0,
                            'D6':-5.5 }         
        #,'delta' :{[ 0.040 , 0.00167, 0.00167 ]}                   # Vurgaftman2
        ,'delta':{'delta_': 0.010,'delta_so':0.00567,'delta_cr':0.00567}
    }

    ,'kp_6_bands':{  
        'A1' :-8.21     ,'A2' :-0.68      ,'A3' : 7.57              # Vurgaftman2
        ,'A4' :-5.23     ,'A5' :-5.11      ,'A6' :-5.96              # Vurgaftman2
    }
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :0.78 eV   (Vurgaftman2)
        'S1'   :-4.432    ,'S2'   :-4.432
      # 'S1'   : 1        ,'S2'   : 1                             # S :1 + 2F :1 + 2 (0) :1 (Vurgaftman1)
        ,'E_P1' :14.6      ,'E_P2' :14.6                           # Vurgaftman1
        ,'B1'   : 0        ,'B2'   : 0         ,'B3' : 0             # Vurgaftman2 A7 :0
        ,'A1'   :10.508    ,'A2'   :-0.68      ,'A3' :-11.148        # based on Vurgaftman2
        ,'A4'   : 4.129    ,'A5'   : 4.249     ,'A6' :  7.276        # based on Vurgaftman2
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         #
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         #
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :1.1e-8                                  # GaN value !!!
        }
    }
}
#

######### zinc oxide ##############################################
# [W.R.L. Lambrecht et al., MRS Internet J. Nitride Semicond. Res. 4S1, G6.8 (1999)]
,'ZnO':{
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'II_VI'   
   
    ,'lattice_consts':{
        'a'           :3.250                                    # (300 K) http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
        ,'a_expansion' :6.51e-5                                  # (300 K) http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
        ,'c'           :5.205                                    # (300 K) http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
        ,'c_expansion' :3.02e-5                                  # (300 K) http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
    }
    
    ,'dielectric_consts':{    
        'static_a'  :7.77    ,'static_c'  :8.91                   # N. Ashkenov et al., JAP 93, 126 (2003)
        ,'optical_a' :3.70    ,'optical_c' :3.75                   # http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
    }

     ,'elastic_consts':{
      # 'c11 :190       ,'c12' :110       ,'c13' :90               # http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
      # ,'c33' :196       ,'c44' :39 
        'c11' :206       ,'c12' :118       ,'c13' :118              # experiment of G. Carlotti et al., J. Phys.: Condens. Matter 7, 9147 (1995) - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'c33' :211       ,'c44' :44 
      # 'c11' :217       ,'c12' :117       ,'c13' :121              # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
      # ,'c33' :225       ,'c44' :50 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :-0.57     ,'e33' :1.34                             # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'e15' :-0.30                                            # GaN value
    }  
    
    ,'pyroelectric_const' :-0.0220                               # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
  # ,'pyroelectric_const' :-0.0340                               # K. Shimada et al., PRB 88, 075203 (2013) (theory)
  # ,'pyroelectric_const' :-0.050                                # S.-H. Park, D. Ahn, APL 87, 253509 (2007)
  # ,'pyroelectric_const' :-0.03220                              # A. Malashevich, D. Vanderbilt, PRB 75, 045106 (2007) (theory, extrapolation)
       
    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.24                         # [Lambrecht] polaronic mass
          ,'mass_t'                :0.28                         # [Lambrecht] polaronic mass
          ,'bandgap'               :3.436                        # (0 K) [Zhang, 4 K]
          ,'bandgap_alpha'         :0.8e-3                       # http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
          ,'bandgap_beta'          :0                            # 
          ,'defpot_absolute_l'     :-2.30                        # assumption of isotropic deformation potential (Janotti et al., PRB 74, 045202 (2006)
          ,'defpot_absolute_t'     :-2.30                        # assumption of isotropic deformation potential (Janotti et al., PRB 74, 045202 (2006)
       }
    }

    ,'valence_bands':{
        'bandoffset'        :-0.03783                           # 
        
        ,'HH':{ 'mass_l' :0.54      ,'mass_t' :2.74  }                # [Lambrecht]
        ,'LH':{ 'mass_l' :0.55      ,'mass_t' :3.03  }                # [Lambrecht]
        ,'SO':{ 'mass_l' :1.12      ,'mass_t' :0.27  }                # [Lambrecht]
        
        # (5 K, Landoldt-Börnstein, II-VI compounds)
        ,'defpotentials' :{'D1':-3.90,
                            'D2':-4.13,
                            'D3':-1.15,
                            'D4':-1.22,
                            'D5':-1.53,
                            'D6':-2.88 }         
        #,'delta' :{[ 0.050 , 0.0016666, 0.0016666 ]}               # Delta_so :0.005 eV (Claus Klingshirn et al., Physik Journal, Vol. 1, 37 (2006))
        ,'delta':{'delta_': 0.010,'delta_so':0.00567,'delta_cr':0.00567}
    }

    ,'kp_6_bands':{  
        'A1' :-6.68036     ,'A2' :-0.45388     ,'A3' : 6.1275       # W.J. Fan, J. of Crystal Growth 287, 28 (2006)
        ,'A4' :-2.70374     ,'A5' :-2.7669      ,'A6' :-4.62566      # W.J. Fan, J. of Crystal Growth 287, 28 (2006)
    }
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :0.78 eV   (Vurgaftman2)
        'S1'   : 1        ,'S2'   : 1                             # ??? S :1 + 2F :1 + 2 (0) :1 (Vurgaftman1)
        ,'E_P1' :14.0      ,'E_P2' :14.0                           # GaN value
        ,'B1'   : 0        ,'B2'   : 0         ,'B3' : 0             # 
        ,'A1' :-6.68036     ,'A2' :-0.45388     ,'A3' : 6.1275       # kp_6_bands value!!!
        ,'A4' :-2.70374     ,'A5' :-2.7669      ,'A6' :-4.62566      # kp_6_bands value!!!
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :1.1e-8                                  # GaN value !!!
        }
    }
}
#

######### magnesium oxide ##############################################
# [W.R.L. Lambrecht et al., MRS Internet J. Nitride Semicond. Res. 4S1, G6.8 (1999)]
,'MgO':{
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'II_VI'   
   
    ,'lattice_consts':{
        'a'           :3.413                                    # (300 K) W.R.L. Lambrecht et al., MRS Internet J. Nitride Semicond. Res. 4S1, G6.8 (1999)
        ,'a_expansion' :11.2e-5                                  # (300 K) rocksalt value !!! units ? Angstrom? [nm/K] ,'a',,'a',c at 300 K  -  Thermal expansion coefficient: 11.2 x 10-6/°C [http://www.mtberlin.com/frames_cryst/descriptions/substrates.htm]
        ,'c'           :4.109252                                 # (300 K) W.R.L. Lambrecht et al., MRS Internet J. Nitride Semicond. Res. 4S1, G6.8 (1999)
        ,'c_expansion' :11.2e-5                                  # (300 K) http://www.onr.navy.mil/sci_tech/31/312/ncsr/materials/zno.asp
    }
    
    ,'dielectric_consts':{    
        'static_a'  :7.77    ,'static_c'  :8.91                   # ZnO value
        ,'optical_a' :2.95    ,'optical_c' :2.95                   # rocksalt - A.R. Oganov et al., J. of Chem. Phys. 118, 10174 (2003)
    }

     ,'elastic_consts':{
        'c11' :222       ,'c12' :90        ,'c13' :58               # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'c33' :109       ,'c44' :105 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :-0.58     ,'e33' :1.64                             # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'e15' :-0.30                                            # GaN value
    }  
    
    ,'pyroelectric_const' :-0.060                                # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
  # ,'pyroelectric_const' :-0.135                                # theory - K. Shimada et al., PRB 88, 075203 (2013)
  # ,'pyroelectric_const' :-0.111                                # theory - A. Malashevich, D. Vanderbilt, PRB 75, 045106 (2007) (extrapolation, assuming no bowing factor between ZnO and MgO)
  # ,'pyroelectric_const' :-0.070                                # S.-H. Park, D. Ahn, APL 87, 253509 (2007)

    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.24                         # ZnO value
          ,'mass_t'                :0.28                         # ZnO value
          ,'bandgap'               :5.076                        # (0 K) [Zhang]
          ,'bandgap_alpha'         :0.8e-3                       # ZnO value
          ,'bandgap_beta'          :0                            # ZnO value
          ,'defpot_absolute_l'     :-2.30                        # ZnO value
          ,'defpot_absolute_t'     :-2.30                        # ZnO value
       }
    }

    ,'valence_bands':{
        'bandoffset'        :0.61183                            # VBO :-0.574 [Zhang], shift so that VBO(ZnO/MgO) :0.364 eV  (35 % of band offset)
                                                               # E_v,av(MgO) : 1.6793     VBO :-0.03783 - (-1.75783) :1.72 eV
                                                               # 'bandoffset' MgO/ZnO: CBO(MgO-ZnO) :2.68 eV [Kilic and Zunger, Appl. Phys. Lett. 81, 73 (2002)]
                                                               #                     VBO(MgO-ZnO) :1.72 eV [Kilic and Zunger, Appl. Phys. Lett. 81, 73 (2002)]
        
        ,'HH':{ 'mass_l' :2.77      ,'mass_t' :1.6   }                # Y.-N. Xu et al., PRB 43, 4461 (1991)
        ,'LH':{ 'mass_l' :0.31      ,'mass_t' :0.32  }                # Y.-N. Xu et al., PRB 43, 4461 (1991)
        ,'SO':{ 'mass_l' :1.12      ,'mass_t' :0.27  }                # ZnO value
        
        # ZnO values
        ,'defpotentials' :{'D1':-3.90,
                            'D2':-4.13,
                            'D3':-1.15,
                            'D4':-1.22,
                            'D5':-1.53,
                            'D6':-2.88 }         
        #,'delta' :{[ 0.050 , 0.0016666, 0.0016666 ]}               # ZnO values
        ,'delta':{'delta_': 0.010,'delta_so':0.00567,'delta_cr':0.00567}
    }

    ,'kp_6_bands':{  
        'A1' :-6.68036     ,'A2' :-0.45388     ,'A3' : 6.1275       # ZnO values
        ,'A4' :-2.70374     ,'A5' :-2.7669      ,'A6' :-4.62566      # ZnO values
    }
 
    ,'kp_8_bands':{                                                # ,'bandgap'(Gamma) :0.78 eV   (Vurgaftman2)
        'S1'   : 1        ,'S2'   : 1                             # ??? S :1 + 2F :1 + 2 (0) :1 (Vurgaftman1)
        ,'E_P1' :14.0      ,'E_P2' :14.0                           # GaN value
        ,'B1'   : 0        ,'B2'   : 0         ,'B3' : 0             # 
        ,'A1' :-6.68036     ,'A2' :-0.45388     ,'A3' : 6.1275       # kp_6_bands value!!!
        ,'A4' :-2.70374     ,'A5' :-2.7669      ,'A6' :-4.62566      # kp_6_bands value!!!
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :1.1e-8                                  # GaN value !!!
        }
    }
}
#

######### cadmium oxide ##############################################
# Note: Usually CdO has the rocksalt crystal structure:
#       (cubic Fm3m, rocksalt ,'B1')
#       However, it is possible to grow Cd(x)Zn(1-x)O alloys with the wurtzite structure.
,'CdO':{
    'mat_crys_strc'    :'Wurtzite'
    ,'valence' :'II_VI'   
   
    ,'lattice_consts':{
        'a'           :3.60                                     # (300 K) P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'a_expansion' :6.51e-5                                  # ZnO value
        ,'c'           :5.58                                     # (300 K) P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'c_expansion' :3.02e-5                                  # ZnO value
    }

    ,'dielectric_consts':{    
        'static_a'  :9.65    ,'static_c'  :9.65                   # MgO value
        ,'optical_a' :2.95    ,'optical_c' :2.95                   # MgO value
    }

     ,'elastic_consts':{
        'c11' :150       ,'c12' :108       ,'c13' :61               # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'c33' :105       ,'c44' :47 
    }
                                             
    ,'piezoelectric_consts':{ 
        'e31' :-0.48     ,'e33' :1.67                             # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006)
        ,'e15' :-0.30                                            # GaN value
    }  
    
  # ,'pyroelectric_const' :-0.106                               # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006) (Table II)
    ,'pyroelectric_const' :-0.10                                # theory - P. Gopal, N.A. Spaldin, Journal of Electronic Materials 35, 538 (2006) (Abstract)
       
    ,'conduction_bands':{
        'Gamma':{
          'mass_l'                :0.24                         # ZnO value
          ,'mass_t'                :0.28                         # ZnO value
          ,'bandgap'               :3.436                        # ZnO value
          ,'bandgap_alpha'         :0.8e-3                       # ZnO value
          ,'bandgap_beta'          :0                            # 
          ,'defpot_absolute_l'     :-2.30                        # ZnO value
          ,'defpot_absolute_t'     :-2.30                        # ZnO value
       }
    }

    ,'valence_bands':{
        'bandoffset'        :-0.03783                           # ZnO value
        
        ,'HH':{ 'mass_l' :0.54      ,'mass_t' :2.74  }                # ZnO value
        ,'LH':{ 'mass_l' :0.55      ,'mass_t' :3.03  }                # ZnO value
        ,'SO':{ 'mass_l' :1.12      ,'mass_t' :0.27  }                # ZnO value
        
   
        ,'defpotentials' :{'D1':-3.90,# ZnO value
                            'D2':-4.13,
                            'D3':-1.15,
                            'D4':-1.22,
                            'D5':-1.53,
                            'D6':-2.88 }        
        #,'delta' :{[ 0.050 , 0.0016666, 0.0016666 ]}               # ZnO value
        ,'delta':{'delta_': 0.010,'delta_so':0.00567,'delta_cr':0.00567}
    }

    ,'kp_6_bands':{  
        'A1' :-6.68036     ,'A2' :-0.45388     ,'A3' : 6.1275       # ZnO value
        ,'A4' :-2.70374     ,'A5' :-2.7669      ,'A6' :-4.62566      # ZnO value
    }
 
    ,'kp_8_bands':{                                                # ZnO value
        'S1'   : 1        ,'S2'   : 1                             # ZnO value
        ,'E_P1' :14.0      ,'E_P2' :14.0                           # ZnO value
        ,'B1'   : 0        ,'B2'   : 0         ,'B3' : 0             # 
        ,'A1' :-6.68036     ,'A2' :-0.45388     ,'A3' : 6.1275       # ZnO value
        ,'A4' :-2.70374     ,'A5' :-2.7669      ,'A6' :-4.62566      # ZnO value
    }

    ,'mobility_constant':{
        'electrons':{  'mumax' :100      ,'exponent' :1.0  }         # ???
        ,'holes':{      'mumax' :100      ,'exponent' :1.0  }         # ???
    }

    ,'recombination':{  
        'SRH':{       'tau_n' :1.0e-9      ,'nref_n' :1.0e19         # InP value !!!
                   ,'tau_p' :1.0e-9      ,'nref_p' :1.0e18         # InP value !!!
        }
        ,'Auger':{     'c_n' :0             ,'c_p' :0                 # InP value !!!
        }
        ,'radiative':{ 'c' :1.1e-8                                  # GaN value !!!
        }
    }
}
#    
}

alloyproperty ={

        ##########################################################################################
##########################################################################################
#          T E R N A R Y     A L L O Y S     --    IV - IV       V A L E N C E
##########################################################################################
##########################################################################################




######### silicon-germanium (SiGe) ####################################
'SiGe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'IV_IV'
    ,'binary_x'   :'Ge'
    ,'binary_1_x' :'Si'    

    ,'conduction_bands':{
        'Delta':{
          'bandgap' :0.206                                      #   4.2 K, J. Weber et al., PRB 40, 5683 (1989)
        }
    }

    ,'kp_6_bands':{  
        'L' :0   ,'M' :0   ,'N' :0              # !!! L,N !!!       # M varies linearly with x - M. Rieger, P. Vogl, PRB 48, 14276 (1993) - for L,N see eq. (35) in this paper
    }                                   

}
#
##########################################################################################
##########################################################################################
#          T E R N A R Y     A L L O Y S     --    III - V       V A L E N C E
##########################################################################################
##########################################################################################




######### aluminum gallium arsenide (AlGaAs) ##########################
#                                                ,'g' :0.4       # (g factor of Al0.3Ga0.7As)
,'AlGaAs_Bowing_x':{
    'mat_crys_strc'       :'Zincblende'    
    ,'valence'    :'III_V'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :-0.127 + 1.310 * 1 }                  # Vurgaftman1: -0.127 + 1.310 * x
        ,'L':    { 'bandgap' : 0     }                              # Vurgaftman1
        ,'X':    { 'bandgap' : 0.055 }                              # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO' :0                                           # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :14000 }    	                       # http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/hall.html
        ,'holes':{ 'mumax' : 1000 }	                           # geschätzt S. Ziegler (E26)
    }
}
#


,'AlGaAs_Bowing_1_x':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :-0.127 + 1.310 * 0 }                  # Vurgaftman1: -0.127 + 1.310 * x
        ,'L':    { 'bandgap' : 0     }                              # Vurgaftman1
        ,'X':    { 'bandgap' : 0.055 }                              # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO' :0                                           # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :-1000 }			                   # http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/hall.html
        ,'holes':{ 'mumax' : 1000 }	    		               # geschätzt S. Ziegler (E26)
    }
}
#


,'AlGaAs':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'

    ,'binary_x'   :'AlAs'
    ,'binary_1_x' :'GaAs' 
    ,'bowing_x'   :'AlGaAs_Bowing_x'
    ,'bowing_1_x' :'AlGaAs_Bowing_1_x'
}
#

######### indium gallium arsenide (InGaAs) ############################
,'InGaAs':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InAs'
    ,'binary_1_x' :'GaAs'

    ,'conduction_bands':{
        'Gamma':{ 'mass'            :0.0091                        # Vurgaftman1
               ,'bandgap'         :0.477                         # Vurgaftman1
               ,'defpot_absolute' :2.61 }                        # Vurgaftman1
        ,'L':    { 'bandgap'         :0.33 }                        # Vurgaftman1
        ,'X':    { 'bandgap'         :1.4  }                        # Vurgaftman1
    }

    ,'valence_bands':{
        'bandoffset' :-0.38                                     # Vurgaftman1 (Does this value refer to the ,'valence' band edge or to the average ,'valence' band edge?)
        ,'HH':{ 'mass'   :-0.145  }                                 # Vurgaftman1: hh along [001]
        ,'LH':{ 'mass'   : 0.0202 }                                 # Vurgaftman1: lh along [001]
        ,'delta_SO'   : 0.15                                     # Vurgaftman1
    }
                                                               # gamma1 (In0.53Ga0.47As) :11.01  ==>  ,'L' :-28.73
                                                               # gamma2 (In0.53Ga0.47As) : 4.18  ==>  ,'M' :-3.65
                                                               # gamma3 (In0.53Ga0.47As) : 4.84  ==>  ,'N' :-29.04
    ,'kp_6_bands':{                                                # ? gamma3 - gamma2 :0.481 (Vurgaftman1)
      # gamma1 :   gamma2 :   gamma3 : 
        'L' :-32.28984344                                       # => ,'L' :-28.73 (for In0.53Ga0.47As)
        ,'M' :-1.140907266                                       # => ,'M' :-3.65  (for In0.53Ga0.47As)
        ,'N' :-34.03693296                                       # => ,'N' :-29.04 (for In0.53Ga0.47As)
        ,'kappa' :4.96                                           # approximation ,'kappa' :-N/6+M/3-1/3
      # ,'kappa' :4.0                                            # Traynor, PRB 51, 7361 (1995) <- uses kappa=1.1 for GaAs 
    }
    ,'kp_8_bands':{
        'S' :3.54                                               # S :2 * F :2 * 1.77 (Vurgaftman1)
        ,'E_P' :-1.48                                            # Vurgaftman1
        ,'L' :-25.063                                            # consistent to 6-band (for In0.53Ga0.47As)
        ,'M' :-1.141                                             # consistent to 6-band (for In0.53Ga0.47As)
        ,'N' :-26.809                                            # consistent to 6-band (for In0.53Ga0.47As)
        ,'kappa' :3.75                                           # approximation ,'kappa' :-N/6+M/3-1/3
      # ,'kappa' :4.47                                           # consistent to 6-band (Traynor)
                                                               # ? gamma3 - gamma2 :0.481 (Vurgaftman1)                         
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :41000 }			                   # http://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/hall.html
        ,'holes':{ 'mumax' :    0 }			                   # ?
    }
}
#

######### indium aluminium arsenide (InAlAs) ##########################
,'InAlAs':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InAs'
    ,'binary_1_x' :'AlAs' 

    ,'conduction_bands':{
        'Gamma':{ 'mass'            : 0.049                        # Vurgaftman1
               ,'bandgap'         : 0.70                         # Vurgaftman1
               ,'defpot_absolute' :-1.4 }                        # Vurgaftman1
        ,'X':    { 'bandgap'         : 0   }                        # Vurgaftman1
    }

    ,'valence_bands':{
        'bandoffset' :-0.64                                     # Vurgaftman1 (Does this value refer to the ,'valence' band edge or to the average ,'valence' band edge?)
        ,'delta_SO'   : 0.15                                     # Vurgaftman1
    }

    ,'kp_8_bands':{
        'S'   :-8.88                                            # S :2 * F :2 * (-4.44) (Vurgaftman1)
        ,'E_P' :-4.81                                            # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' : 0 }			                       # 
        ,'holes':{ 'mumax' : 0 }			                       #  
    }
}
#

######### gallium indium phosphide (GaInP) ############################
,'GaInP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'GaP'
    ,'binary_1_x' :'InP'

    ,'conduction_bands':{
        'Gamma':{ 'mass'       : 0.051                             # Vurgaftman1
               ,'bandgap'    : 0.65 }                            # Vurgaftman1
        ,'L':    { 'bandgap'    : 1.03 }                            # Vurgaftman1
        ,'X':    { 'bandgap'    : 0.20 }                            # Vurgaftman1
    }

    ,'valence_bands':{
        'defpot_uniaxial_d' : 0                                 # Vurgaftman1
        ,'delta_SO'          : 0                                 # Vurgaftman1
    }

    ,'kp_8_bands':{
        'S' :1.56                                               # S :2 * F :2 * 0.78 (Vurgaftman1)
    }
}
#

######### aluminium indium phosphide (AlInP) ##########################
,'AlInP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlP'
    ,'binary_1_x' :'InP'

    ,'conduction_bands':{
        'Gamma':{ 'mass'    : 0.22                                 # Vurgaftman1
               ,'bandgap' :-0.48 }                               # Vurgaftman1
        ,'X':    { 'bandgap' : 0.38 }                               # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO'       :-0.19                                 # Vurgaftman1
    }
}
#

######### aluminium gallium phosphide (AlGaP) #########################
,'AlGaP':{
    'mat_crys_strc'       :'Zincblende' 
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlP'
    ,'binary_1_x' :'GaP'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' : 0    }                               # Vurgaftman1
        ,'X':    { 'bandgap' : 0.13 }                               # Vurgaftman1
    }
}
#

######### gallium indium antimonide (GaInSb) ##########################
,'GaInSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'GaSb'
    ,'binary_1_x' :'InSb'

    ,'conduction_bands':{
        'Gamma':{ 'mass'    :0.0092                                # Vurgaftman1
               ,'bandgap' :0.415 }                               # Vurgaftman1
        ,'L':    { 'bandgap' :0.4   }                               # Vurgaftman1
        ,'X':    { 'bandgap' :0.33  }                               # Vurgaftman1
    }

    ,'valence_bands':{
        'LH':{ 'mass'       :0.011 }                               # Vurgaftman1: lh along [001]
        ,'delta_SO'       :0.1                                   # Vurgaftman1
    }

    ,'kp_8_bands':{
        'S' :-13.68                                             # S :2 * F :2 * (-6.84) (Vurgaftman1)
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' : 60000 }			                   # http://www.ioffe.ru/SVA/NSM/Semicond/GaInSb/electric.html
        ,'holes':{ 'mumax' : 0 }			                       # 
    }
}
#

######### aluminum indium antimonide (AlInSb) #########################
,'AlInSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlSb'
    ,'binary_1_x' :'InSb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.43 }                                # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO' :0.25                                        # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' : 0 }			                       #  
        ,'holes':{ 'mumax' : 0 }			                       #  
    }
}
#

######### aluminum gallium antimonide (AlGaSb) ########################

,'AlGaSb_Bowing_x':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :-0.044 + 1.22 * 1 }                   # Vurgaftman1: -0.044 + 1.22 * x
        ,'L':    { 'bandgap' :0    }                                # Vurgaftman1
        ,'X':    { 'bandgap' :0    }                                # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO' :0.3                                         # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :0 }			                       #  
        ,'holes':{ 'mumax' : 1200 }			                   # A.H. Ramelan & E.M. Goldys - Hole mobility in Al(x)Ga(1-x)Sb grown by metalorganic chemical vapor deposition
    }
}
#


,'AlGaSb_Bowing_1_x':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :-0.044 + 1.22 * 0 }                   # Vurgaftman1
        ,'L':    { 'bandgap' : 0    }                               # Vurgaftman1
        ,'X':    { 'bandgap' : 0    }                               # Vurgaftman1
    }
    ,'valence_bands':{
        'delta_SO' :0.3                                         # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :0 }			                       #  
        ,'holes':{ 'mumax' : 2000 }			                   # A.H. Ramelan & E.M. Goldys - Hole mobility in AlxGa(1-x)Sb grown by metalorganic chemical vapor deposition
    }
}
#


,'AlGaSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlSb'
    ,'binary_1_x' :'GaSb'
    ,'bowing_x'   :'AlGaSb_Bowing_x'
    ,'bowing_1_x' :'AlGaSb_Bowing_1_x'
}
#

######### gallium arsenide antimonide (GaAsSb) ########################
,'GaAsSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'GaSb'
    ,'binary_1_x' :'GaAs'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :1.43 }                                # Vurgaftman1
        ,'L':    { 'bandgap' :1.2  }                                # Vurgaftman1
        ,'X':    { 'bandgap' :1.2  }                                # Vurgaftman1
    }

    ,'valence_bands':{
        'bandoffset' :-1.06                                     # Vurgaftman1 (Does this value refer to the ,'valence' band edge or to the average ,'valence' band edge?)
        ,'delta_SO'   : 0.6                                      # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :0 }			                       #  
        ,'holes':{ 'mumax' : 0 }			                       #  
    }
}
#

######### indium arsenide antimonide (InAsSb) #########################
,'InAsSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InAs'
    ,'binary_1_x' :'InSb'

    ,'conduction_bands':{
        'Gamma':{ 'mass'    :0.035                                 # Vurgaftman1
               ,'bandgap' :0.67 }                                # Vurgaftman1
        ,'L':    { 'bandgap' :0.6  }                                # Vurgaftman1
        ,'X':    { 'bandgap' :0.6  }                                # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO' :1.2                                         # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :-80000 }			                   # http://www.ioffe.ru/SVA/NSM/Semicond/InAsSb/electric.html
        ,'holes':{ 'mumax' : 0 }			                       # 
    }
}
#

######### aluminum arsenide antimonide (AlAsSb) #######################
,'AlAsSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlAs'
    ,'binary_1_x' :'AlSb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.8  }                                # Vurgaftman1
        ,'L':    { 'bandgap' :0.28 }                                # Vurgaftman1
        ,'X':    { 'bandgap' :0.28 }                                # Vurgaftman1
    }

    ,'valence_bands':{
        'bandoffset' :-1.71                                     # Vurgaftman1 (Does this value refer to the ,'valence' band edge or to the average ,'valence' band edge?)
        ,'delta_SO'   : 0.15                                     # Vurgaftman1
    }
    
    ,'mobility_constant':{
        'electrons':{ 'mumax' :0 }			                       # 
        ,'holes':{ 'mumax' : 0 }			                       # 
    }
}
#

######### gallium arsenide phosphide (GaAsP) ##########################
,'GaAsP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'GaP'
    ,'binary_1_x' :'GaAs'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.19 }                                # Vurgaftman1
        ,'L':    { 'bandgap' :0.16 }                                # Vurgaftman1
        ,'X':    { 'bandgap' :0.24 }                                # Vurgaftman1
    }
}
#

######### indium arsenide phosphide (InAsP) ###########################
,'InAsP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InAs'
    ,'binary_1_x' :'InP'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.10 }                                # Vurgaftman1
        ,'L':    { 'bandgap' :0.27 }                                # Vurgaftman1
        ,'X':    { 'bandgap' :0.27 }                                # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO' :0.16                                        # Vurgaftman1
    }
}
#

######### aluminium arsenide phosphide (AlAsP) ########################
,'AlAsP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlAs'
    ,'binary_1_x' :'AlP'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.22 }                                # Vurgaftman1
        ,'L':    { 'bandgap' :0.22 }                                # Vurgaftman1
        ,'X':    { 'bandgap' :0.22 }                                # Vurgaftman1
    }
}
#

######### gallium phosphide antimonide (GaPSb) ########################
,'GaPSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'GaP'
    ,'binary_1_x' :'GaSb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :2.7 }                                 # Vurgaftman1
        ,'L':    { 'bandgap' :1.7 }                                 # Note: Error in Vurgaftman1, see note in Ref. 10 of Vurgaftman2.
        ,'X':    { 'bandgap' :1.7 }                                 # Note: Error in Vurgaftman1, see note in Ref. 10 of Vurgaftman2.
    }
}
#

######### indium phosphide antimonide (InPSb) #########################
,'InPSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InP'
    ,'binary_1_x' :'InSb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :1.9 }                                 # Vurgaftman1
        ,'L':    { 'bandgap' :1.9 }                                 # Vurgaftman1
        ,'X':    { 'bandgap' :1.9 }                                 # Vurgaftman1
    }

    ,'valence_bands':{
        'delta_SO' :0.75                                        # Vurgaftman1
    }
}
#

######### aluminium phosphide antimonide (AlPSb) ######################
,'AlPSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlP'
    ,'binary_1_x' :'AlSb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :2.7 }                                 # Vurgaftman1
        ,'L':    { 'bandgap' :2.7 }                                 # Vurgaftman1
        ,'X':    { 'bandgap' :2.7 }                                 # Vurgaftman1
    }
}
#

######### indium gallium nitride (InGaN) (zincblende) #################
,'InGaN_zb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InN_zb'
    ,'binary_1_x' :'GaN_zb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :1.4 }                                 # Vurgaftman2
        ,'L':    { 'bandgap' :1.84 }                                # Vurgaftman2
        ,'X':    { 'bandgap' :0.69 }                                # Vurgaftman2
    }
}
#

######### aluminum gallium nitride (AlGaN) (zincblende) ###############
,'AlGaN_zb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlN_zb'
    ,'binary_1_x' :'GaN_zb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.7 }                                 # Vurgaftman2
        ,'L':    { 'bandgap' :0.80 }                                # Vurgaftman2
        ,'X':    { 'bandgap' :0.61 }                                # Vurgaftman2
    }
}
#

######### aluminum indium nitride (AlInN) (zincblende) ################
,'AlInN_zb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlN_zb'
    ,'binary_1_x' :'InN_zb'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :2.5 }                                 # Vurgaftman2
        ,'L':    { 'bandgap' :0.80 }                                # Vurgaftman2
        ,'X':    { 'bandgap' :0.61 }                                # Vurgaftman2
    }
}
#

######### gallium arsenide nitride (GaAsN) ############################

,'GaAsN_Bowing_x':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'

    ,'conduction_bands':{                                          # The bowing factor for zinc blende GaAsN reads 20.4-100x [Vurgaftman2] rather than 120.4-100x [Vurgaftman1] (Table XXX).
        'Gamma':{ 'bandgap' :20.4 - 100 * 1 }                      # Vurgaftman2: 20.4 - 100 * x  ==> only valid for x < 0.15
    }

    ,'valence_bands':{
        'delta_SO' :0                                           # Vurgaftman1
    }
}
#


,'GaAsN_Bowing_1_x':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'

    ,'conduction_bands':{                                          # The bowing factor for zinc blende GaAsN reads 20.4-100x [Vurgaftman2] rather than 120.4-100x [Vurgaftman1] (Table XXX).
        'Gamma':{ 'bandgap' :20.4 - 100 * 0 }                      # Vurgaftman2: 20.4 - 100 * x  ==> only valid for x < 0.15
    }

    ,'valence_bands':{
        'delta_SO' :0                                           # Vurgaftman1
    }
}
#


,'GaAsN':{
    'mat_crys_strc'       :'Zincblende'

    ,'valence'    :'III_V'
    ,'binary_x'   :'GaN_zb'
    ,'binary_1_x' :'GaAs'
    ,'bowing_x'   :'GaAsN_Bowing_x'
    ,'bowing_1_x' :'GaAsN_Bowing_1_x'
}
#

######### indium arsenide nitride (InAsN) #############################
,'InAsN':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InN_zb'
    ,'binary_1_x' :'InAs'
}
#

######### indium phosphide nitride (InPN) #############################
,'InPN':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InN_zb'
    ,'binary_1_x' :'InP'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :15 }                                # Vurgaftman1
    }
}
#

######### gallium phosphide nitride (GaPN) ############################
,'GaPN':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'GaN_zb'
    ,'binary_1_x' :'GaP'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :3.9 }                               # Vurgaftman1
    }
}
#

######### indium antimonide nitride (InSbN) ###########################
,'InSbN':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InN_zb'
    ,'binary_1_x' :'InSb'
}
#



######### zinc beryllium selenide (ZnBeSe) ##########################
,'ZnBeSe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'BeSe'
    ,'binary_1_x' :'ZnSe'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.679 }                               # 0.679 corresponds to b=0.97 in Egap and 70 % CBO
        ,'L':    { 'bandgap' :0.679 }                               # 0.679 corresponds to b=0.97 in Egap and 70 % CBO
        ,'X':    { 'bandgap' :0.679 }                               # 0.679 corresponds to b=0.97 in Egap and 70 % CBO
    }

    ,'valence_bands':{
        'bandoffset'  :-0.291                                   # -0.291 corresponds to b=0.97 in Egap and 70 % CBO
    }
}
#

######### zinc sulfide selenide (ZnSSe) ##########################
,'ZnSSe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'ZnS'
    ,'binary_1_x' :'ZnSe'
}
#

######### zinc magnesium selenide (ZnMgSe) ##########################
,'ZnMgSe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'MgSe'
    ,'binary_1_x' :'ZnSe'

    ,'lattice_consts':{
        'a'           :-0.7                                   # [Angstrom]    B. Jobst et al., APL 69, 97 (1996)
    }
    
    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.4 }                               # B. Jobst et al., APL 69, 97 (1996)
    }
}
#

######### cadmium zinc selenide (CdZnSe) ##########################
,'CdZnSe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'CdSe'
    ,'binary_1_x' :'ZnSe'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.245 }                            # ???
        ,'L':    { 'bandgap' :0.245 }                            # ???
        ,'X':    { 'bandgap' :0.245 }                            # ???
    }

    ,'valence_bands':{
        'bandoffset'  :-0.105                                # ???
    }

}
#

######### cadmium zinc telluride (CdZnTe) ##########################
,'CdZnTe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'CdTe'
    ,'binary_1_x' :'ZnTe'
}
#

######### mercury cadmium telluride (HgCdTe) ##########################
,'HgCdTe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'CdTe'
    ,'binary_1_x' :'HgTe'
}
#

######### cadmium magnesium telluride (CdMgTe) ##########################
,'CdMgTe':{
    'mat_crys_strc'       :'Cd(1-x)Mg(x)Te'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'MgTe'
    ,'binary_1_x' :'CdTe'
}
#

######### cadmium manganese telluride (CdMnTe) ##########################
,'CdMnTe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'MnTe'
    ,'binary_1_x' :'CdTe'
}
#

######### cadmium manganese selenide (CdMnSe) ##########################
,'CdMnSe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'MnSe'
    ,'binary_1_x' :'CdSe'
}
#
######### zinc manganese telluride (ZnMnTe) ##########################
,'ZnMnTe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'MnTe'
    ,'binary_1_x' :'ZnTe'
}
#

######### zinc manganese selenide (ZnMnSe) ##########################
,'ZnMnSe':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'MnSe'
    ,'binary_1_x' :'ZnSe'
}
#

######### zinc cadmium sulfide (ZnCdS) ##########################
,'ZnCdS':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'ZnS'
    ,'binary_1_x' :'CdS'
}
#

##########################################################################################
##########################################################################################
#          WURTZITE
##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
#          I N S U L A T O R S    A N D    M E T A L S
##########################################################################################
##########################################################################################

######### air #############################################



##########################################################################################
##########################################################################################
#          T E R N A R Y     A L L O Y S     --    III - V       V A L E N C E
##########################################################################################
##########################################################################################




######### indium gallium nitride (InGaN) ##############################
,'InGaN':{
    'mat_crys_strc'       :'Wurtzite'
    ,'valence'    :'III_V'
    ,'binary_x'   :'InN'
    ,'binary_1_x' :'GaN'

    ,'pyroelectric_const' :-0.037                                # Vurgaftman2 and O. Ambacher (Ambacher has different sign in bowing formula and bowing parameter which is okay.)

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :1.4 }                                 # Vurgaftman2
    }
}
#

######### aluminum gallium nitride (AlGaN) ############################
,'AlGaN':{
    'mat_crys_strc'       :'Wurtzite'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlN'
    ,'binary_1_x' :'GaN'

    ,'pyroelectric_const' :-0.021                                # Vurgaftman2 and O. Ambacher (Ambacher has different sign in bowing formula and bowing parameter which is okay.)

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.7 }                                 # Vurgaftman2
    }
}
#

######### aluminum indium nitride (AlInN) #############################
,'AlInN':{
    'mat_crys_strc'       :'Wurtzite'
    ,'valence'    :'III_V'
    ,'binary_x'   :'AlN'
    ,'binary_1_x' :'InN'

    ,'pyroelectric_const' :-0.070                                # Vurgaftman2 and O. Ambacher (Ambacher has different sign in bowing formula and bowing parameter which is okay.)

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :2.5 }                                 # Vurgaftman2
    }
}
#

######### magnesium zinc oxide (MgZnO) #############################
,'MgZnO':{
    'mat_crys_strc'       :'Wurtzite'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'MgO'
    ,'binary_1_x' :'ZnO'

    ,'conduction_bands':{
        'Gamma':{ 'bandgap' :0.56 }                                 # W.R.L. Lambrecht et al., MRS Internet J. Nitride Semicond. Res. 4S1, G6.8 (1999)
    }

}
#

######### cadmium zinc oxide (CdZnO) #############################
,'CdZnO':{
    'mat_crys_strc'       :'Wurtzite'
    ,'valence'    :'II_VI'
    ,'binary_x'   :'CdO'
    ,'binary_1_x' :'ZnO'
}
#

}

alloyproperty4 ={
##########################################################################################
##########################################################################################
#          Q U A T E R N A R Y     A L L O Y S     --    III - V       V A L E N C E
##########################################################################################
##########################################################################################


######### aluminum gallium indium nitride (AlGaInN) ###################
'AlGaInN_zb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlN_zb'
    ,'binary2'    :'GaN_zb'
    ,'binary3'    :'InN_zb'
    ,'ternary12'  :'AlGaN_zb'
    ,'ternary13'  :'AlInN_zb'
    ,'ternary23'  :'InGaN_zb'
}
#

######### aluminum gallium indium phosphide (AlGaInP) #################
,'AlGaInP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlP'
    ,'binary2'    :'GaP'
    ,'binary3'    :'InP'
    ,'ternary12'  :'AlGaP'
    ,'ternary13'  :'AlInP'
    ,'ternary23'  :'GaInP'
}
#

######### aluminum gallium indium arsenide (AlGaInAs) #################
,'AlGaInAs':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlAs'
    ,'binary2'    :'GaAs'
    ,'binary3'    :'InAs'
    ,'ternary12'  :'AlGaAs'
    ,'ternary13'  :'InAlAs'
    ,'ternary23'  :'InGaAs'
}
#

######### aluminum gallium indium antimonide (AlGaInSb) ###############
,'AlGaInSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlSb'
    ,'binary2'    :'GaSb'
    ,'binary3'    :'InSb'
    ,'ternary12'  :'AlGaSb'
    ,'ternary13'  :'AlInSb'
    ,'ternary23'  :'GaInSb'
}
#

######### aluminum arsenide antimonide phosphide (AlAsSbP) ############
,'AlAsSbP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlAs'
    ,'binary2'    :'AlSb'
    ,'binary3'    :'AlP'
    ,'ternary12'  :'AlAsSb'
    ,'ternary13'  :'AlAsP'
    ,'ternary23'  :'AlPSb'
}
#

######### gallium arsenide antimonide phosphide (GaAsSbP) #############
,'GaAsSbP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'GaAs'
    ,'binary2'    :'GaSb'
    ,'binary3'    :'GaP'
    ,'ternary12'  :'GaAsSb'
    ,'ternary13'  :'GaAsP'
    ,'ternary23'  :'GaPSb'
}
#

######### indium arsenide antimonide phosphide (InAsSbP) ##############
,'InAsSbP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'InAs'
    ,'binary2'    :'InSb'
    ,'binary3'    :'InP'
    ,'ternary12'  :'InAsSb'
    ,'ternary13'  :'InAsP'
    ,'ternary23'  :'InPSb'
}
#

######### gallium indium arsenide phosphide (GaInAsP) #################
,'GaInAsP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'GaAs'
    ,'binary2'    :'InAs'
    ,'binary3'    :'InP'
    ,'binary4'    :'GaP'
    ,'ternary12'  :'InGaAs'
    ,'ternary23'  :'InAsP'
    ,'ternary34'  :'GaInP'
    ,'ternary14'  :'GaAsP'
}
#

######### indium gallium arsenide nitride (InGaAsN) ###################
,'InGaAsN':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'InAs'
    ,'binary2'    :'GaAs'
    ,'binary3'    :'GaN_zb'
    ,'binary4'    :'InN_zb'
    ,'ternary12'  :'InGaAs'
    ,'ternary23'  :'GaAsN'
    ,'ternary34'  :'InGaN_zb'
    ,'ternary14'  :'InAsN'
}
#

######### gallium indium arsenide antimonide (GaInAsSb) ###############
,'GaInAsSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'GaAs'
    ,'binary2'    :'InAs'
    ,'binary3'    :'InSb'
    ,'binary4'    :'GaSb'
    ,'ternary12'  :'InGaAs'
    ,'ternary23'  :'InAsSb'
    ,'ternary34'  :'GaInSb'
    ,'ternary14'  :'GaAsSb'
}
#

######### aluminum gallium arsenide antimonide (AlGaAsSb) #############
,'AlGaAsSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlAs'
    ,'binary2'    :'GaAs'
    ,'binary3'    :'GaSb'
    ,'binary4'    :'AlSb'
    ,'ternary12'  :'AlGaAs'
    ,'ternary23'  :'GaAsSb'
    ,'ternary34'  :'AlGaSb'
    ,'ternary14'  :'AlAsSb'
}
#

######### indium aluminum arsenide antimonide (InAlAsSb) ##############
,'InAlAsSb':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'InAs'
    ,'binary2'    :'AlAs'
    ,'binary3'    :'AlSb'
    ,'binary4'    :'InSb'
    ,'ternary12'  :'InAlAs'
    ,'ternary23'  :'AlAsSb'
    ,'ternary34'  :'AlInSb'
    ,'ternary14'  :'InAsSb'
}
#

######### aluminum gallium arsenide phosphide (AlGaAsP) ##############
,'AlGaAsP':{
    'mat_crys_strc'       :'Zincblende'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlAs'
    ,'binary2'    :'GaAs'
    ,'binary3'    :'GaP'
    ,'binary4'    :'AlP'
    ,'ternary12'  :'AlGaAs'
    ,'ternary23'  :'GaAsP'
    ,'ternary34'  :'AlGaP'
    ,'ternary14'  :'AlAsP'
}
#




##########################################################################################
##########################################################################################
#          Q U A T E R N A R Y     A L L O Y S     --    III - V       V A L E N C E
##########################################################################################
##########################################################################################




######### aluminum gallium indium nitride (AlGaInN) ###################
,'AlGaInN':{
    'mat_crys_strc'       :'Wurtzite'
    ,'valence'    :'III_V'
    ,'binary1'    :'AlN'
    ,'binary2'    :'GaN'
    ,'binary3'    :'InN'
    ,'ternary12'  :'AlGaN'
    ,'ternary13'  :'AlInN'
    ,'ternary23'  :'InGaN'
}
#
}

