


type_of_method = {
                0 : "Hartree-Fock"                                     , 
                1 : "Local and gradient corrected functionals"         , 
                2 : "Hybrid functionals"                               ,
                3 : "Meta-GGA and hybrid meta-GGA's"                   ,
                4 : "Perturbatively corrected double hybrid functional",
                5 : "Second Order Many Body Perturbation Theory"       ,
                6 : "High-level Single Reference Methods"              ,
                7 : "Semiempirical Methods"                            ,
                }


type_of_method_short = {

        0 : "HF", 
        1 : "DFT", 
        2 : "Hybrid",
        3 : "Meta-GGA",
        4 : "D-Hybrid",
        5 : "MP2",
        6 : "HL-SR",
        7 : "SemiEmp",
        }
        

methods = {
        0 : ['HF'], 
        
        1 : [
            "HFS Hartree-Fock-Slater Exchange only functional",
            "LDA Local density approximation",
            "VWN5 Vosko-Wilk-Nusair local density approx",
            "VWN3 Vosko-Wilk-Nusair local density approx.",
            "PWLDA Perdew-Wang parameterization of LDA",
            "BP86 Becke '88 exchange and Perdew '86 correlation",
            "BLYP Becke '88 exchange and Lee-Yang-Parr correlation",
            "OLYP Handy's 'optimal' exchange and Lee-Yang-Parr correlation",
            "GLYP Gill's '96 exchange and Lee-Yang-Parr correlation",
            "XLYP The Xu and Goddard exchange and Lee-Yang-Parr correlation",
            "PW91 Perdew-Wang '91 GGA functional",
            "mPWPW Modified PW exchange and PW correlation",
            "mPWLYP Modified PW exchange and LYP correlation",
            "PBE Perdew-Burke-Erzerhoff GGA functional",
            "RPBE 'Modified' PBE",
            "REVPBE 'Revised' PBE",
            "PWP Perdew-Wang '91 exchange and Perdew 86 correlation"
            ],

        2 : [
            "B1LYP The one-parameter hybrid Becke'88 exchange and LYP correlation (25% HF exchange)",
            "B3LYP The popular B3LYP functional (20% HF exchange)",
            "O3LYP The Handy hybrid functional",
            "X3LYP The Xu and Goddard hybrid functional",
            "B1P The one parameter hybrid version of BP86",
            "B3P The three parameter hybrid version of BP86",
            "B3PW The three parameter hybrid version of PW91",
            "PW1PW One parameter hybrid version of PW91",
            "mPW1PW One parameter hybrid version of mPWPW",
            "mPW1LYP One parameter hybrid version of mPWLYP",
            "PBE0 One parameter hybrid version of PBE",
            ],
        
        3: [
        
            "TPSS The TPSS meta-GGA functional",
            "TPSSh The hybrid version of TPSS",
            "TPSS0 A 25% exchange version of TPSShi",
            ],
        
        4: [
            "B2PLYP The new mixture of MP2 and DFT from Grimme",
            "RI-B2PLYP B2PLYP with RI applied to the MP2 part",
            "B2PLYP-D B2PLYP with Van der Waals correction",
            "RI-B2PLYP RIJONX The same but with RI also applied in the SCF part",
            "mPW2PLYP mPW exchange instead of B88 (also with RI and RIJONX as above for B2PYLP)",
            "mPW is supposed to improve on weak interactions",
            "B2GP-PLYP Gershom Martin's 'general purpose' reparameterization",
            "B2K-PLYP Gershom Martin's 'kinetic' reparameterization",
            "B2T-PLYP Gershom Martin's 'thermochemistry' reparameterization",
            ],            
            
        5: [
            "MP2",
            "RI-MP2",
            "SCS-MP2 Spin-component scaled MP2",
            "RI-SCS-MP2 Spin-component scaled RI-MP2",
            ],
    
            
           

         6: [           
            "CCSD  Coupled cluster singles and doubles",
            "CCSD(T) Same with perturbative triples correction",
            "QCISD Quadratic Configuration interaction",
            "QCISD(T) Same with perturbative triples correction",
            "CPF/1 Coupled pair functional",
            "NCPF/1 A 'new' modified coupled pair functional",
            "CEPA/1 Coupled electron pair approximation",
            "NCEPA/1 The CEPA analogue of NCPF/1",
            "MP3 MP3 energies",
            "SCS-MP3 Grimme's refined version of MP3",
            ],
            
            

         7: [             
            "ZINDO/S",
            "ZINDO/1",
            "ZINDO/2",
            "NDDO/1",
            "NDDO/2",
            "MNDO",
            "AM1",
            "PM3",
            ],            
}

list_scf_convergence = {
        0 : "Default", 
        1 : "TightSCF : tight SCF convergence",
        2 : "LooseSCF : loose SCF convergence",
        3 : "SloppySCF: sloppy SCF convergence",
        4 : "StrongSCF: strong SCF convergence",
        5 : "VeryTightSCF :very tight SCF convergence",
        6 : "ScfConv6 : energy convergence check and ETol=10-6",
        7 : "ScfConv7 : energy convergence check and ETol=10-7",
        8 : "ScfConv8 : energy convergence check and ETol=10-8",
        9 : "ScfConv9 : energy convergence check and ETol=10-9",
        10: "ScfConv10: energy convergence check and ETol=10-10",
}



list_semiemp_excited =  { 
                       0 : "Nothing",
                       1 : "CIS",
                         } 


list_HF_excited = { 
                    0 : "Nothing",
                    1 : "CIS",
                    2 : "CIS(D)",
                    }


list_DFT_excited = { 
                   0 : "Nothing",
                   1 : "TD-DFT",
                    }


basis_set_group_dict = { 
                          0 : "Pople Style basis sets", 
                          1 : "Pople with one diffuse function on non-hydrogen atoms", 
                          2 : "Pople with one diffuse function on all atoms", 
                          3 : "Dunning basis sets", 
                          4 : "Ahlrichs basis sets", 
                          5 : "Def2 Ahlrichs basis sets", 
                          6 : "Jensen Basis Sets", 
                          7 : "Atomic Natural Orbital Basis Sets", 
                          8 : "Miscellenous and Specialized Basis Sets", 
                        }



list_Pople_Basis = [
    "3-21G Pople 3-21G",
    "3-21GSP Buenker 3-21GSP",
    "4-22GSP Buenker 4-22GSP",
    "6-31G Pople 6-31G and its modifications",
    "6-311G Pople 6-311G and its modifications",
    
    "3-21G*  3-21G plus one polarisation function all non-hydrogens atoms",
    "3-21GSP* 3-21GSP plus one polarisation function all non-hydrogens atoms",
    "4-22GSP* 4-22GSP plus one polarisation function all non-hydrogens atoms",
    "6-31G* 6-31G plus one polarisation function all non-hydrogens atoms",
    "6-311G* 6-311G plus one polarisation function all non-hydrogens atoms",
    
    "3-21G**  3-21G plus one polarisation function all atoms",
    "3-21GSP** 3-21GSP plus one polarisation function all atoms",
    "4-22GSP** 4-22GSP plus one polarisation function all atoms",
    "6-31G** 6-31G plus one polarisation function all atoms",
    "6-311G** 6-311G plus one polarisation function all atoms",
    
    "3-21G(2d)  3-21G plus two polarisation functions all non-hydrogens atoms",
    "3-21GSP(2d) 3-21GSP plus two polarisation functions all non-hydrogens atoms",
    "4-22GSP(2d) 4-22GSP plus two polarisation functions all non-hydrogens atoms",
    "6-31G(2d) 6-31G plus two polarisation functions all non-hydrogens atoms",
    "6-311G(2d) 6-311G plus two polarisation functions all non-hydrogens atoms",
    
    "3-21G(2d,2p)  3-21G plus two polarisation functions all atoms",
    "3-21GSP(2d,2p) 3-21GSP plus two polarisation functions all atoms",
    "4-22GSP(2d,2p) 4-22GSP plus two polarisation functions all atoms",
    "6-31G(2d,2p) 6-31G plus two polarisation functions all atoms",
    "6-311G(2d,2p) 6-311G plus two polarisation functions all atoms",
    
    "3-21G(2df)  3-21G plus three polarisation functions all non-hydrogens atoms",
    "3-21GSP(2df) 3-21GSP plus three polarisation functions all non-hydrogens atoms",
    "4-22GSP(2df) 4-22GSP plus three polarisation functions all non-hydrogens atoms",
    "6-31G(2df) 6-31G plus three polarisation functions all non-hydrogens atoms",
    "6-311G(2df) 6-311G plus three polarisation functions all non-hydrogens atoms",
    
    "3-21G(2df,2pd)  3-21G plus three polarisation functions all atoms",
    "3-21GSP(2df,2pd) 3-21GSP plus three polarisation functions all atoms",
    "4-22GSP(2df,2pd) 4-22GSP plus three polarisation functions all atoms",
    "6-31G(2df,2pd) 6-31G plus three polarisation functions all atoms",
    "6-311G(2df,2pd) 6-311G plus three polarisation functions all atoms",
    
    "3-21G(3df)  3-21G plus four polarisation functions all non-hydrogens atoms",
    "3-21GSP(3df) 3-21GSP plus four polarisation functions all non-hydrogens atoms",
    "4-22GSP(3df) 4-22GSP plus four polarisation functions all non-hydrogens atoms",
    "6-31G(3df) 6-31G plus four polarisation functions all non-hydrogens atoms",
    "6-311G(3df) 6-311G plus four polarisation functions all non-hydrogens atoms",
    
    "3-21G(3df,3pd)  3-21G plus four polarisation functions all atoms",
    "3-21GSP(3df,3pd) 3-21GSP plus four polarisation functions all atoms",
    "4-22GSP(3df,3pd) 4-22GSP plus four polarisation functions all atoms",
    "6-31G(3df,3pd) 6-31G plus four polarisation functions all atoms",
    "6-311G(3df,3pd) 6-311G plus four polarisation functions all atoms",
]

list_Pople_Diffuse_Non_Hydrogen_Basis = [
    "3-21+G  3-21G plus diffuse functions on all non-hydrogens atoms",
    "3-21+GSP 3-21GSP plus diffuse functions on all non-hydrogens atoms",
    "4-22+GSP 4-22GSP plus diffuse functions on all non-hydrogens atoms",
    "6-31+G 6-31G plus diffuse functions on all non-hydrogens atoms",
    "6-311+G 6-311G plus diffuse functions on all non-hydrogens atoms",

    "3-21+G*  3-21G + diff. non-hydrogens + 1 pol. non-hydrogens",
    "3-21+GSP* 3-21GSP + diff. non-hydrogens + 1 pol. non-hydrogens",
    "4-22+GSP* 4-22GSP + diff. non-hydrogens + 1 pol. non-hydrogens",
    "6-31+G* 6-31G + diff. non-hydrogens + 1 pol. non-hydrogens",
    "6-311+G* 6-311G + diff. non-hydrogens + 1 pol. non-hydrogens",

    "3-21+G**  3-21G + diff. non-hydrogens + 1 pol. ",
    "3-21+GSP** 3-21GSP + diff. non-hydrogens + 1 pol. ",
    "4-22+GSP** 4-22GSP + diff. non-hydrogens + 1 pol. ",
    "6-31+G** 6-31G + diff. non-hydrogens + 1 pol. ",
    "6-311+G** 6-311G + diff. non-hydrogens + 1 pol. ",

    "3-21+G(2d)  3-21G + diff. non-hydrogens + 2 pol. non-hydrogens",
    "3-21+GSP(2d) 3-21GSP + diff. non-hydrogens + 2 pol. non-hydrogens",
    "4-22+GSP(2d) 4-22GSP + diff. non-hydrogens + 2 pol. non-hydrogens",
    "6-31+G(2d) 6-31G + diff. non-hydrogens + 2 pol. non-hydrogens",
    "6-311+G(2d) 6-311G + diff. non-hydrogens + 2 pol. non-hydrogens",

    "3-21+G(2d,2p)  3-21G + diff. non-hydrogens + 2 pol.",
    "3-21+GSP(2d,2p) 3-21GSP + diff. non-hydrogens + 2 pol.",
    "4-22+GSP(2d,2p) 4-22GSP + diff. non-hydrogens + 2 pol.",
    "6-31+G(2d,2p) 6-31G + diff. non-hydrogens + 2 pol.",
    "6-311+G(2d,2p) 6-311G + diff. non-hydrogens + 2 pol.",

    "3-21+G(2df)  3-21G + diff. non-hydrogens + 3 pol. non-hydrogens",
    "3-21+GSP(2df) 3-21GSP + diff. non-hydrogens + 3 pol. non-hydrogens",
    "4-22+GSP(2df) 4-22GSP + diff. non-hydrogens + 3 pol. non-hydrogens",
    "6-31+G(2df) 6-31G + diff. non-hydrogens + 3 pol. non-hydrogens",
    "6-311+G(2df) 6-311G + diff. non-hydrogens + 3 pol. non-hydrogens",

    "3-21+G(2df,2pd)  3-21G + diff. non-hydrogens + 3 pol.",
    "3-21+GSP(2df,2pd) 3-21GSP + diff. non-hydrogens + 3 pol.",
    "4-22+GSP(2df,2pd) 4-22GSP + diff. non-hydrogens + 3 pol.",
    "6-31+G(2df,2pd) 6-31G + diff. non-hydrogens + 3 pol.",
    "6-311+G(2df,2pd) 6-311G + diff. non-hydrogens + 3 pol.",

    "3-21+G(3df)  3-21G + diff. non-hydrogens + 4 pol. non-hydrogens",
    "3-21+GSP(3df) 3-21GSP + diff. non-hydrogens + 4 pol. non-hydrogens",
    "4-22+GSP(3df) 4-22GSP + diff. non-hydrogens + 4 pol. non-hydrogens",
    "6-31+G(3df) 6-31G + diff. non-hydrogens + 4 pol. non-hydrogens",
    "6-311+G(3df) 6-311G + diff. non-hydrogens + 4 pol. non-hydrogens",

    "3-21+G(3df,3pd)  3-21G + diff. non-hydrogens + 4 pol.",
    "3-21+GSP(3df,3pd) 3-21GSP + diff. non-hydrogens + 4 pol.",
    "4-22+GSP(3df,3pd) 4-22GSP + diff. non-hydrogens + 4 pol.",
    "6-31+G(3df,3pd) 6-31G + diff. non-hydrogens + 4 pol.",
    "6-311+G(3df,3pd) 6-311G + diff. non-hydrogens + 4 pol.",
]

list_Pople_Diffuse_All_Atoms_Basis = [
    "3-21++G  3-21G plus diffuse functions on all atoms",
    "3-21++GSP 3-21GSP plus diffuse functions on all atoms",
    "4-22++GSP 4-22GSP plus diffuse functions on all atoms",
    "6-31++G 6-31G plus diffuse functions on all atoms",
    "6-311++G 6-311G plus diffuse functions on all atoms",
    
    "3-21++G*  3-21G + diff. + 1 pol. non-hydrogens",
    "3-21++GSP* 3-21GSP + diff. + 1 pol. non-hydrogens",
    "4-22++GSP* 4-22GSP + diff. + 1 pol. non-hydrogens",
    "6-31++G* 6-31G + diff. + 1 pol. non-hydrogens",
    "6-311++G* 6-311G + diff. + 1 pol. non-hydrogens",
    
    "3-21++G**  3-21G + diff. + 1 pol. ",
    "3-21++GSP** 3-21GSP + diff. + 1 pol. ",
    "4-22++GSP** 4-22GSP + diff. + 1 pol. ",
    "6-31++G** 6-31G + diff. + 1 pol. ",
    "6-311++G** 6-311G + diff. + 1 pol. ",
    
    "3-21++G(2d)  3-21G + diff. + 2 pol. non-hydrogens",
    "3-21++GSP(2d) 3-21GSP + diff. + 2 pol. non-hydrogens",
    "4-22++GSP(2d) 4-22GSP + diff. + 2 pol. non-hydrogens",
    "6-31++G(2d) 6-31G + diff. + 2 pol. non-hydrogens",
    "6-311++G(2d) 6-311G + diff. + 2 pol. non-hydrogens",
    
    "3-21++G(2d,2p)  3-21G + diff. + 2 pol.",
    "3-21++GSP(2d,2p) 3-21GSP + diff. + 2 pol.",
    "4-22++GSP(2d,2p) 4-22GSP + diff. + 2 pol.",
    "6-31++G(2d,2p) 6-31G + diff. + 2 pol.",
    "6-311++G(2d,2p) 6-311G + diff. + 2 pol.",
    
    "3-21++G(2df)  3-21G + diff. + 3 pol. non-hydrogens",
    "3-21++GSP(2df) 3-21GSP + diff. + 3 pol. non-hydrogens",
    "4-22++GSP(2df) 4-22GSP + diff. + 3 pol. non-hydrogens",
    "6-31++G(2df) 6-31G + diff. + 3 pol. non-hydrogens",
    "6-311++G(2df) 6-311G + diff. + 3 pol. non-hydrogens",
    
    "3-21++G(2df,2pd)  3-21G + diff. + 3 pol.",
    "3-21++GSP(2df,2pd) 3-21GSP + diff. + 3 pol.",
    "4-22++GSP(2df,2pd) 4-22GSP + diff. + 3 pol.",
    "6-31++G(2df,2pd) 6-31G + diff. + 3 pol.",
    "6-311++G(2df,2pd) 6-311G + diff. + 3 pol.",
    
    "3-21++G(3df)  3-21G + diff. + 4 pol. non-hydrogens",
    "3-21++GSP(3df) 3-21GSP + diff. + 4 pol. non-hydrogens",
    "4-22++GSP(3df) 4-22GSP + diff. + 4 pol. non-hydrogens",
    "6-31++G(3df) 6-31G + diff. + 4 pol. non-hydrogens",
    "6-311++G(3df) 6-311G + diff. + 4 pol. non-hydrogens",
    
    "3-21++G(3df,3pd)  3-21G + diff. + 4 pol.",
    "3-21++GSP(3df,3pd) 3-21GSP + diff. + 4 pol.",
    "4-22++GSP(3df,3pd) 4-22GSP + diff. + 4 pol.",
    "6-31++G(3df,3pd) 6-31G + diff. + 4 pol.",
    "6-311++G(3df,3pd) 6-311G + diff. + 4 pol.",
]

list_Dunning_Basis  = [
    "cc-pVDZ Dunning correlation concisistent polarized double zeta",
    "cc-(p)VDZ Same but no polarization on hydrogens",
    "Aug-cc-pVDZ Same but including diffuse functions",
    "cc-pVTZ Dunning correlation concisistent polarized triple zeta",
    "cc-(p)VTZ Same but no polarization on hydrogen",
    "Aug-cc-pVTZ Same but including diffuse functions(g-functions deleted!)",
    "cc-pVQZ Dunning correlation concisistent polarized quadruple zeta",
    "Aug-cc-pVQZ with diffuse functions",
    "cc-pV5Z Dunning correlation concisistent polarized quintuple zeta",
    "Aug-cc-pV5Z with diffuse functions",
    "cc-pV6Z Dunning correlation concisistent polarized sextuple zeta",
    "Aug-cc-pV6Z ... with diffuse functions",
    "cc-pCVDZ Core-polarized double-zeta correlation consistent basis set",
    "cc-pCVTZ Same for triple zeta",
    "cc-pCVQZ Same for quadruple zeta",
    "cc-pCV5Z Same for quintuple zeta",
    "cc-pV6Z Same for sextuple zeta",
    "Aug-pCVDZ Same double zeta with diffuse functions augmented",
    "Aug-pCVTZ Same for triple zeta",
    "Aug-pCVQZ Same for quadruple zeta",
    "Aug-pCV5Z Same for quintuple zeta",
    "Aug-cc-pV6Z Same for sextuple zeta",
    "DUNNING-DZP Dunning's original double zeta basis set",
    ]


list_Ahlrichs_Basis = [
    "SV Ahlrichs split valence basis set",
    "VDZ Ahlrichs split valence basis set",
    "VTZ Ahlrichs Valence triple zeta basis set",
    "TZV Ahlrichs triple-zeta valence basis set. NOT identical to VTZ",
    "QZVP Ahlrichs quadruple-zeta basis set. P is already polarized",
    "DZ Ahlrichs double zeta basis set",
    "QZVPP(-g,-f) QZVPP with highest polarization functions deleted",

    "SV(P) SV + One polar set on all non-hydrogens atoms",
    "VDZ(P) VDZ + One polar set on all non-hydrogens atoms",
    "VTZ(P) VTZ  + One polar set on all non-hydrogens atoms",
    "TZV(P) TZV  + One polar set on all non-hydrogens atoms",
    "DZ(P)  DZ   + One polar set on all non-hydrogens atoms",

    "SVP SV + One polar set on all atoms",
    "VDZP VDZ + One polar set on all atoms",
    "VTZP VTZ  + One polar set on all atoms",
    "TZVP TZV  + One polar set on all atoms",
    "DZP  DZ   + One polar set on all atoms",

    "SV(2D) SV + Two polar set on all non-hydrogens atoms",
    "VDZ(2D) VDZ + Two polar set on all non-hydrogens atoms",
    "VTZ(2D) VTZ  + Two polar set on all non-hydrogens atoms",
    "TZV(2D) TZV  + Two polar set on all non-hydrogens atoms",
    "DZ(2D)  DZ   + Two polar set on all non-hydrogens atoms",

    "SV(2D,2P) SV + Two polar set on all atoms",
    "VDZ(2D,2P) VDZ + Two polar set on all atoms",
    "VTZ(2D,2P) VTZ  + Two polar set on all atoms",
    "TZV(2D,2P) TZV  + Two polar set on all atoms",
    "DZ(2D,2P)  DZ   + Two polar set on all atoms",

    "SV(2df) SV + Three polar set on all non-hydrogens atoms",
    "VDZ(2df) VDZ + Three polar set on all non-hydrogens atoms",
    "VTZ(2df) VTZ  + Three polar set on all non-hydrogens atoms",
    "TZV(2df) TZV  + Three polar set on all non-hydrogens atoms",
    "DZ(2df)  DZ   + Three polar set on all non-hydrogens atoms",

    "SV(2df,2pd) SV + Three polar set on all atoms",
    "VDZ(2df,2pd) VDZ + Three polar set on all atoms",
    "VTZ(2df,2pd) VTZ  + Three polar set on all atoms",
    "TZV(2df,2pd) TZV  + Three polar set on all atoms",
    "DZ(2df,2pd)  DZ   + Three polar set on all atoms",

    "SV(PP) SV + Three polar set on all non-hydrogens atoms",
    "VDZ(PP) VDZ + Three polar set on all non-hydrogens atoms",
    "VTZ(PP) VTZ  + Three polar set on all non-hydrogens atoms",
    "TZV(PP) TZV  + Three polar set on all non-hydrogens atoms",
    "DZ(PP)  DZ   + Three polar set on all non-hydrogens atoms",

    "SVPP SV + Three polar set on all atoms",
    "VDZPP VDZ + Three polar set on all atoms",
    "VTZPP VTZ  + Three polar set on all atoms",
    "TZVPP TZV  + Three polar set on all atoms",
    "DZPP  DZ   + Three polar set on all atoms",

    "SV(P)+ SV plus Pople diff. func. + 1 polar on non-hydrogens",
    "VDZ(P)+ VDZ plus Pople diff. func. + 1 polar on non-hydrogens",
    "VTZ(P)+ VTZ plus Pople diff. func. + 1 polar on non-hydrogens",
    "TZV(P)+ TZV  plus Pople diff. func. + 1 polar on non-hydrogens",

    "SVP++ SV plus Pople diff. func. + One polar set on all atoms",
    "TZVP++ TZV  plus Pople diff. func. + One polar set on all atoms",
    "TZV(2D) TZV  lus Pople diff. func. + One polar set on all atoms",

    "aug-SV(P) SV plus Dunning diff. One polar set on all non-hydrogens atoms",
    "aug-VDZ(P) VDZ plus Dunning diff. One polar set on all non-hydrogens atoms",
    "aug-VTZ(P) VTZ  plus Dunning diff. One polar set on all non-hydrogens atoms",
    "aug-TZV(P) TZV  plus Dunning diff. One polar set on all non-hydrogens atoms",
    "aug-DZ(P)  DZ   plus Dunning diff. One polar set on all non-hydrogens atoms",

    "aug-SVP SV plus Dunning diff. One polar set on all atoms",
    "aug-VDZP VDZ plus Dunning diff. One polar set on all atoms",
    "aug-VTZP VTZ  plus Dunning diff. One polar set on all atoms",
    "aug-TZVP TZV  plus Dunning diff. One polar set on all atoms",
    "aug-DZP  DZ   plus Dunning diff. One polar set on all atoms",

    "aug-SV(2D) SV plus Dunning diff. Two polar set on all non-hydrogens atoms",
    "aug-VDZ(2D) VDZ plus Dunning diff. Two polar set on all non-hydrogens atoms",
    "aug-VTZ(2D) VTZ  plus Dunning diff. Two polar set on all non-hydrogens atoms",
    "aug-TZV(2D) TZV  plus Dunning diff. Two polar set on all non-hydrogens atoms",
    "aug-DZ(2D)  DZ   plus Dunning diff. Two polar set on all non-hydrogens atoms",

    "aug-SV(2D,2P) SV plus Dunning diff. Two polar set on all atoms",
    "aug-VDZ(2D,2P) VDZ plus Dunning diff. Two polar set on all atoms",
    "aug-VTZ(2D,2P) VTZ  plus Dunning diff. Two polar set on all atoms",
    "aug-TZV(2D,2P) TZV  plus Dunning diff. Two polar set on all atoms",
    "aug-DZ(2D,2P)  DZ   plus Dunning diff. Two polar set on all atoms",

    "aug-SV(2df) SV plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-VDZ(2df) VDZ plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-VTZ(2df) VTZ  plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-TZV(2df) TZV  plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-DZ(2df)  DZ   plus Dunning diff. Three polar set on all non-hydrogens atoms",

    "aug-SV(2df,2pd) SV plus Dunning diff. Three polar set on all atoms",
    "aug-VDZ(2df,2pd) VDZ plus Dunning diff. Three polar set on all atoms",
    "aug-VTZ(2df,2pd) VTZ  plus Dunning diff. Three polar set on all atoms",
    "aug-TZV(2df,2pd) TZV  plus Dunning diff. Three polar set on all atoms",
    "aug-DZ(2df,2pd)  DZ   plus Dunning diff. Three polar set on all atoms",

    "aug-SV(PP) SV plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-VDZ(PP) VDZ plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-VTZ(PP) VTZ  plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-TZV(PP) TZV  plus Dunning diff. Three polar set on all non-hydrogens atoms",
    "aug-DZ(PP)  DZ   plus Dunning diff. Three polar set on all non-hydrogens atoms",

    "aug-SVPP SV plus Dunning diff. Three polar set on all atoms",
    "aug-VDZPP VDZ plus Dunning diff. Three polar set on all atoms",
    "aug-VTZPP VTZ  plus Dunning diff. Three polar set on all atoms",
    "aug-TZVPP TZV  plus Dunning diff. Three polar set on all atoms",
    "aug-DZPP  DZ   plus Dunning diff. Three polar set on all atoms",
    ]


list_Def2Ahlrichs_Basis =[
    "Def2-SV(P) SV basis set with 'new' polarization functions",
    "Def2-SVP",
    "Def2-TZVP TZVP basis set with 'new' polarization functions",
    "Def2-TZVP(-f) Delete the f-polarization functions from def2-TZVP",
    "Def2-TZVP(-df) delete the double d-function and replace it by the older single d-function.",
    "Def2-TZVPP TZVPP basis set with 'new' polarization functions",
    "Def2-aug-TZVPP Same but with diffuse functions from aug-cc-pVTZ",
    "Def2-QZVPP Very accurate quadruple-zeta basis.",
    "Def2-QZVPP(-g,-f) highest angular momentum polarization functions deleted",
]

list_Jensen_Basis =  [
	"PC-1 Polarization consistent basis sets (H-Ar) optimized for DFT",
	"PC-2 double zeta polarization consistent basis sets (H-Ar) optimized for DFT",
	"PC-3 triple zeta polarization consistent basis sets (H-Ar) optimized for DFT",
	"PC-4 quadruple zeta polarization consistent basis sets (H-Ar) optimized for DFT",
	"Aug-PC-1 PC-1 with augmentations by diffuse functions",
	"Aug-PC-2 PC-2 with augmentations by diffuse functions",
	"Aug-PC-3 PC-3 with augmentations by diffuse functions",
	"Aug-PC-4 PC-4 with augmentations by diffuse functions",
]

list_ANO_Basis = [ 
    "ano-pVDZ better than the cc-pVDZ (but much larger number of primitives of course)",
    "ano-pVTZ",
    "ano-pVQZ",
    "ano-pV5Z",
    "saug-ano-pVDZ ano-pVDZ augmentation with a single set of s,p functions.",
    "saug-ano-pVTZ ano-pVTZ augmentation with a single set of s,p functions.",
    "saug-ano-pVQZ ano-pVQZ augmentation with a single set of s,p functions. ",
    "saug-ano-pV5Z ano-pV5Z augmentation with a single set of s,p functions.",
    "aug-ano-pVDZ ano-pVDZ full augmentation with spd",
    "aug-ano-pVTZ ano-pVTZ full augmentation with spdf",
    "aug-ano-pVQZ ano-pVQZ full augmentation with spdfg",
    "BNANO-DZP (Bonn-ANO-DZP), small DZP type ANO basis set from the Bonn group",
    "BNANO-TZ2P (Bonn-ANO-TZ2P), slightly larger triple-zeta ANO with two pol. Bonn group",
    "BNANO-TZ3P Same but with a contracted set of f-polarization functions on the heavy atoms",
    "NASA-AMES-ANO The original NASA/AMES ANO basis set (quadruple-zeta type)",
    "BAUSCHLICHER ANO First row transition metal ANO sets",
    "ROOS-ANO-DZP A fairly large DZP basis set from Roos, same size as aug-ano-pVDZ",
    "ROOS-ANO-TZP A fairly large TZP basis from Roos, same size as aug-ano-pVTZ",
]

list_Miscellenous_Basis = [
    "DGAUSS DGauss polarized valence double zeta basis set",
    "DZVP-DFT DGauss polarized valence double zeta basis set",
    "SADLEJ-PVTZ Sadlej's polarized triple zeta basis for poarlizability and related calculations",
    "EPR-II Barone's Basis set for EPR calculations (double zeta)",
    "EPR-III Barone's Basis set for EPR calculations (triple-zeta)",
    "IGLO-II Kutzelniggs basis set for NMR and EPR calculations",
    "IGLO-III Kutzelniggs basis set for NMR and EPR calculations (accurate)",
    "Partridge-1 Accurate uncontracted basis set",
    "Partridge-2 Accurate uncontracted basis set",
    "Partridge-3 Accurate uncontracted basis set",
    "Wachters Good first row transition metal basis set",
]

listAuxBasisView = [

    "AutoAux Automatic construction of a general purpose fitting basis",
    "DEMON/J The DeMon/J Coulomb fitting basis",
    "DGAUSS/J The DGauss A1 Coulomb fitting basis",
    "SV/J (=VDZ/J) Ahlrichs Coulomb fitting basis for the SVP basis",
    "TZV/J (=VTZ/J) Ahlrichs Coulomb fitting basis for the TZV or TZVP basis",
    "QZVPP/J Ahlrichs Coulomb fitting for the QZVPP basis",
    "Def2-SVP/J Ahlrichs Coulomb fitting for def-SVP",
    "Def2-TZVPP/J Ahlrichs Coulomb fitting for def2-TZVPP/J",
    "Def2-QZVPP/J Ahlrichs Coulomb fitting for def2-QZVPP/J",
    "SV/J(-f) Same as SV/J but with the highest angular momentum aux-function deleted",
    "TZV/J(-f) Same as TZV/J but with the highest angular momentum aux-function deleted",
    "SV/C (=VDZ/C) The Ahlrichs correlation fitting basis for MP2-RI with SVP",
    "TZV/C (=VTZ/C) The Ahlrichs correlation fitting basis for MP2-RI with TZVP",
    "TZVPP/C (=VTZPP/C) The Ahlrichs correlation fitting basis for MP2-RI with extended triple-z bases",
    "QZVP/C Correlation fitting for the QZVP basis",
    "QZVPP/C Correlation fitting for the QZVPP basis",
    "Def2-SVP/C Correlation fitting for the def2-SVP basis",
    "Def2-TZVP/C Correlation fitting for the def2-TZVP basis",
    "Def2-TZVPP/C Correlation fitting for the def2-TZVPP basis",
    "Def2-QZVPP/C Correlation fitting for the def2-QZVPP basis",
    "cc-pVDZ/C Aux-basis for the cc-pVDZ orbital basis",
    "cc-pVTZ/C Aux-basis for the cc-pVTZ orbital basis",
    "cc-pVQZ/C Aux-basis for the cc-pVQZ orbital basis",
    "cc-pV5Z/C Aux-basis for the cc-pV5Z orbital basis",
    "cc-pV6Z/C Aux-basis for the cc-pV6Z orbital basis,",
    "Aug-cc-pVDZ/C Aux-basis for the aug-cc-pVDZ orbital basis",
    "Aug-cc-pVTZ/C Aux-basis for the aug-cc-pVTZ orbital basis",
    "Aug-SV/C Aux basis for SVP and related bases but with diffuse functions",
    "Aug-TZV/C Aux basis for TZVP and related bases but with diffuse functions",
    "Aug-TZVPP/C Aux basis for TZVPP and related bases but with diffuse functions",
    "SVP/JK Coulomb+Exchange fitting for SVP",
    "TZVPP/JK Coulomb+Exchange fitting for TZVPP",
    "QZVPP/JK Coulomb+Exchange fitting for QZVPP",
    "Def2-SVP/JK Coulomb+Exchange fitting for def2-SVP",
    "Def2-TZVPP/JK Coulomb+Exchange fitting for def2-TZVPP",
    "Def2-QZVPP/JK Coulomb+Exchange fitting for def2-QZVPP",
    "cc-pVDZ/JK Coulomb+Exchange fitting for cc-pVDZ",
    "cc-pVTZ/JK Coulomb+Exchange fitting for cc-pVTZ",
    "cc-pVQZ/JK Coulomb+Exchange fitting for cc-pVQZ",
    "cc-pV5Z/JK Coulomb+Exchange fitting for cc-pV5Z",
    "cc-pV6Z/JK Coulomb+Exchange fitting for cc-pV6Z",
]


type_basis_dict = {
                    0 : list_Pople_Basis,
                    1 : list_Pople_Diffuse_Non_Hydrogen_Basis,
                    2 : list_Pople_Diffuse_All_Atoms_Basis,
                    3 : list_Dunning_Basis,
                    4 : list_Ahlrichs_Basis,
                    5 : list_Def2Ahlrichs_Basis,
                    6 : list_Jensen_Basis,
                    7 : list_ANO_Basis,
                    8 : list_Miscellenous_Basis,
                    }











