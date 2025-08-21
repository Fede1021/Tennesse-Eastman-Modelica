model MonliticComplete

//Constant Global Parameters

constant Real RkJ(unit="J/(mol.K)") = 8.31451;
constant Real Rkcal( unit="cal/(mol.K)") = 1.987;
parameter Real Cp_vap[nComp](each unit="kJ/(kg*degC)") = {14.6, 2.04, 1.05, 1.85, 1.87, 2.02, 0.712, 0.628} "Vapor heat capacity";
parameter Real Cp_liq[nComp](each unit="kJ/(kg.degC)") = {0, 0, 0, 7.66, 4.17, 4.45, 2.55, 2.45} "Liquid heat capacity";
constant Integer nComp = 8 "Number of components";
parameter Real M[nComp](each unit="g/mol") = {2, 25.4, 28.0, 32.0, 46.0, 48.0, 62.0, 76.0} "Molecular weight";
constant Integer nChemReact = 3 "Number of chemical reactions";
parameter Real v[nComp,nChemReact] = [-1, -1, -1/3; 0,  0,  0; -1, -1,  0; -1,  0,  -1; 0, -1, -1/3;
                                      0,  0,  1; 1,  0,  0; 0,  1,  0] "Stoichiometric coefficients";
parameter Real Tref( unit="K") = 273.15;
parameter Real alpha[nChemReact](each unit="1") = {1.0399157,1.010373129,1} "Adjustable parameter of reaction rate";
parameter Real cp_cw( unit="kJ/(kg.degC)") = 4.18;
parameter Real rho[nComp](each unit="kg/m3") = {0,0,0,299,365,328,612,617} "Density of liquids";
parameter Real H_vap[nComp](each unit="kJ/kg") = {0,0,0,202,372,372,523,486}"Heat of vaporization";
parameter Real k_comp= 0.7166374645; //Compressor coefficient
parameter Real A[nComp] = {0,0,0,20.81, 21.24, 21.24, 21.32, 22.10};
parameter Real B[nComp] = {0,0,0,-1444, -2114, -2144, -2748, -3318};
parameter Real C[nComp] = {0,0,0,259, 265.5, 265.5, 232.9, 249.6};



//MixingUnit//

parameter Real Vmix(unit="m3") = 141.53 "Mixing zone volume";
Real N_mix[nComp](each unit="kmol",
              start={48.83796012,13.49310238,40.03019459,10.44186782,
              28.48830431,2.514120166,5.403006586,2.517154711},
              each fixed=true)  "Total molar holdup";
Real pMix(unit="MPa");
Real Tmix(unit="K", start=359.25, fixed=true);
Real concentrationMix[nComp];

//Valve 1//

parameter Real Cv1(unit = "kmol/s/MPa**(1/2)") = 0.83337;

//Reactor//

Real P_rea_tot(unit="MPa");
Real N_reactor[nComp](each unit="kmol",
              start={5.2046,2.228951,4.65083,0.11534,7.44601,1.14649,56.0529,59.82564},
              each fixed=true)  "Total molar holdup";
Real Rate[nChemReact](each unit="kmol/s") "Reaction Rate";
Real T_rea( unit="K", start=393.55, fixed=true) "Reactor temperature";
Real Q_rea "Cooling heat flow";
Real delt_Hr[nChemReact];
parameter Real Vr( unit="m3") = 36.8 "Reactor volume";
Real Vliq_rea( unit="m3") "Liquid volume";
Real Vvap_rea( unit="m3") "Vapor volume";
Real p_rea[nComp](each unit="MPa") "Reactor partial pressure";
parameter Real m_CWS_rea( unit="mol/s") = 93.37/3600*1000/18;
parameter Real T_CWSs_rea_in( unit="K") = 308;
parameter Real T_CWSs_rea_out( unit="K") = 367.599;
parameter Real gamma_rea[nComp] = {0,0,0,0.996011,1,1.078,0.999,0.999} "Activity coefficient reactor";
Real x_liq_rea[nComp](each unit="1") "Molar fraction of liquid in the reactor";
Real x_vap_rea[nComp](each unit="1") "Molar fraction of liquid in the reactor";
Real pSat_rea[nComp](each unit="MPa") "Partial pressure at reactor temperature";
parameter Real rho_liq_reactor( unit = "kmol/m3")= 9.337145754 "Reactor Liquid density";


//Valve 2//

parameter Real Cv2(unit = "kmol/s/MPa**(1/2)") = 1.5344;

//Separator//

Real N_separator[nComp](each unit="kmol",
              start={27.497976,12.09634,24.57206,0.0836049,
              5.864057,0.9011249,24.137664,19.752554})"Total molar holdup";
Real T_sep( unit="K", start=353.25, fixed=true) "Separator temperature";
parameter Real Vs(unit="m3")= 99.1 "Separator volume";
Real Q_sep "Cooling heat flow";
Real HoVsep "Vaporization heat";
parameter Real m_CWS_sep( unit="mol/s") = 49.37/3600*1000/18;
parameter Real T_CWSs_sep_in( unit="K") = 313;
parameter Real T_CWSs_sep_out( unit="K") = 350.447;
Real Vliq_sep( unit="m3") "Liquid volume";
Real Vvap_sep( unit="m3") "Vapor volume";
Real P_sep_tot(unit="MPa");
parameter Real rho_liq_separator( unit = "kmol/m3")= 10.29397546 "Reactor Liquid density";
parameter Real gamma_sep[nComp] = {0,0,0,1.001383,1.001383,1.091383,1.001383,0.992188}"Activity coefficient separator";
Real x_sep[nComp](each unit="1") "Molar fraction of liquid in the separator";
Real pSat_sep[nComp](each unit="MPa") "Partial pressure at separator temperature";
Real p_sep[nComp](each unit="MPa") "Separator partial pressure";


//Stripper

Real NG(unit="kmol", start = 20.38135);
Real NH( unit = "kmol", start = 17.313);
Real phi[nComp]"Split factor stripper";
Real T_str( unit="K", start=338.85, fixed=true) "Stripper temperature";
Real Q_str "Cooling heat flow";
Real HoVstr "Vaporization heat";
parameter Real m_CWS_str(unit="mol/s") = 230.31/3600;
parameter Real rho_liq_stripper( unit = "kmol/m3")= 8.6496 "Stripper Liquid density";
Real Vliq_str( unit="m3") "Liquid volume";
parameter Real PC_1[nComp] = {0,0,0,3.51e-7, 3.69e-7, 3.69e-7, 1.84e-9, 8.32e-8};
parameter Real PC_2[nComp] = {0,0,0,-0.000111, -0.0001, -0.0001, 1.37e-5, -7.64e-6};
parameter Real PC_3[nComp] = {0,0,0,0.011351, 0.010197, 0.010049, 0.000217, 0.000976};
parameter Real PC_4[nComp] = {0,0,0,0.548012, 0.620794, 0.628854, 0.001393, -0.01568};


//Fluxes Parameters
Real F1 (unit="kmol/h");
Real F2 (unit="kmol/h");
Real F3 (unit="kmol/h");
Real F4 (unit="kmol/h");
Real F5 (unit="kmol/h");
Real F6 (unit="kmol/h");
Real F7 (unit="kmol/h");
Real F8 (unit="kmol/h");
Real F9 (unit="kmol/h");
Real F10 (unit="kmol/h");
Real F11 (unit="kmol/h");        

Real y_F1[nComp];
Real y_F2[nComp];
Real y_F3[nComp];
Real y_F4[nComp];
Real y_F5[nComp];
Real y_F6[nComp];
Real y_F7[nComp];
Real y_F8[nComp];
Real y_F9[nComp];
Real x_F10[nComp];
Real x_F11[nComp];

Real T_F1 (unit = "K");
Real T_F2 (unit = "K");
Real T_F3 (unit = "K");
Real T_F4 (unit = "K");
Real T_F5 (unit = "K");
Real T_F6 (unit = "K");
Real T_F7 (unit = "K");
Real T_F8 (unit = "K");
Real T_F9 (unit = "K");
Real T_F10 (unit = "K");
Real T_F11 (unit = "K");

Real p_F1 (unit = "MPa");
Real p_F2 (unit = "MPa");
Real p_F3 (unit = "MPa");
Real p_F4 (unit = "MPa");
Real p_F5 (unit = "MPa");
Real p_F6 (unit = "MPa");
Real p_F7 (unit = "MPa");
Real p_F8 (unit = "MPa");
Real p_F9 (unit = "MPa");
Real p_F10 (unit = "MPa");
Real p_F11 (unit = "MPa");




equation







end MonliticComplete;
