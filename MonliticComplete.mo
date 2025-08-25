model MonliticComplete

// GLOBAL PARAMETERS //

constant Integer A = 1;
constant Integer B = 2;
constant Integer C = 3;
constant Integer D = 4;
constant Integer E = 5;
constant Integer F = 6;
constant Integer G = 7;
constant Integer H = 8;
constant Integer nComp = 8 "Number of components";
constant Integer VaporComponent[:] = {1,2,3};
constant Integer LiquidComponent[:] = {4,5,6,7,8};
constant Integer AllComponent[:] = 1:8;
constant Integer FeedComponent[:] = 1:6;
constant Integer nChemReact = 3 "Number of chemical reactions";
constant Real RkJ(unit="J/(mol.K)") = 8.31451;
constant Real Rkcal( unit="cal/(mol.K)") = 1.987;

parameter Real Cp_vap[nComp](each unit="kJ/(kg*degC)") = {14.6, 2.04, 1.05, 1.85, 1.87, 2.02, 0.712, 0.628} "Vapor heat capacity";
parameter Real Cp_liq[nComp](each unit="kJ/(kg.degC)") = {0, 0, 0, 7.66, 4.17, 4.45, 2.55, 2.45} "Liquid heat capacity";
parameter Real M[nComp](each unit="g/mol") = {2, 25.4, 28.0, 32.0, 46.0, 48.0, 62.0, 76.0} "Molecular weight";
parameter Real v[nComp,nChemReact] = [-1, -1, -1/3; 0,  0,  0; -1, -1,  0; -1,  0,  -1; 0, -1, -1/3; 0,  0,  1; 1,  0,  0; 0,  1,  0] "Stoichiometric coefficients";
parameter Real Tref( unit="K") = 273.15;
parameter Real cp_cw( unit="kJ/(kg.degC)") = 4.18;
parameter Real rho[nComp](each unit="kg/m3") = {0,0,0,299,365,328,612,617} "Density of liquids";
parameter Real H_vap[nComp](each unit="kJ/kg") = {0,0,0,202,372,372,523,486} "Heat of vaporization";
parameter Real A_ant[nComp] = {0,0,0,20.81, 21.24, 21.24, 21.32, 22.10};
parameter Real B_ant[nComp] = {0,0,0,-1444, -2114, -2144, -2748, -3318};
parameter Real C_ant[nComp] = {0,0,0,259, 265.5, 265.5, 232.9, 249.6};


parameter Real alpha[nChemReact](each unit="1") = {1.0399157,1.010373129,1} "Adjustable parameter of reaction rate";
parameter Real k_comp= 0.7166374645; //Compressor coefficient


// MIXING UNIT //

parameter Real Vmix(unit="m3") = 141.53 "Mixing zone volume";
Real N_mix[nComp](each unit="kmol",
              start={48.83796012,13.49310238,40.03019459,10.44186782,
              28.48830431,2.514120166,5.403006586,2.517154711},
              each fixed=true)  "Total molar holdup";
Real pMix(unit="MPa");
Real Tmix(unit="K", start=359.25, fixed=true);
Real y_Mix[nComp];

// VALVE 1 //

parameter Real Cv1(unit = "kmol/s/MPa**(1/2)") = 0.83337;

// REACTOR //

parameter Real rho_liq_reactor( unit = "kmol/m3")= 9.337145754 "Reactor Liquid density";
parameter Real m_CWS_rea( unit="mol/s") = 93.37/3600*1000/18;
parameter Real T_CWSs_rea_in( unit="K") = 308;
parameter Real T_CWSs_rea_out( unit="K") = 367.599;
parameter Real gamma_rea[nComp] = {0,0,0,0.996011,1,1.078,0.999,0.999} "Activity coefficient reactor";
parameter Real Vrea( unit="m3") = 36.8 "Reactor volume";
Real P_rea_tot(unit="MPa");
Real N_reactor[nComp](each unit="kmol",
              start={5.2046,2.228951,4.65083,0.11534,7.44601,1.14649,56.0529,59.82564},
              each fixed=true)  "Total molar holdup";
Real Rate[nChemReact](each unit="kmol/s") "Reaction Rate";
Real T_rea( unit="K", start=393.55, fixed=true) "Reactor temperature";
Real Q_rea "Cooling heat flow";
Real delta_Hr[nChemReact];
Real Vliq_rea( unit="m3") "Liquid volume";
Real Vvap_rea( unit="m3") "Vapor volume";
Real p_rea[nComp](each unit="MPa") "Reactor partial pressure";
Real x_rea[nComp](each unit="1") "Molar fraction of liquid in the reactor";
Real y_rea[nComp](each unit="1") "Molar fraction of vapor in the reactor";
Real pSat_rea[nComp](each unit="MPa") "Partial pressure at reactor temperature";
Real deltaHrSum "Reaction heat";

// VALVE 2 //

parameter Real Cv2(unit = "kmol/s/MPa**(1/2)") = 1.5344;

// SEPARATOR //

parameter Real rho_liq_separator( unit = "kmol/m3")= 10.29397546 "Reactor Liquid density";
parameter Real gamma_sep[nComp] = {0,0,0,1.001383,1.001383,1.091383,1.001383,0.992188}"Activity coefficient separator";
parameter Real m_CWS_sep( unit="mol/s") = 49.37/3600*1000/18;
parameter Real T_CWSs_sep_in( unit="K") = 313;
parameter Real T_CWSs_sep_out( unit="K") = 350.447;
parameter Real Vsep(unit="m3")= 99.1 "Separator volume";
Real N_separator[nComp](each unit="kmol",
              start={27.497976,12.09634,24.57206,0.0836049,
              5.864057,0.9011249,24.137664,19.752554})"Total molar holdup";
Real T_sep( unit="K", start=353.25, fixed=true) "Separator temperature";
Real Q_sep "Cooling heat flow";
Real HoVsep "Vaporization heat";
Real Vliq_sep( unit="m3") "Liquid volume";
Real Vvap_sep( unit="m3") "Vapor volume";
Real P_sep_tot(unit="MPa");
Real x_sep[nComp](each unit="1") "Molar fraction of liquid in the separator";
Real y_sep[nComp](each unit="1") "Molar fraction of vapor in the separator";
Real pSat_sep[nComp](each unit="MPa") "Partial pressure at separator temperature";
Real p_sep[nComp](each unit="MPa") "Separator partial pressure";



// STRIPPER //

parameter Real PC_1[nComp] = {0,0,0,3.51e-7, 3.69e-7, 3.69e-7, 1.84e-9, 8.32e-8};
parameter Real PC_2[nComp] = {0,0,0,-0.000111, -0.0001, -0.0001, 1.37e-5, -7.64e-6};
parameter Real PC_3[nComp] = {0,0,0,0.011351, 0.010197, 0.010049, 0.000217, 0.000976};
parameter Real PC_4[nComp] = {0,0,0,0.548012, 0.620794, 0.628854, 0.001393, -0.01568};
parameter Real m_CWS_str(unit="mol/s") = 230.31/3600;
parameter Real rho_liq_stripper( unit = "kmol/m3")= 8.6496 "Stripper Liquid density";
Real NG(unit="kmol", start = 20.38135);
Real NH( unit = "kmol", start = 17.313);
Real phi[nComp]"Split factor stripper";
Real T_str( unit="K", start=338.85, fixed=true) "Stripper temperature";
Real Q_str "Cooling heat flow";
Real HoVstr "Vaporization heat";
Real Vliq_str( unit="m3") "Liquid volume";
Real x_str[nComp](each unit="1") "Molar fraction of liquid in the stripper";
Real y_str[nComp](each unit="1") "Molar fraction of vapor in the stripper";
Real pStr(unit="MPa");

// FLUXES PARAMETERS //

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
//Real p_F6 (unit = "MPa");
//Real p_F7 (unit = "MPa");
Real p_F8 (unit = "MPa");
Real p_F9 (unit = "MPa");
Real p_F10 (unit = "MPa");
Real p_F11 (unit = "MPa");

equation

// BOUNDARY CONDITIONS //

F1=11.2/3600;
F2=114.5/3600;
F3=98.0/3600;
F4=417.5/3600;
F8= 1200.5/3600;
F9=15.1/3600;
F10= 259.5/3600;
F11=211.3/3600;

y_F1={0.9999,0.0001,0,0,0,0,0,0};
y_F2={0,0.0001,0,0.9999,0,0,0,0};
y_F3={0,0,0,0,0.9999,0.0001,0,0};
y_F4={0.485,0.005,0.51,0,0,0,0,0};

p_F1=2.7;
p_F2=2.7;
p_F3=2.7;
p_F4=2.7;

T_F1=45+273.15;
T_F2=45+273.15;
T_F3=45+273.15;
T_F4=45+273.15;


// MIXING UNIT //

for i in AllComponent loop
    der(N_mix[i]) = F1*y_F1[i] + F2*y_F2[i] + F3*y_F3[i] + F5*y_F5[i] + F8*y_F8[i] - F6*y_F6[i];
end for;
  
for i in AllComponent loop
    y_Mix[i] = N_mix[i]/sum(N_mix);
end for;
  
pMix*Vmix = 1e-3*sum(N_mix)*RkJ*Tmix;

sum(N_mix[i]*Cp_vap[i]*M[i] for i in 1:nComp)*der(Tmix) =
  F1*sum(y_F1[i]*Cp_vap[i]*M[i] for i in 1:nComp)*(T_F1-Tmix) +
  F2*sum(y_F2[i]*Cp_vap[i]*M[i] for i in 1:nComp)*(T_F2-Tmix) +
  F3*sum(y_F3[i]*Cp_vap[i]*M[i] for i in 1:nComp)*(T_F3-Tmix) +
  F5*sum(y_F5[i]*Cp_vap[i]*M[i] for i in 1:nComp)*(T_F5-Tmix) +
  F8*sum(y_F8[i]*Cp_vap[i]*M[i] for i in 1:nComp)*(T_F8-Tmix);

p_F5=pMix;
//p_F6=pMix;
p_F8=pMix;

y_F6 = y_Mix;

T_F6=Tmix;

// VALVE 1 //

F6 = Cv1*sqrt(pMix-P_rea_tot);

// REACTOR //

for i in 1:nComp loop
    der(N_reactor[i]) = F6*y_F6[i] - F7*y_F7[i] + sum(v[i,j]*Rate[j] for j in 1:nChemReact);
end for;

(sum(N_reactor[i]*Cp_vap[i]*M[i] for i in VaporComponent)+
   sum(N_reactor[i]*Cp_liq[i]*M[i] for i in LiquidComponent))*der(T_rea) =
   1e-3*F6*sum(y_F6[i]*Cp_vap[i]*M[i] for i in FeedComponent)*(T_F6-T_rea)
   - Q_rea - deltaHrSum;

deltaHrSum=sum(Rate[i]*delta_Hr[i] for i in 1:nChemReact);
delta_Hr[1]=1e-3*(sum(v[i,1]*Cp_vap[i]*M[i]*(T_rea-Tref) for i in AllComponent) - 136033.04);
delta_Hr[2]=1e-3*(sum(v[i,2]*Cp_vap[i]*M[i]*(T_rea-Tref) for i in AllComponent) - 93337.9616);
delta_Hr[3]=1e-3*(sum(v[i,3]*Cp_vap[i]*M[i]*(T_rea-Tref) for i in AllComponent));

Rate[1] = if noEvent(p_rea[A]>0 and p_rea[C]>0 and p_rea[D]>0) then
         (1/3600)*alpha[1]*Vvap_rea*exp(44.06-42600/(Rkcal*T_rea))*
         (p_rea[A]*1e3)^1.08*(p_rea[C]*1e3)^0.311*(p_rea[D]*1e3)^0.874
         else 0;
Rate[2] = if noEvent(p_rea[A]>0 and p_rea[C]>0 and p_rea[E]>0) then
         (1/3600)*alpha[2]*Vvap_rea*exp(10.27-19500/(Rkcal*T_rea))*
         (p_rea[A]*1e3)^1.15*(p_rea[C]*1e3)^0.370*(p_rea[E]*1e3)^1.00
         else 0;
Rate[3] = if noEvent(p_rea[A]>0 and p_rea[D]>0 and p_rea[E]>0) then
         (1/3600)*alpha[3]*Vvap_rea*exp(59.50-59500/(Rkcal*T_rea))*
         (p_rea[A]*1e3)*(0.77*(p_rea[D]*1e3)+(p_rea[E]*1e3))
         else 0;
         
Q_rea = 1e-3*m_CWS_rea*cp_cw*(T_CWSs_rea_out-T_CWSs_rea_in);

P_rea_tot = sum(p_rea);

for i in VaporComponent loop
    p_rea[i] = 1e-3*N_reactor[i]*RkJ*T_rea/Vvap_rea;
end for;
for i in LiquidComponent loop
    p_rea[i] = gamma_rea[i]*x_rea[i]*pSat_rea[i];
end for;

for i in VaporComponent loop
     pSat_rea[i] = 0;
end for;
for i in LiquidComponent loop
     pSat_rea[i] = 1e-6*exp(A_ant[i]+B_ant[i]/(T_rea-Tref+C_ant[i]));
end for;

Vliq_rea = sum(N_reactor[i]/rho_liq_reactor for i in LiquidComponent);
Vrea = Vliq_rea + Vvap_rea;

for i in VaporComponent loop
    x_rea[i] = 0;
end for;
for i in LiquidComponent loop
    x_rea[i] = N_reactor[i]/sum(N_reactor[i] for i in LiquidComponent);
end for;

for i in AllComponent loop
    y_rea[i] = p_rea[i]/P_rea_tot;
end for;

y_F7 = y_rea;

T_F7 = T_rea;

// VALVE 2 //

F7 = Cv2*sqrt(P_rea_tot-P_sep_tot);

// SEPARATOR //

for i in AllComponent loop
    der(N_separator[i]) = y_F7[i]*F7 - y_F8[i]*F8 - y_F9[i]*F9 - x_F10[i]*F10;
end for;

(sum(N_separator[i]*Cp_vap[i]*M[i] for i in VaporComponent)+
   sum(N_separator[i]*Cp_liq[i]*M[i] for i in LiquidComponent))*der(T_sep) =
   1e-3*F7*sum(y_F7[i]*Cp_vap[i]*M[i] for i in AllComponent)*(T_F7-T_sep)
   - Q_sep + HoVsep;
   
HoVsep =  1e-3*F10*(sum(x_sep[i]*M[i]*H_vap[i] for i in LiquidComponent));

Q_sep = 1e-3*m_CWS_sep*cp_cw*(T_CWSs_sep_out-T_CWSs_sep_in);

T_F8=T_sep*(p_F8/P_sep_tot)^((1-k_comp)/k_comp);

P_sep_tot = sum(p_sep);

for i in VaporComponent loop
    p_sep[i] = 1e-3*N_separator[i]*RkJ*T_sep/Vvap_sep;
end for;
for i in LiquidComponent loop
    p_sep[i] = gamma_sep[i]*x_sep[i]*pSat_sep[i];
end for;

for i in VaporComponent loop
     pSat_sep[i] = 0;
  end for;
for i in LiquidComponent loop
     pSat_sep[i] = 1e-6*exp(A_ant[i]+B_ant[i]/(T_sep-Tref+C_ant[i]));
end for;

Vliq_sep = sum(N_separator[i]/rho_liq_separator for i in LiquidComponent);

Vsep = Vliq_sep + Vvap_sep;

for i in VaporComponent loop
    x_sep[i] = 0;
end for;
for i in LiquidComponent loop
    x_sep[i] = N_separator[i]/sum(N_separator[i] for i in LiquidComponent);
end for;

for i in AllComponent loop
    y_sep[i] = p_sep[i]/P_sep_tot;
end for;

x_F10=x_sep;

y_F8=y_sep;
y_F9=y_sep;

T_F9=T_sep;
T_F10=T_sep; 
p_F9 = P_sep_tot;
p_F10 = P_sep_tot;


// STRIPPER //

der(NG) = (1-phi[7])*(x_F10[7]*F10 + y_F4[7]*F4) - x_F11[7]*F11;
der(NH) = (1-phi[8])*(x_F10[8]*F10 + y_F4[8]*F4) - x_F11[8]*F11;

(NG*Cp_liq[7]*M[7]+NH*Cp_liq[8]*M[8])*der(T_str) =
   1e-3*F10*sum(x_F10[i]*Cp_vap[i]*M[i] for i in AllComponent)*(T_F10-T_str)+
   1e-3*F4*sum(y_F4[i]*Cp_vap[i]*M[i] for i in AllComponent)*(T_F4-T_str)   + Q_str - HoVstr;

HoVstr=1e-3*(sum( (y_F5[i]*F5 - y_F4[i]*F4)*M[i]*H_vap[i] for i in LiquidComponent));
Q_str = 2.258717*m_CWS_str;

Vliq_str = (NG+NH)/rho_liq_stripper;

for i in VaporComponent loop
    phi[i] = 1;
end for;
for i in LiquidComponent loop
    phi[i] = PC_1[i]*(T_str-Tref)^3+PC_2[i]*(T_str-Tref)^2+PC_3[i]*(T_str-Tref)+PC_4[i];
end for;

F5 = F10 + F4 - F11 - (der(NG) + der(NH));

for i in AllComponent loop
    y_str[i] = phi[i]*(y_F4[i]*F4+ x_F10[i]*F10)/F5;
end for;

for i in VaporComponent loop
    x_str[i] = 0;
end for;

for i in 4:6 loop
    x_str[i] = (y_F4[i]*F4 + x_F10[i]*F10 - y_F5[i]*F5)/F11;
end for;

x_str[7] = ((1- x_F11[4]-x_F11[5]-x_F11[6])*NG)/(NG + NH);
x_str[8] = ((1- x_F11[4]-x_F11[5]-x_F11[6])*NH)/(NG + NH);

y_F5=y_str;
x_F11=x_str;

T_F11=T_str;
T_F5=T_str;

p_F11=pStr;


pStr=pMix; //hypothesis for semplification


end MonliticComplete;
