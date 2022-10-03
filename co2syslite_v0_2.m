function output = ...
    co2syslite_v0_2(tempIn,SalIn,pressIn,PO4In,SIn,TAlkIn,DICin,Ca,Mg,SO4)
% input values
temp = tempIn;
pressure = pressIn;
salinity = SalIn;
PO4 = PO4In;
SI = SIn;
CaComp = Ca;
MgComp = Mg;
SO4Comp = SO4;

% initial carbon system conditions
TAlk = TAlkIn/1e6;
DIC = DICin/1e6;


constants = CalculateConstants(temp,pressure,salinity,PO4,SI,CaComp,MgComp,SO4Comp);
pH = CalculatepHfromTATC(TAlk,DIC,constants);
pCO2 = CalculatepCO2fromTCpH(DIC,pH,constants);
[OmegaCa,CO3] = CaSolubility(temp,pressure,salinity,DIC,pH,constants,CaComp,MgComp,SO4Comp);

output = [pH pCO2 OmegaCa CO3.*1e6];

end


% shit hot stripped down co2sys (~200 times faster than passing to the function)
function constants = CalculateConstants(temp,pressure,salinity,PO4,SI,CaComp,MgComp,SO4Comp)
TempC       = temp;
Pdbar       = pressure;
Sal          = salinity;
sqrSal       = sqrt(salinity);
TP           = PO4;
TSi          = SI;
RGasConstant = 83.1451; 
IonS = 19.924 .* Sal ./ (1000 - 1.005   .* Sal); % ionic strength

% calculate values of constants
TempK    = TempC + 273.15;
RT       = RGasConstant.*TempK;
logTempK = log(TempK);
Pbar     = Pdbar ./ 10;

TB =  0.0004157.*Sal./35; % in mol/kg-SW   total boron
TF = (0.000067./18.998).*(Sal./1.80655); % total F
TS = (0.14./96.062).*(Sal./1.80655); % total S

% K0
TempK100  = TempK./100;
lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + Sal .*...
    (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
K0   = exp(lnK0);

% calculate KS
lnKS = -4276.1./TempK + 141.328 - 23.093.*logTempK +...             
      (-13856./TempK + 324.57 - 47.986.*logTempK).*sqrt(IonS) +...     
      (35474./TempK - 771.54 + 114.723.*logTempK).*IonS +...           
      (-2698./TempK).*sqrt(IonS).*IonS + (1776./TempK).*IonS.^2; 
KS = exp(lnKS)...            % this is on the free pH scale in mol/kg-H2O
    .* (1 - 0.001005 .* Sal);   % convert to mol/kg-SW

% KF
lnKF = 1590.2./TempK - 12.641 + 1.525.*IonS.^0.5;
KF   = exp(lnKF)...                 % this is on the free pH scale in mol/kg-H2O
    .*(1 - 0.001005.*Sal);          % convert to mol/kg-SW

SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);    % sws to total
FREEtoTOT =  1 + TS./KS;    % convert free to total scale

% fH
fH = 1.2948 - 0.002036.*TempK + (0.0004607 -...
        0.000001475.*TempK).*Sal.^2;

% KB
lnKBtop = -8966.9 - 2890.53.*sqrSal - 77.942.*Sal +...
        1.728.*sqrSal.*Sal - 0.0996.*Sal.^2;
lnKB = lnKBtop./TempK + 148.0248 + 137.1942.*sqrSal +...
    1.62142.*Sal + (-24.4344 - 25.085.*sqrSal - 0.2474.*...
    Sal).*logTempK + 0.053105.*sqrSal.*TempK;
KB = exp(lnKB)...    % this is on the total pH scale in mol/kg-SW
    ./SWStoTOT;         % convert to SWS pH scale

% KW
lnKW = 148.9802 - 13847.26./TempK - 23.6521.*logTempK +...
    (-5.977 + 118.67./TempK + 1.0495.*logTempK).*...
    sqrSal - 0.01615.*Sal;
KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2

% KPs, KSi
lnKP1 = -4576.752./TempK + 115.54 - 18.453.*logTempK + (-106.736./TempK +...
    0.69171).*sqrSal + (-0.65643./TempK - 0.01844).*Sal;
KP1 = exp(lnKP1);
lnKP2 = -8814.715./TempK + 172.1033 - 27.927.*logTempK + (-160.34./TempK +...
    1.3566).*sqrSal + (0.37335./TempK - 0.05778).*Sal;
KP2 = exp(lnKP2);
lnKP3 = -3070.75./TempK - 18.126 + (17.27039./TempK + 2.81197).*sqrSal +...
    (-44.99486./TempK - 0.09984).*Sal;
KP3 = exp(lnKP3);
lnKSi = -8904.2./TempK + 117.4 - 19.334.*logTempK + (-458.79./TempK +...
    3.5913).*sqrt(IonS) + (188.74./TempK - 1.5998).*IonS +...
    (-12.1652./TempK + 0.07871).*IonS.^2;
KSi = exp(lnKSi)...                % this is on the SWS pH scale in mol/kg-H2O
    .*(1 - 0.001005.*Sal);        % convert to mol/kg-SW

    
% K1, K2
pK1 = 3670.7./TempK - 62.008 + 9.7944.*logTempK...
         - 0.0118.*Sal + 0.000116.*Sal.^2;
K1 = 10.^(-pK1); % this is on the SWS pH scale in mol/kg-SW
% This is in Table 4 on p. 1739.
pK2 = 1394.7./TempK + 4.777 - 0.0184.*Sal + 0.000118.*Sal.^2;
K2 = 10.^(-pK2); % this is on the SWS pH scale in mol/kg-SW



%%%% pressure corrections
deltaV  = -25.5 + 0.1271.*TempC;  % K1
Kappa   = (-3.08 + 0.0877.*TempC)./1000;
lnK1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV  = -15.82 - 0.0219.*TempC; % K2
Kappa   = (1.13 - 0.1475.*TempC)./1000;
lnK2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV  = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2; % KB
Kappa   = -2.84./1000; 
lnKBfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV  = -20.02 + 0.1119.*TempC - 0.001409.*TempC.^2; % KW
Kappa   = (-5.13 + 0.0794.*TempC)./1000; 
lnKWfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV = -9.78 - 0.009.*TempC - 0.000942.*TempC.^2;     % KF
Kappa = (-3.91 + 0.054.*TempC)./1000;
lnKFfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV = -18.03 + 0.0466.*TempC + 0.000316.*TempC.^2;   % kS
Kappa = (-4.53 + 0.09.*TempC)./1000;
lnKSfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV = -14.51 + 0.1211.*TempC - 0.000321.*TempC.^2;   % KP1
Kappa  = (-2.67 + 0.0427.*TempC)./1000;
lnKP1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV = -23.12 + 0.1758.*TempC - 0.002647.*TempC.^2;   % KP2
Kappa  = (-5.15 + 0.09  .*TempC)./1000;
lnKP2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV = -26.57 + 0.202 .*TempC - 0.003042.*TempC.^2;   % KP3
Kappa  = (-4.08 + 0.0714.*TempC)./1000;
lnKP3fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

deltaV = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2;   % KSi
Kappa  = -2.84./1000;
lnKSifac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

% CorrectKsForPressureHere:
K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
KWfac  = exp(lnKWfac);  KW  = KW .*KWfac;
KBfac  = exp(lnKBfac);  KB  = KB .*KBfac;
KFfac  = exp(lnKFfac);  KF  = KF .*KFfac;
KSfac  = exp(lnKSfac);  KS  = KS .*KSfac;
KP1fac = exp(lnKP1fac); KP1 = KP1.*KP1fac;
KP2fac = exp(lnKP2fac); KP2 = KP2.*KP2fac;
KP3fac = exp(lnKP3fac); KP3 = KP3.*KP3fac;
KSifac = exp(lnKSifac); KSi = KSi.*KSifac;

SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;
pHfactor = SWStoTOT;

% ConvertFromSWSpHScaleToChosenScale:
K1  = K1.* pHfactor; K2  = K2.* pHfactor;
KW  = KW.* pHfactor; KB  = KB.* pHfactor;
KP1 = KP1.*pHfactor; KP2 = KP2.*pHfactor;
KP3 = KP3.*pHfactor; KSi = KSi.*pHfactor;

Delta = (57.7 - 0.118.*TempK);
b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
% For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
P1atm = 1.01325; % in bar
FugFac = exp((b + 2.*Delta).*P1atm./RT);

VPWP = exp(24.4543 - 67.4509.*(100./TempK) - 4.8489.*log(TempK./100));
VPCorrWP = exp(-0.000544.*Sal);
VPSWWP = VPWP.*VPCorrWP;
VPFac = 1 - VPSWWP;

% palaeo Ks following Zeebe & Tyrell [2019]
K1rat = 1 + 5/1000.*(CaComp./10.3-1) + 17/1000.*(MgComp./53-1) + ...
    208./1000.*(SO4Comp./(290.46/10.3)-1);
K2rat = 1 + 157/1000.*(CaComp./10.3-1) + 420/1000.*(MgComp./53-1) + ...
    176./1000.*(SO4Comp./(290.46/10.3)-1);

K1 = K1*K1rat;
K2 = K2*K2rat;
    
% for passing to local functions
constants = [sqrSal Pbar RT ...
    K0 fH FugFac TempK logTempK ... 
    K1 K2 KW KB KF KS KP1 KP2 KP3 KSi .... 
    TB TF TS TP TSi VPFac];
end


function [pH] = CalculatepHfromTATC(TAx,TCx,constants)

K1 = constants(9);
K2 = constants(10);
KW = constants(11);
KB = constants(12);
KF = constants(13);
KS = constants(14);
KP1 = constants(15);
KP2 = constants(16);
KP3 = constants(17);
KSi = constants(18);
TB = constants(19);
TF = constants(20);
TS = constants(21);
TP = constants(22);
TSi = constants(23);

pHGuess     = 8;       % this is the first guess
pHTol       = 0.0001;  % tolerance for iterations end
ln10        = log(10); %
pHx         = pHGuess; %
deltapH     = pHTol+1;
while any(abs(deltapH) > pHTol)
    H         = 10.^(-pHx);
    Denom     = (H.*H + K1.*H + K1.*K2);
    CAlk      = TCx.*K1.*(H + 2.*K2)./Denom;
    BAlk      = TB.*KB./(KB + H);
    OH        = KW./H;
    PhosTop   = KP1.*KP2.*H + 2.*KP1.*KP2.*KP3 - H.*H.*H;
    PhosBot   = H.*H.*H + KP1.*H.*H + KP1.*KP2.*H + KP1.*KP2.*KP3;
    PAlk      = TP.*PhosTop./PhosBot;
    SiAlk     = TSi.*KSi./(KSi + H);
    FREEtoTOT = (1 + TS./KS); % pH scale conversion factor
    Hfree     = H./FREEtoTOT; % for H on the total scale
    HSO4      = TS./(1 + KS./Hfree); % since KS is on the free scale
    HF        = TF./(1 + KF./Hfree); % since KF is on the free scale
    Residual  = TAx - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF;
    % find Slope dTA/dpH;
    % (this is not exact, but keeps all important terms);
    Slope     = ln10.*(TCx.*K1.*H.*(H.*H + K1.*K2 + 4.*H.*K2)./Denom./Denom + BAlk.*H./(KB + H) + OH + H);
    deltapH   = Residual./Slope; % this is Newton's method
    % to keep the jump from being too big;
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1;
        deltapH(FF)=deltapH(FF)./2;
    end
    pHx       = pHx + deltapH; % Is on the same scale as K1 and K2 were calculated...
end
pH=pHx;


end


function [pCO2] = CalculatepCO2fromTCpH(TCx,pHx,constants)
K0 = constants(4);
K1 = constants(9);
K2 = constants(10);
FugFac = constants(6);

H            = 10.^(-pHx);
fCO2x        = TCx.*H.*H./(H.*H + K1.*H + K1.*K2)./K0;
fCO2 = fCO2x;

pCO2 = fCO2./FugFac.*1000000;  % Generate pCO2 from fCO2, convert to uatm



end % end nested function


function [OmegaCa,CO3] = CaSolubility(temp,pressure,salinity,TC,pH,constants,CaComp,MgComp,SO4Comp)

TempC = temp;
Pdbar = pressure;
Sal = salinity;
sqrSal = sqrt(salinity);

K1 = constants(9);
K2 = constants(10);
TempK    = TempC + 273.15;
RGasConstant = 83.1451; 
RT       = RGasConstant.*TempK;
logTempK = log(TempK);
Pbar     = Pdbar ./ 10;


% CalculateCa:
% '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
% '       this is .010285.*Sali./35
Ca = 0.02128./40.087.*(Sal./1.80655);% ' in mol/kg-SW
% CalciteSolubility:
% '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
logKCa = -171.9065 - 0.077993.*TempK + 2839.319./TempK;
logKCa = logKCa + 71.595.*logTempK./log(10);
logKCa = logKCa + (-0.77712 + 0.0028426.*TempK + 178.34./TempK).*sqrSal;
logKCa = logKCa - 0.07711.*Sal + 0.0041249.*sqrSal.*Sal;
% '       sd fit = .01 (for Sal part, not part independent of Sal)
KCa = 10.^(logKCa);% ' this is in (mol/kg-SW)^2
% AragoniteSolubility:
% '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
logKAr = -171.945 - 0.077993.*TempK + 2903.293./TempK;
logKAr = logKAr + 71.595.*logTempK./log(10);
logKAr = logKAr + (-0.068393 + 0.0017276.*TempK + 88.135./TempK).*sqrSal;
logKAr = logKAr - 0.10018.*Sal + 0.0059415.*sqrSal.*Sal;
% '       sd fit = .009 (for Sal part, not part independent of Sal)
KAr    = 10.^(logKAr);% ' this is in (mol/kg-SW)^2
% PressureCorrectionForCalcite:
% '       Ingle, Marine Chemistry 3:301-319, 1975
% '       same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
% '       has typos (-.5304, -.3692, and 10^3 for Kappa factor)
deltaVKCa = -48.76 + 0.5304.*TempC;
KappaKCa  = (-11.76 + 0.3692.*TempC)./1000;
lnKCafac  = (-deltaVKCa + 0.5.*KappaKCa.*Pbar).*Pbar./RT;
KCa       = KCa.*exp(lnKCafac);
% PressureCorrectionForAragonite:
% '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
% '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
% '       and 10^3 for Kappa factor)
deltaVKAr = deltaVKCa + 2.8;
KappaKAr  = KappaKCa;
lnKArfac  = (-deltaVKAr + 0.5.*KappaKAr.*Pbar).*Pbar./RT;
KAr       = KAr.*exp(lnKArfac);


Kspcrat = 1 + 185/1000.*(CaComp./10.3-1) + 518/1000.*(MgComp./53-1) + ...
    106./1000.*(SO4Comp./(290.46/10.3)-1);
KCa = KCa*Kspcrat;
KAr = KAr*Kspcrat;

% CalculateOmegasHere:
H = 10.^(-pH);
HCO3 = TC.*K1.*H  ./(K1.*H + H.*H + K1.*K2);
CO3 = TC.*K1.*K2./(K1.*H + H.*H + K1.*K2);
OmegaCa = CO3.*Ca./KCa; % OmegaCa, dimensionless
OmegaAr = CO3.*Ca./KAr; % OmegaAr, dimensionless

end % end nested function