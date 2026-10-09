clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. 
%  For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200;Tamb=300;P3overP2=7;Pamb=100*kPa;mfurate=0.58*kg/s;AF=170.35;             % These are the ones from the book
cFuel='H2';                                                           % Pick Gasoline as the fuel (other choices check Sp.Name)
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;

%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the Bisection method
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using bisection (https://en.wikipedia.org/wiki/Bisection_method)
cMethod = 'Bisection Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
TL = T1;
TH = 1000;                                                                  % A guess for the TH (must be too high)
iter = 0;
while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;
    for i=1:NSp
        hi2(i)    = HNasa(Ti,SpS(i));
    end
    h2i = Yair*hi2';                                                        % Single value (1x5 times 5x1). Intermediate value
    if h2i > h2
        TH = Ti; % new right boundary
    else
        TL = Ti; % new left boundary
    end
end
T2 = (TH+TL)/2;
T2bis = T2;
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total entropy stage 1
S2  = s2thermal - Rg*log(P2/Pref);                                          % Total entropy stage 2
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
%% Difference between two approaches: so close but not identical
fprintf('----------------------------------------------\n%8s| %9.4f %9.4f  [K]\n----------------------------------------------\n','T2-int vs T2-bis',T2int,T2bis);
%% Here starts your part (compressor,combustor,turbine and nozzle). ...

%bisection [2-3] compressor, isentropic :
cMethod = 'Bisection Method';
sPart = 'Compressor';

P3 = P2*P3overP2;

s3thermal = s2thermal + Rg*log(P3/P2);          %isentropic compression (no change in entropy)

% Bisection to determine T3
TL = T2;
TH = 1000;
iter = 0;

while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;

    for i=1:NSp
        si3(i) = SNasa(Ti,SpS(i));
    end

    s3i = Yair*si3';

    if s3i > s3thermal
        TH = Ti;
    else
        TL = Ti;
    end
end

T3 = (TH+TL)/2;

for i=1:NSp
    hi3(i) = HNasa(T3,SpS(i));
    si3(i) = SNasa(T3,SpS(i));
end

h3 = Yair*hi3';
s3thermal = Yair*si3';

S3 = s3thermal - Rg*log(P3/Pref);

v3 = v2;


wc = h3-h2;                 %compressor work

%% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T2,T3);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P2/kPa,P3/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v2,v3);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h2/kJ,h3/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S2/kJ,S3/kJ);
fprintf('-------------------------------------\n');
fprintf('Compressor work: %9.2f [kJ/kg]\n',wc/kJ);



%interpolation [2-3] compressor, isentropic :
cMethod = 'Interpolation Method';
sPart = 'Compressor';

P3 = P2*P3overP2;           %pressure at exit of compressor

s3thermal = s2thermal + Rg*log(P3/P2);

T3 = interp1(sair_a,TR,s3thermal);

for i=1:NSp
    hi3(i) = HNasa(T3,SpS(i));
    si3(i) = SNasa(T3,SpS(i));
end

h3 = Yair*hi3';
s3thermal = Yair*si3';

S3 = s3thermal - Rg*log(P3/Pref);

v3=v2;

wc = h3-h2;         %compressor work input


fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T2,T3);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P2/kPa,P3/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v2,v3);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h2/kJ,h3/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S2/kJ,S3/kJ);
fprintf('-------------------------------------\n');
fprintf('Compressor work: %9.2f [kJ/kg]\n',wc/kJ);



%Combustion--------------------------------------------------------------------

MO2 = Mi(2);
MN2 = Mi(5);
MH2 = Mi(1); 
MH2O = Mi(4);              %Molar masses of air components
N2overO2 = 0.79/0.21;                    %Proportion of N2 to O2 in the air
nAF = AF/((MO2*0.21 + MN2*0.79)/MH2); %Molar AF ratio

%Mixture Compositions:
%a*H2 + b*O2 + c*N2 --> d*H2O + e*O2 + f*N2
%Only H2 disappears and is converted to H2O, taking 1/2 of one O2 molecule,
%N2 stays untouched

a = 1;
b = 0.21*nAF;
c = 0.79*nAF;
d = 1;
e = b - 1/2;
f = c;

%Mass fractions (A,B,C,D,E,F):

mtot1 = a*MH2 + b*MO2 + c*MN2;       %Total mass of reactants
A = a*MH2/mtot1;
B = b*MO2/mtot1;
C = c*MN2/mtot1;
Reactants_fractions = [A;B;C;0];

mtot2 = d*MH2O + e*MO2 + f*MN2;    %Total mass of products
D = d*MH2O/mtot2;
E = e*MO2/mtot2;
F = f*MN2/mtot2;

Products_fractions = [0;E;F;D];

%Composition table

Components = ["H2";"O2";"N2";"H2O"];
Mixture_compositions = table(Components,Reactants_fractions,Products_fractions);

%Equivalence ratio:

nAFstoec = 0.5/0.21;              %Stoechiometric molar Air/fuel ratio
AFstoec = nAFstoec*((0.21*MO2+0.79*MN2)/MH2); %Stoechiometric Air/Fuel ratio (mass)
Equivalence_ratio = AFstoec/AF;

%Specific gas constant:
Yprod = [0, E, 0, D, F];
Yair = [0 0.21 0 0 0.79];

Rg3 = Runiv*sum(Yair./Mi);
Rg4 = Runiv*sum(Yprod./Mi);
%Thermodynamics part :

mairrate = AF*mfurate;
Hfuel = HNasa(Tamb,SpS(1)); %[J/kg]

h4 = (mairrate.*h3 + mfurate.*Hfuel)/(mairrate + mfurate); %Enthalpy as a function of mass flux and initial enthalpies 
fcth = @(T) D*HNasa(T,SpS(4)) + E*HNasa(T,SpS(2)) + F*HNasa(T,SpS(5)) - h4;
if T3 < 1000    %Making sure Tmax never exceeds the Nasa polynomials range
    Tmax = T3 + 2000;
else
    Tmax = 3000;
end
intT = [T3, Tmax];                      
T4 = fzero(fcth,intT);




%----------------------------------------------------------------------------------------
%% Turbine [4-5]: adiabatic and isentropic

% Product mass fractions, matching the order in SpS
Yprod = [0 E 0 D F];
Rgprod = Runiv * sum(Yprod ./ Mi);

% Constant-pressure combustor
P4 = P3;

% Turbine power equals compressor power
mgasrate = mairrate + mfurate;
turbinePower = mairrate * wc;             % [W]
h5 = h4 - turbinePower / mgasrate;        % [J/kg gas]

% Find T5 using the combustor's existing enthalpy function
% fcth(T) = product enthalpy at T minus h4
fcth5 = @(T) fcth(T) + h4 - h5;
T5 = fzero(fcth5, [TR(1), T4]);

% Temperature-dependent entropy at states 4 and 5
for i = 1:NSp
    si4(i) = SNasa(T4, SpS(i));
    si5(i) = SNasa(T5, SpS(i));
    hi5(i) = HNasa(T5, SpS(i));
end

s4thermal = Yprod * si4';
s5thermal = Yprod * si5';

% Isentropic expansion: S5 = S4
P5 = P4 * exp((s5thermal - s4thermal) / Rgprod);

% Mole fractions to Mass fractions Entropy of mixing; 
Xprod = (Yprod ./ Mi) / sum(Yprod ./ Mi);
present = Yprod > 0;
smix = -Runiv * sum((Yprod(present) ./ Mi(present)) ...
    .* log(Xprod(present)));

% Total specific entropy [J/kg/K]
S4 = s4thermal - Rgprod * log(P4/Pref) + smix;
S5 = s5thermal - Rgprod * log(P5/Pref) + smix;

% Negligible velocity at turbine outlet
v5 = 0;

% Check the enthalpy obtained from the calculated temperature
h5check = Yprod * hi5';
powerError = mgasrate * (h4 - h5check) - turbinePower;

% Results
fprintf('\nTurbine [4-5]\n');
fprintf('T5:             %.2f K\n', T5);
fprintf('P5:             %.2f kPa\n', P5/kPa);
fprintf('v5:             %.2f m/s\n', v5);
fprintf('h5:             %.2f kJ/kg\n', h5/kJ);
fprintf('Turbine power:  %.2f kW\n', turbinePower/kJ);
fprintf('S4:             %.6f kJ/kg/K\n', S4/kJ);
fprintf('S5:             %.6f kJ/kg/K\n', S5/kJ);
fprintf('Power residual: %.6f W\n', powerError);
