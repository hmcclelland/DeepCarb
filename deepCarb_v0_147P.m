%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% This is DeepCarb version 0.147 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DeepCarb is a simple carbon cycle box model consisting of three boxes
% (atmosphere, surface ocean, deep ocean) and two infinite reservoirs
% (terrestrial silicate rock and deep ocean CaCO3). It is designed to
% assess the impact of carbon cycle perturbations on atmospheric CO2 and
% seawater carbonate chemistry, particularly over i) long timescales and
% ii) in a Monte Carlo fashion with random permutations of key model
% parameterisations.
%
% By default, the model can be run in one of two ways. Setting the variable
% nMC = 1 will result in a simulation run only once using the default set
% of model parameterisations. Setting nMC = n will result in n model
% simulations in which key parameterisations are randomly varied over wide
% ranges. These parameterisations are i) the relationship between the rate 
% of carbon drawdown from silicate weathering and temperature, ii) the
% slope of the relationship between CaCO3 dissolution/preservation to the 
% saturation state of the deep ocean box, and iii) the slope of the
% relationship between CaCO3 production in the surface ocean and the
% saturation state of the surface ocean box.
%
% To maintain flexibility, this script is not a function. Before running
% the model, at least four parameters need to be manually changed in the
% script. These are:
%   1) nMC - the number of simulations to run, as described above.
%   2) parfor/for - use the parallel computing toolbox to speed the model
%      up when running multiple simulations. To do so, manually change
%      for to parfor on line 152.
%   3) ageEnd - the time in years at which the simulation should stop. Note
%      that rather than preallocating all matrices, the model loads a set of
%      stable simulations which includes the necessary variables. These
%      have a length of 15 Myr. If longer simulations are required, the
%      matrices need to be modified in the accompanying data files from 
%      which the initial stable conditions are loaded.
%   4) The model needs to be perturbed in some way. Examples of how this
%      can be done are given below andthroughout the script. For example, 
%      atmospheric CO2 could be modified or prescribed (atCO2in), weathering 
%      rates could be modified (w12), seawater [Ca] could be changed either 
%      directly (oceanCa), or via the Ca/total cation ratio of the
%      weathering product (weatherCa).
%
% In addition, the script produces a plot of the results. Edit linlog to 
% control whether the data are shown on a linear or logarithmic scale 
% (1 = linear, 0  = log).
%
% The results are contained in the cell array dataOut which has dimensions 
% n*1. The format of each cell is an ageEnd*15 matrix, with the following 
% columns:
%   1) Total C in the atmosphere box (Gt)
%   2) Total C in the surface ocean box (Gt)
%   3) TAlk in the surface ocean box (umol/kg)
%   4) DIC in the surface ocean box (umol/kg)
%   5) pH in the surface ocean box (total scale)
%   6) Omega calcite in the surface ocean box at 10.3 mM [Ca]
%   7) TAlk in the deep ocean box (umol/kg)
%   8) DIC in the deep ocean box (umol/kg)
%   9) pH in the deep ocean box (total)
%   10) Omega calcite in the deep ocean box at 10.3 mM [Ca]
%   11) Seawater [Ca] (mM)
%   12) C diffusion between the atmosphere and surface ocean box (Gt/year)
%       Negative values represent net diffusion to the atmosphere
%   13) Net CaCO3 burial/dissolution rate (moles/year)
%       Negative values represent net burial
%   14) Carbon drawdown as a result of silicate weathering (Gt/year)
%   15) Change in CaCO3 export from the surface ocean box (moles/year)
%
%
% This script depends on the following functions which are packaged here
% but are not my work. Please give appropriate credit to the original
% authors.
% - A stripped-down version of co2sys ('co2syslite') is used, which runs
% approximately 1000 times faster than co2sys. Please credit (e.g.) van
% Heuven et al. (2011) CO2SYS v 1.1, MATLAB program developed for CO2 
% system calculations.
% - The function parforwait is required to provide a waitbar within a
% parfor loop. See:
% https://uk.mathworks.com/matlabcentral/fileexchange/71083-waitbar-for-parfor
% - The parallel computing toolbox is required if parfor is used.
% - The random number generator (normrnd) requires the statistics and
% machine learning toolbox
%
%
% ------------------------------- EXAMPLES -------------------------------
%
% Where these example recommend changing the values of nMC, linlog, ageEnd,
% and the use or otherwise of the parallel computing toolbox, these can be
% found on lines 135, 140, 136, and 156 respectively. In the latter case,
% 'for' simply needs to be replaced with 'parfor'.
%
% 1. Simulate the anthropogenic C release and run the model once using the
% default settings. Set nMC = 1, linlog = 0 (to view the long and short
% term effects, ageEnd = 1e6, and use a for loop rather than parfor.
% Uncomment lines 274-275 to load the anthropogenic C release data, and
% comment line 280 specifying no atmospheric CO2 perturbation.
%
% 2. Simulate the PETM and run the model with 100 random permutations of
% the key model parameterisations. Set nMC = 100, linlog = 1, ageEnd =
% 2.5e5, and use the parallel computing toolbox (change the for loop on
% line 156 to a parfor loop). Uncomment line 263 to perturb the model with
% one possible PETM C-release scenario, and comment line 280 specifying no
% atmospheric CO2 perturbation.
%
% 3. Prescribe an atmospheric CO2 doubling over a range of doubling rates.
% Set nMC = 50, linlog = 0, ageEnd = 1e7, and use the parallel computing
% toolbox. Uncomment the code on lines 339-347, which will vary the
% doubling rate from ~10 to 1e6 years. Note that setting nMC>1 will result
% in random selections of the key model parameterisations described above.
% To maintain consistency between the doubling rate experiments, remove
% this feature by uncommenting line 235.
%
% 4. Decrease seawater [Ca] starting from Eocene boundary conditions and
% run a full Monte Carlo simulation. Set nMC = 1000, linlog = 1, ageEnd =
% 15e6, and use the parallel computing toolbox. Uncomment line 293 to drive
% [Ca] down by 10 mM over 8 Myr, starting at model year 1 million, and then
% let the model equilibrate. Comment out line 294 which specifies no
% seawater [Ca] perturbation. Comment out lines 162-163 which load the 
% modern stable conditions, and uncomment lines 165-166 to load Eocene
% boundary conditions.
%
% 5. Drive CO2 down by increasing the weathering rate independently of
% climate, and run a full Monte Carlo simulation. Set nMC = 1000, linlog =
% 1, ageEnd = 15e6, and use the parallel computing toolbox. Uncomment lines
% 383-389 to increase the silicate weathering rate by 5% over 4 Myr,
% starting at model year 1 million, and then hold it at that point. Comment
% out lines 162-163 which load the modern stable conditions, and uncomment
% lines 165-166 to load Eocene boundary conditions.
%
% ------------------------------------------------------------------------

clear vars

% edit key model setup parameters
nMC = 1;   % number of MC simulations. change for to parfor if >~10!
ageEnd = 2.5e5;  % how many model years to run?

% edit figure properties
figN = 1; % figure number
linlog = 1; % edit here, log = 0
if nMC<11    % set  plot transparency
    trans = 1;
else
    trans = 0.97668*nMC^-0.54185;
end

% preallocate cell arrays for output files
dataOut = cell(nMC,1);
randOut = cell(nMC,1);
dataOutStep = cell(1,1);

if nMC>1    % use either parfor waitbar or normal waitbar
    WaitMessage = parfor_wait(nMC, 'ReportInterval', 1,  'Waitbar', true);
end
tic
for j = 1:nMC
    
    % the model is packaged with two sets of initial stable conditions
    % Modern - 280 ppm atmospheric CO2
    % Eocene - 800 ppm atmospheric CO2
    % modern
      loadStable = ...
          load('initial_stable_v0_147_modern_280ppm_10stepsize_200mSurface_cut.mat');
    % Eocene
%     loadStable = ...
%        load('initial_stable_v0_147_Eo_800ppm_10stepsize_200mSurface_cut.mat');

    % load matrices containing stable conditions
    age = loadStable.age;
    sT = loadStable.sT;
    atCO2in = loadStable.atCO2in;
    atCO2out = loadStable.atCO2out;
    co2sysDeep = loadStable.co2sysDeep;
    co2sysOutSurf = loadStable.co2sysOutSurf;
    dm = loadStable.dm;
    m1 = loadStable.m1;
    m2 = loadStable.m2;
    k12 = loadStable.k12;
    oceanCa = loadStable.oceanCa;
    absChangeSurfOmega = loadStable.absChangeSurfOmega;
    
    oceanVol = loadStable.oceanVol;
    surfVol = loadStable.surfVol;
    deepVol = loadStable.deepVol;
    weatherCa = loadStable.weatherCa;
    DCa = loadStable.DCa;
    ocDICin = loadStable.ocDICin;
    ocTAlkin = loadStable.ocTAlkin;
    orgC = loadStable.orgC;
    k = loadStable.k;
    
    TAlkSurfMix = loadStable.TAlkSurfMix;
    DICSurfMix = loadStable.DICSurfMix;
    TAlkDeepOut = loadStable.TAlkDeepOut;
    DICDeepOut = loadStable.DICDeepOut;
    TAlkDeepMix = loadStable.TAlkDeepMix;
    DICDeepMix = loadStable.DICDeepMix;
    TAlkRelDeep = loadStable.TAlkRelDeep;
    DICRelDeep = loadStable.DICRelDeep;
       
    sT(1) = 10;     % not required, included for completeness
    runL = ageEnd/sT(1);     % output matrix size

    % preallocate arrays
    absChangeSurfOmega(3:runL,1) = NaN;
    atCO2in(3:runL,1) = NaN;
    atCO2out(3:runL,1) = NaN;
    co2sysDeep(3:runL,4) = NaN;      %4
    co2sysOutSurf(3:runL,4) = NaN;   %4
    DICDeepOut(3:runL,1) = NaN;
    DICSurfMix(3:runL,1) = NaN;
    dm(3:runL,1) = NaN;
    k12(3:runL,1) = NaN;
    m1(3:runL,1) = NaN;
    m2(3:runL,1) = NaN;
    oceanCa(3:runL,1) = NaN;
    sT(3:runL,1) = NaN;
    TAlkDeepOut(3:runL,1) = NaN;
    TAlkSurfMix(3:runL,1) = NaN;
    age(3:runL,1) = NaN;
    
    % These random numbers set the sensitivities of the key model
    % parameterisations. The first two values control the shape and slope
    % of the SiW-T relationship, the third sets the slope of the
    % relationship between CaCO3 burial/preservation and Omega in the deep
    % box, the fourth sets the slope between the change in CaCO3 production
    % in the surface ocean and Omega
    if nMC==1   % if no MC simulation, use default parameters
        MCrand = [0.2 0.1 1 0];
    else        % otherwise perturb these
        MCrand = ...
            [abs(normrnd(0,0.1,1,1)) rand(1,1)*0.1 rand(1,1)*0.5+0.75 rand(1,1).*10];
        % Uncomment this to run perturbed the model in multiple different
        % ways with the same parameterisation sensitivities
        %MCrand = [0.2 0.1 1 0];
    end
    % This would alternatively give a similar SiW-CO2 slope to 'standard' 
    % LOSCAR(Penman et al. 2021)
    % MCrand = [0.5 4 1 0];
    % This would produce a simulation with a very weak SiW-CO2 sensitivity
    % MCrand = [0.005 0.005 1 0];   
    
    %%% Example CO2 perturbations here %%%
    % perturbations are performed in the loop to avoid broadcast variables
    % CO2 is added in Pg/Gt, where the modern atmosphere
    % has 597 Pg @ 280 ppm. In the model, Pg/Gt is converted to concentration 
    % before the perturbation is applied.
    %
    % CO2add requires three numbers - 1) start year, 2) duration, and 3) 
    % C release per year (Gt). There is no limit to the length of the matrix,
    % facilitating multiple releases.
    
    % If all anthropogenic carbon emissions stopped in 2022. 
    % Total CO2 emissions to 2022 are 1.5 trillion tons, released at a 
    % constant rate here for simplicity
    %CO2add = [200 200 1500*12/44/200];
    % If CO2 emissions stop at 2100 and then all anthropogenic C is removed
    % by 2200
    %CO2add = [200 280 1500*12/44/200 ; ...
    %    480 100 -1500*12/44/100];

    % The PETM, approxiamtely according to Komar & Zeebe (2010)
    %CO2add = [1000+5e4 3000 3000/3000 ; 10000+5e4 40000 1480/40000];

    %%%%%% CO2 perturbation similar to Zeebe 2008 - vary both release time
    %%%%%% and release magnitude
    % Release times from 40 to 1210 years, 600 to 8000 GtC total
    % release magnitudes
    %        CO2add = [200 10*(rem(j-1,10)+2)^2 (800*(floor((j-1)/10)+1)-200)/(10*(rem(j-1,10)+2)^2)];

    % Or load from file (note that in this example file the approximate CO2
    % removal by forests has been subtracted out of the yearly C emissions
    % data)
    %CO2add = load("CO2_release_anthropogenic");
    %CO2add = CO2add.a;

    % No atmospheric CO2 perturbation. Note that this variable is required
    % in the main loop below. If no CO2 perturbation is desired, this line
    % must be uncommented
    CO2add = [0 0 0];
    
    % calculate CO2 to add in each time step based on the above
    CO2in = zeros(ageEnd/sT(1)+2,1);
    for m = 1:size(CO2add,1)
        for l = 1:size(CO2in,1)
            if l*sT(1)>CO2add(m,1) && l*sT(1)<=CO2add(m,1) + CO2add(m,2)
                CO2in(l,1) = CO2in(l,1) + CO2add(m,3)*sT(1);
            end
        end
    end
    
    % Uncomment this to force a change in seawater [Ca]
    % addCa = [1e6 9e6 -10/8e6];  % drive seawater Ca externally
    addCa = [0 1e6 0];        % no external change in seawater Ca
    
    % Additional CaCO3 burial from the surface ocean. Multiplier is used in
    % a parameterisation below.
    Caburial = 1e11;
    % Proportion of carbonate alkalinity in the weathering product
    TAlkFac = 0.8168;
    
    % preallocate arrays that aren't loaded with previous stable
    % configuration
    c12 = NaN(size(m1,1),1);
    w12 = NaN(size(m1,1),1);


  if nMC==1
        w = waitbar(0,'running...');
  end
    i = 2;
    % save current time and break while loop if iteration is too long
    time0 = tic;
    % should be quicker than 20s per Myr. If not, stop the simulation
    % (usually required if the model drifts into near-impossible carbonate
    % chemistry space)
    timeLimit = ageEnd/1e6*20;
    % stop the simulation if 1) the end has been reached, 2) it is taking
    % too long, or 3) if the results have gotten a bit unreal
    while age(i)<ageEnd && i<=runL && toc(time0)<timeLimit && isreal(m1(i-1))
        i = i+1;

        % time increment in this step. Not required - only included in case
        % a variable or different timestep is desired in future
        sT(i) = 10;    
        age(i) = age(i-1) + sT(i);  % age progression in this iteration
        
        % ocean mixing time = 2 kyr by default
        k23 = oceanVol/4000*2*sT(i);   % transfer rate between surface & deep ocean boxes for Tr = 4ka (m3/year)  
        
        % for responsive atmospheric CO2 as prescribed above
        % Convert GtC to CO2 using a linear conversionfactor of 400/597 for
        % simplicity.
        atCO2in(i) = atCO2out(i-1) + CO2in(i)*280/597;  % CO2 at start of this step
 
        % something like this could be used to prescribed an atmospheric CO2
        % doubling over a range of timescales
        atCO2add = 0;
%         atCO2add = [100 10^(1+j/10) 280];
%         if age(i)>atCO2add(1) && age(i)<(atCO2add(1)+atCO2add(2))
%             atCO2in(i) = 280 + atCO2add(3)/(atCO2add(2))*...
%                 (age(i)-atCO2add(1));            
%         elseif age(i)>=(atCO2add(1)+atCO2add(2))
%             atCO2in(i) = 560;
%         else
%             atCO2in(i) = 280;
%         end
        
        % surface ocean-atmosphere carbon transfer rate
        % Broecker & Peng transfer rate = 0.06 following Zeebe 2012 
        % multiplied by surface area and CO2 disequilibrium
        % result is in moles/m2/year
        k12(i) = 0.06*360000000*1000^2*(atCO2in(i) - co2sysOutSurf(i-1,2));

        % The mass of carbon in the atmosphere is equal to the mass in the
        % previous step multiplied by the CO2 ratio between this step and
        % the last. Make the simplifying assumption that mass C and CO2
        % directly scale:
        m1(i) = m1(i-1) + CO2in(i);   % change in total C in atmosphere (Gt)
        
        % Carbon transport prop. to amount in atmosphere
        dm(i) = k12(i)/(12*1e15)*sT(i);   % convert moles back to GtC
        
        % update total C in atmosphere and surface ocean according to CO2 exchange
        m1(i) = m1(i) - dm(i);
        m2(i) = m2(i-1) + dm(i);
        
        % CO2 degassing 3.2*10^12 mol/year (Fischer et al. 2019)
        m1(i) = m1(i) + 0.0384*sT(i);   %0.0384 = 3.2*1e12*12/1e15
        
        % CO2 drawdown via weathering; CO2 removed from atmosphere & added as HCO3
            % MC rand 1 & 2 here. Low 1 and high 2 = no relationship to CO2 & high
            % rate, high 1 and low 2 = strongly sensitive to CO2
        w12(i) = ((0.0384)/(1+exp(-(MCrand(1,1)*0.02)*...
            (atCO2out(i-1)-atCO2in(1))))*(atCO2out(i-1)/atCO2in(1))...
            ^(MCrand(1,2)*0.25)+0.0384*1.5)*sT(i);
        % uncomment this to use prescribe no weathering sensitivity to T
        %w12(i) = 0.0384*2*sT(i);
        
        % uncomment below to drive weathering rates independently of
        % temperature - this example isequivalent to gradually moving the 
        % SiW-T relationship upwards by 5% 
%         if age(i)<1e6
%             w12(i) = w12(i);        
%         elseif age(i)<5e6 && age(i)>=1e6 
%             w12(i) = w12(i)*(1+0.05*(age(i)-1e6)/4e6); % 5% in weathering only run
%         else
%             w12(i) = w12(i)*1.05;
%         end
        
        m1(i) = m1(i) - w12(i);
        m2(i) = m2(i) + w12(i);
        if sum(atCO2add)~=0
            m1(i) = atCO2in(i)/280*597;
        end

        % atmospheric CO2 after these fluxes must equal the ratio of the
        % mass of C in the atmosphere box at the end of these calculations
        % to that when the model was initialised. Note that all first model
        % steps were initialised with 597 GtC in the atmosphere
        % irrespective of ultimately desired boundary conditions, hence CO2
        % is normalised to 597 GtC in the model prestep.
        atCO2out(i) = m1(i)/597*280;
        
        % temperatures for co2sys, scaled directly to CO2
        if i<=4000/sT(i)
            deepT = 7*log(mean(m1(1:i,1),'omitnan')/597*280) - 37;
            if deepT<0 % prevent    T from dropping below 0oC
                deepT = 0;
            end
        else
            % temperature is integrated over the last 4 kyr, to account for
            % mixing of warmer/cooler surface water to depth
            deepT = 7*log(mean(m1(i-round(3999/sT(i),0):i,1),'omitnan')/597*280) - 37;
            if deepT<0
                deepT = 0;
            end
        end
        surfT = 5.24*log(m1(i)/597*280) - 6.56;
        % alternatively, hold T constant
        %deepT = 7*log(597/597*280) - 37;
        %surfT = 5.24*log(597/597*280) - 6.56;
        
        % TAlk conc. change must equal mass delivered from weathering,
        % converted from GtC to moles
        TAlkSurf = TAlkSurfMix(i-1) + TAlkSurfMix(i-1)*...
            (w12(i)*1e15/12)/(surfVol*1000*TAlkSurfMix(i-1)/1e6)*TAlkFac;
        
        % DIC in the surface ocean box is equal to the initial DIC
        % multiplied by the proportional change of the total C mass in this
        % box
        DICSurf = ocDICin*m2(i)/m2(1);

        % Slowly decrease the proportion of Ca in the weathering product
%         if age(i)<5e6
%             weatherCa = weatherCa/1.000002;
%         else
%             weatherCa = 0.5;
%         end

        % Delta-Ca must equal Delta-TAlk multipled by weatherCa, the
        % proportion of cations in the weathering product that are Ca2+
        % assume no need to mix Ca between surface and deep boxes due to long
        % residence time
        oceanCa(i) = oceanCa(i-1) + ...
            weatherCa*(TAlkSurfMix(i-1)*...
            (w12(i)*1e15/12)/(surfVol*1000*TAlkSurfMix(i-1)/1e6)*TAlkFac)./1000.*...
            surfVol/oceanVol;
        % What happens if weathering cations are recycled before reaching
        % the ocean?
        % (!) no delivery of cations from enhanced weathering (!)
        % note that 0.7680 is the base weathering rate
%         oceanCa(i) = oceanCa(i-1) + ...
%             weatherCa*((0.7680+m2(i))/m2(1)*ocDICin - m2(i)/m2(1)*ocDICin)./1000.*...
%             surfVol/oceanVol;

        % This is one simple way in which C export &remineralisation in 
        % the deep ocean could be perturbed
%         if i==1e6/sT(i)
%             orgC = orgC - 20;
%         elseif i==2e6/sT(i)
%             orgC = orgC + 20;
%         end

        % Calculate concentration of TAlk and DIC in the deep box from the
        % previous model step (both mix conservatively)
        TAlkDeep = TAlkDeepMix*TAlkRelDeep;
        DICDeep = DICDeepMix*DICRelDeep;    

        % mix these together in proportions defined by water exchange rate (k23)
        TAlkSurfMix(i) = TAlkSurf*(surfVol-k23)/surfVol + TAlkDeep*k23/surfVol;
        DICSurfMix(i) = DICSurf*(surfVol-k23)/surfVol + DICDeep*k23/surfVol;
           
        TAlkDeepMix = TAlkDeep*(deepVol-k23)/deepVol + TAlkSurf*k23/deepVol;
        DICDeepMix = DICDeep*(deepVol-k23)/deepVol + DICSurf*k23/deepVol;

        % CaCO3 flux from surface ocean to sediment
        TAlkAbsSurf = surfVol*1000*TAlkSurfMix(i)/1e6; % abs. TAlk in deep ocean in moles
        DICAbsSurf = surfVol*1000*DICSurfMix(i)/1e6;
        % parameterisation changes CaCO3 relative to the start conditions
        % of this simulation rather than to an absolute value
        surfOmegaFeed = k(co2sysOutSurf(2,3)*oceanCa(2)./10.3.*100,co2sysOutSurf(i-1,3)*oceanCa(i-1)./10.3.*100)...
           ./0.005.*MCrand(1,4);    % feedback 5x more responsive
        % Ca burial varies by +/-20% according to surface ocean Omega
        % note 0.005 is max value of logistic function k i.e. 0.005/0.005*0.5 + 1 = a factor of 1.5
        absChangeSurfOmega(i) = surfOmegaFeed*Caburial*sT(i);
        TAlkSurfMix(i) = ...
            TAlkSurfMix(i)*(1-(2*absChangeSurfOmega(i))/(surfVol*1000*DICSurfMix(i)/1e6)); % CaCO3 burial (/12 to convert Pg into moles)
        DICSurfMix(i) = DICSurfMix(i)*(1-(absChangeSurfOmega(i))/(surfVol*1000*DICSurfMix(i)/1e6));

        % recalculate surface carbonate chemistry after mixing with deep
        co2sysOutSurf(i,:) = ...
            co2syslite_v0_2(surfT,35,0,0,0,TAlkSurfMix(i),DICSurfMix(i),oceanCa(i-1),53,29/oceanCa(i-1));
        m2(i) = m2(i)*DICSurfMix(i)/DICSurf;

        % calculate CaCO3 sediment preservation/dissolution
        OmegaDeep = co2sysDeep(i-1,3)*oceanCa(i-1)/10.3; % Ca corrected Omega
        %c12 = k(OmegaDeep*1000,1*1000);  % transfer function between deep ocean and sediment
        c12(i) = 0.00025/(1+exp(MCrand(1,3)*0.015625*(OmegaDeep*1000-1000)))-0.000125;
        
        % abs. TAlk & DIC in the deep ocean in moles before CaCO3
        % preservation/dissolution
        TAlkAbsDeep = deepVol*1000*TAlkDeepMix/1e6;
        DICAbsDeep = deepVol*1000*DICDeepMix/1e6;

        % The proportional change in TAlk/DIC due to CaCO3
        % preservation/dissolution
        TAlkRelDeep = 1 + (c12(i)*2.5e16)*2/TAlkAbsDeep*sT(i);
        DICRelDeep = 1 + (c12(i)*2.5e16)/DICAbsDeep*sT(i);

        % recalculate deep carbonate chemistry after CaCO3 addition/subtraction
        TAlkDeepOut(i) = TAlkDeepMix*TAlkRelDeep;
        DICDeepOut(i) = DICDeepMix*DICRelDeep+orgC;
        co2sysDeep(i,:) = ...
            co2syslite_v0_2(deepT,35,4000,0,0,TAlkDeepOut(i),DICDeepOut(i),oceanCa(i-1),53,29);

        % calculate change in Ca due to CaCO3 preservation
        % first calculate concentration change in deep ocean DIC in this 
        % step due to CaCO3 burial/dissolution
        DDICdeep = DICDeepMix*DICRelDeep - DICDeepMix;
        % As this is the only process removing DIC from the deep ocean, the
        % Ca concentration change must exactly equal the proportion of DIC
        % removed
        DCa = DDICdeep*deepVol/oceanVol;
        oceanCa(i) = oceanCa(i) + DCa/1000;   % /1000 as DIC in umol, Ca in mmol
        oceanCa(i) = oceanCa(i) - ...
            (surfOmegaFeed*Caburial*sT(i))/DICAbsDeep*DICDeepMix/1000; % change due to CaCO3 burial == abs change in DIC conc/1000 (umol --> mmol)

        % force Ca externally according to addCa variable
        if age(i)>addCa(1) && age(i)<addCa(2)
            oceanCa(i) = oceanCa(i) + addCa(3)*sT(i);
        end
      
        if rem(i,100)==0 && nMC==1
            waitbar(age(i)/ageEnd,w,['current time = ' num2str(round(age(i-1)/1000,0)) 'ka'])
        end
    end

    % collate output data and remove NaNs
    age(2) = 1e-5;   % because top two rows of age both = 0
    % c12*2.5e16 to convert into Pmoles (see code above)
    dataOut{j} = [age m1 m2 TAlkSurfMix DICSurfMix co2sysOutSurf(:,1) ...
        co2sysOutSurf(:,3) TAlkDeepOut DICDeepOut co2sysDeep(:,1) ...
        co2sysDeep(:,3) oceanCa dm c12.*2.5e16 w12 absChangeSurfOmega];

    % remove NaN and non-real rows
    dataOut{j}(sum(isnan(dataOut{j}),2)==size(dataOut{j},2),:) = [];
    Qreal = dataOut{j} == real(dataOut{j});
    dataOut{j} = dataOut{j}(sum(Qreal,2)==size(dataOut{j},2),:);
    % reduce the size of the dataset if the simulation length was very long
    if ageEnd<=1e5
        dataOutS = 10;
    else
        dataOutS = 10^(max(ceil(log10(ageEnd)+1),1))/1000000;
    end
    dataOut{j} = ...
        interp1(dataOut{j}(:,1),dataOut{j}(:,2:size(dataOut{j},2)),...
        1:dataOutS:ageEnd);

    randOut{j} = MCrand;
    if nMC>1
        WaitMessage.Send;
        pause(0.002);
    end
    if j==1 % output step size for plotting
        dataOutStep{j} = dataOutS;
    end
end
toc
if nMC>1
    WaitMessage.Destroy
else
    close(w)  
end
%clearvars -EXCEPT ageEnd dataOut dataOutS dataOutStep figN linlog nMC randOut trans

% If the dataset is very big (>1.5 GB), reduce the size by a factor of 10
dataOutInfo = whos('dataOut');
if ageEnd>1e5 && dataOutInfo.bytes/1e9>1.5
    for i = 1:size(dataOut,1)  
        dataOut{i} = dataOut{i}(1:10:end,:);
    end
    dataOutStep{1} = dataOutStep{1}.*10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The lines below produce a plot of the model results
% First, find axes limits
j = 1;
% initial guess at axis limits
minmaxAx = [min(dataOut{j}(:,11)) max(dataOut{j}(:,11)) ;...
    min(dataOut{j}(:,1))/597*280 max(dataOut{j}(:,1))/597*280 ;...
    min(dataOut{j}(:,3)) max(dataOut{j}(:,3)) ;...
    min(dataOut{j}(:,9)) max(dataOut{j}(:,5)) ;...
   0 max(dataOut{j}(:,6).*dataOut{j}(:,11)./10.3) ;...
    (min(dataOut{j}(:,12)) - dataOut{j}(2,12))/10 (max(dataOut{j}(:,12)) - dataOut{j}(2,12))/10 ;...
    min(dataOut{j}(:,13) - dataOut{j}(:,15)./10) max(dataOut{j}(:,13) - dataOut{j}(:,15)./10);...
    min((dataOut{j}(:,14) - dataOut{j}(2,14))./10) max((dataOut{j}(:,14) - dataOut{j}(2,14))./10)];
for i = 1:nMC
    if isreal(dataOut{i}) && sum(~isnan(dataOut{i}(:,1)))>size(dataOut{i},1)*0.95
        if min(dataOut{i}(:,11))<minmaxAx(1,1)
            minmaxAx(1,1) = min(dataOut{i}(:,11));
        end
        if max(dataOut{i}(:,11))>minmaxAx(1,2)
            minmaxAx(1,2) = max(dataOut{i}(:,11));
        end
        if min(dataOut{i}(:,1))/597*280<minmaxAx(2,1)
            minmaxAx(2,1) = min(dataOut{i}(:,1))/597*280;
        end
        if max(dataOut{i}(:,1))/597*280>minmaxAx(2,2)
            minmaxAx(2,2) = max(dataOut{i}(:,1))/597*280;
        end
        if min(dataOut{i}(:,3))<minmaxAx(3,1)
            minmaxAx(3,1) = min(dataOut{i}(:,3));
        end
        if max(dataOut{i}(:,3))>minmaxAx(3,2)
            minmaxAx(3,2) = max(dataOut{i}(:,3));
        end
        if min(dataOut{i}(:,9))<minmaxAx(4,1)
            minmaxAx(4,1) = min(dataOut{i}(:,9));
        end
        if max(dataOut{i}(:,5))>minmaxAx(4,2)
            minmaxAx(4,2) = max(dataOut{i}(:,5));
        end
        if max(dataOut{i}(:,6).*dataOut{i}(:,11)./10.3)>minmaxAx(5,2)
            minmaxAx(5,2) = max(dataOut{i}(:,6).*dataOut{i}(:,11)./10.3);
        end
        if (min(dataOut{i}(:,12)) - dataOut{i}(2,12))/10<minmaxAx(6,1)
            minmaxAx(6,1) = (min(dataOut{i}(:,12)) - dataOut{i}(2,12))/10;
        end
        if (max(dataOut{i}(:,12)) - dataOut{i}(2,12))/10>minmaxAx(6,2)
            minmaxAx(6,2) = (max(dataOut{i}(:,12)) - dataOut{i}(2,12))/10;
        end
        if min(dataOut{i}(:,13) - dataOut{i}(:,15)./10)<minmaxAx(7,1)
            minmaxAx(7,1) = min(dataOut{i}(:,13) - dataOut{i}(:,15)./10);
        end
        if max(dataOut{i}(:,13) - dataOut{i}(:,15)./10)>minmaxAx(7,2)
            minmaxAx(7,2) = max(dataOut{i}(:,13) - dataOut{i}(:,15)./10);
        end
        if min((dataOut{i}(:,14) - dataOut{i}(2,14))./10)<minmaxAx(8,1)
            minmaxAx(8,1) = min((dataOut{i}(:,14) - dataOut{i}(2,14))./10);
        end
        if max((dataOut{i}(:,14) - dataOut{i}(2,14))./10)>minmaxAx(8,2)
            minmaxAx(8,2) = max((dataOut{i}(:,14) - dataOut{i}(2,14))./10);
        end
    end
end
% define sensible y axis limits
for i = 1:size(minmaxAx,1)
    if minmaxAx(i,1)>minmaxAx(i,2)
        minmaxAx(i,:) = sortrows(minmaxAx(i,:)')';
    end
    if minmaxAx(i,1)/minmaxAx(i,2) < 0 
        minmaxAx(i,:) = minmaxAx(i,:).*[1.05 1.05];
    else
        minmaxAx(i,:) = minmaxAx(i,:).*[0.95 1.05];
    end
    if abs(minmaxAx(i,1))/abs(minmaxAx(i,2)) < 0.1
        minmaxAx(i,1) = minmaxAx(i,1) - minmaxAx(i,2)*0.1;
    elseif abs(minmaxAx(i,1))/abs(minmaxAx(i,2)) > 10
        minmaxAx(i,2) = minmaxAx(i,2) - minmaxAx(i,1)*0.1;
    end
    if minmaxAx(i,1)==minmaxAx(i,2)
        minmaxAx(i,1) = minmaxAx(i,1)-1e-5;
        minmaxAx(i,2) = minmaxAx(i,2)+1e-5;
    end
end
minmaxAx(2,:) = minmaxAx(2,:).*[0.9 1.1];
minmaxAx(4,:) = minmaxAx(4,:)./[0.95 1.05].*[0.99 1.01];
if minmaxAx(7,2)>0
    minmaxAx(7,:) = minmaxAx(7,:).*[1.1 1.1];
else
    minmaxAx(7,:) = minmaxAx(7,:)./[0.95 1.05].*[1.1 0.9];
end

[s,d] = cellfun(@size,dataOut);
maxSize = max([s,d]);
clear s d
% produce sensible x axis ticks
if linlog==1
    linlogT = 'linear';
    if ageEnd<1.05e6
        if ageEnd<5.1e5
            xTickT = 0:0.5e5:1e6;
        else
            xTickT = 0:2e5:1e6;
        end
    else
        if ageEnd>1e7
            xTickT = 0:5e6:10^round(log10(ageEnd),1);
        else
            xTickT = 0:1e6:10^round(log10(ageEnd),0);
        end
    end
else
    linlogT = 'log';
    xTickT = [10 100 1e3 1e4 1e5 1e6 1e7];
end

% plot the figure here
ROtype = whos('randOut');
if ROtype.class ~= "double"
    randOut = cell2mat(randOut);
end

% get screen size information
set(0,'units','centimeter')
PixGet = get(0,'screensize');

close(figure(figN))
H = figure(figN);
set(H,'PaperUnits','centimeter','units','centimeter','papersize',...
    [PixGet(1,3)*1/3 PixGet(1,4)*2/3],'Position',...
    [PixGet(1,3)*0.5 PixGet(1,4)*0.2 PixGet(1,3)*1/3 PixGet(1,4)*2/3],...
    'color',[1 1 1])

H1 = subplot(4,2,1);
hold on
H2 = subplot(4,2,2);
hold on
H3 = subplot(4,2,3);
hold on
H4 = subplot(4,2,4);
hold on
H5 = subplot(4,2,5);
hold on
H6 = subplot(4,2,6);
hold on
H7 = subplot(4,2,7);
hold on
H8 = subplot(4,2,8);
hold on

if nMC<20
    cmap = parula(20);
    cmap = repmat(cmap(1,:),20,1);
else
    cmap = parula(20);
end
for i = 1:nMC
    % only plot simulations that didn't drift into impossible carbonate
    % chemistry space
    if isreal(dataOut{i}) && sum(~isnan(dataOut{i}(:,1)))>size(dataOut{i},1)*0.95
        % plot blue lines if model parameterisations were not perturbed
        % (e.g. if multiple CO2 releases were simulated)
        if size(unique(randOut(:,1)),1)==1
            plotClr = 1;
        else
            plotClr = 20;
        end
        if round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0)>20
            plot(H1,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,11),...
                'color',[cmap(plotClr,:) trans])
            plot(H2,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,1)./597.*280,...
                'color',[cmap(plotClr,:) trans])
            plot(H3,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,3),'color',...
                [cmap(plotClr,:) trans])
            plot(H4,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,5),'color',...
                [cmap(plotClr,:) trans])
            plot(H4,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,9),'color',...
                [cmap(plotClr,:) trans])
            plot(H5,(1:dataOutStep{1}:ageEnd),...
                dataOut{i}(:,6).*dataOut{i}(:,11)./10.3,'color',...
                [cmap(plotClr,:) trans])
            plot(H5,(1:dataOutStep{1}:ageEnd),...
                dataOut{i}(:,10).*dataOut{i}(:,11)./10.3,'color',...
                [cmap(plotClr,:) trans])
            plot(H6,(1:dataOutStep{1}:ageEnd),(dataOut{i}(:,12) - ...
                dataOut{i}(2,12))./10,'color',...
                [cmap(plotClr,:) trans])
            plot(H7,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,13) - ...
                dataOut{i}(:,15)./10,'color',...
                [cmap(plotClr,:) trans])
            plot(H8,(1:dataOutStep{1}:ageEnd),(dataOut{i}(:,14) -...
                dataOut{i}(2,14))./10,'color',...
                [cmap(plotClr,:) trans])
        else
            plot(H1,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,11),...
                'color',[cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H2,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,1)./597.*280,...
                'color',[cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H3,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,3),'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H4,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,5),'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H4,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,9),'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H5,(1:dataOutStep{1}:ageEnd),...
                dataOut{i}(:,6).*dataOut{i}(:,11)./10.3,'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H5,(1:dataOutStep{1}:ageEnd),...
                dataOut{i}(:,10).*dataOut{i}(:,11)./10.3,'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H6,(1:dataOutStep{1}:ageEnd),(dataOut{i}(:,12) - ...
                dataOut{i}(2,12))./10,'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H7,(1:dataOutStep{1}:ageEnd),dataOut{i}(:,13) - ...
                dataOut{i}(:,15)./10,'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
            plot(H8,(1:dataOutStep{1}:ageEnd),(dataOut{i}(:,14) -...
                dataOut{i}(2,14))./10,'color',...
                [cmap(round(randOut(i,1)*60+abs(40/(randOut(i,4)+1)),0),:) trans])
        end
    end
end

ylabel(H1,'[Ca^{2+}] (mM)')
ylabel(H2,'{\it{p}}CO_2')
ylabel(H3,'TAlk (\mumol/kg)')
ylabel(H4,'pH (total)')
ylabel(H5,'\Omega_{calcite}')
ylabel(H6,'\DeltaCO_2 diffusion to ocean (GtC/year)')
ylabel(H7,'CaCO_3 burial (mol/year)')
ylabel(H8,'\Delta SiW (GtC/year)')
xlabel(H7,'time (yr)')
xlabel(H8,'time (yr)')

if nMC>1
    dataOutS = dataOutStep{1};
end

set(H1,'box','on','fontsize',8,'layer','top','position',...
    [0.1 0.78 0.37 0.18],'ylim',minmaxAx(1,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])
set(H2,'box','on','fontsize',8,'layer','top','position',...
    [0.59 0.78 0.37 0.18],'ylim',minmaxAx(2,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])
set(H3,'box','on','fontsize',8,'layer','top','position',...
    [0.1 0.54 0.37 0.18],'ylim',minmaxAx(3,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])
set(H4,'box','on','fontsize',8,'layer','top','position',...
    [0.59 0.54 0.37 0.18],'ylim',minmaxAx(4,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])
set(H5,'box','on','fontsize',8,'layer','top','position',...
    [0.1 0.3 0.37 0.18],'ylim',minmaxAx(5,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])
set(H6,'box','on','fontsize',8,'layer','top','position',...
    [0.59 0.3 0.37 0.18],'ylim',minmaxAx(6,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])
set(H7,'box','on','fontsize',8,'layer','top','position',...
    [0.1 0.06 0.37 0.18],'ylim',minmaxAx(7,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])
set(H8,'box','on','fontsize',8,'layer','top','position',...
    [0.59 0.06 0.37 0.18],'ylim',minmaxAx(8,:),'xscale',linlogT,...
    'xtick',xTickT,'xlim',[dataOutS ageEnd])

% grey backgrounds 
set(H1,'color',[0.9 0.9 0.9])
set(H2,'color',[0.9 0.9 0.9])
set(H3,'color',[0.9 0.9 0.9])
set(H4,'color',[0.9 0.9 0.9])
set(H5,'color',[0.9 0.9 0.9])
set(H6,'color',[0.9 0.9 0.9])
set(H7,'color',[0.9 0.9 0.9])
set(H8,'color',[0.9 0.9 0.9])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the following lines of code to define a new set of stable 
% conditions. First perturb the model and run for >10 Myr, then uncomment
% and run these lines, save the result, and replace loadStable above with
% the new save file name. The clearvars command at the end of the main loop
% will need to be commented out first.
%
%     endLine = size(age,1);
%     co2sysOutSurf(2,:) = co2sysOutSurf(endLine,:);
%     co2sysOutSurf(3:endLine,:) = [];
%     co2sysDeep(2,:) = co2sysDeep(endLine,:);
%     co2sysDeep(3:endLine,:) = [];
%     sT(2) = sT(endLine);
%     sT(3:endLine) = [];
%     dm(2) = dm(endLine);
%     dm(3:endLine) = [];
%     m1(2) = m1(endLine);
%     m1(3:endLine) = [];
%     m2(2) = m2(endLine);
%     m2(3:endLine) = [];
%     k12(2) = k12(endLine);
%     k12(3:endLine) = [];
%     oceanCa(2) = oceanCa(endLine);
%     oceanCa(3:endLine) = [];
%     age(2) = 0;
%     age(3:endLine) = [];
%     atCO2in(2) = atCO2in(endLine);
%     atCO2in(3:endLine) = [];
%     atCO2out(2) = atCO2out(endLine);
%     atCO2out(3:endLine) = [];
%     TAlkSurfMix(2) = TAlkSurfMix(endLine);
%     TAlkSurfMix(3:endLine) = [];
%     DICSurfMix(2) = DICSurfMix(endLine);
%     DICSurfMix(3:endLine) = [];
%     TAlkDeepOut(2) = TAlkDeepOut(endLine);
%     TAlkDeepOut(3:endLine) = [];
%     DICDeepOut(2) = DICDeepOut(endLine);
%     DICDeepOut(3:endLine) = [];
%     absChangeSurfOmega(3:endLine) = [];
%     clearvars -EXCEPT age sT atCO2in atCO2out co2sysDeep co2sysOutSurf ...
%         dm m1 m2 k12 oceanCa absChangeSurfOmega oceanVol surfVol deepVol ...
%         weatherCa DCa ocDICin ocTAlkin orgC k TAlkSurfMix DICSurfMix ...
%         TAlkDeepOut DICDeepOut TAlkDeepMix DICDeepMix TAlkRelDeep ...
%         DICRelDeep 