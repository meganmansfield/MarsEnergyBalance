% Code used to produce graphs in Mansfield et al. (2017) JGR-Planets

%% Given an insolation history, calculate the surface energy balance for an ice surface at T=273.15 K.

for n=[1,2,3,4,5,6,7,8]
    load(strcat('insolation',num2str(n),'d.mat')); %flat surface at equator or 20 latitude
    load(strcat('insolation',num2str(n),'t.mat'));
    load(strcat('insolation',num2str(n),'other.mat'));
    %load(strcat('insolationsloped',num2str(n),'d.mat')); %sloped surface at 20 latitude
    %load(strcat('insolationsloped',num2str(n),'d.mat'));
    %load(strcat('insolationsloped',num2str(n),'other.mat'));
    %load(strcat('insolation40sloped',num2str(n),'dflat.mat')); %surface at 40 latitude
    %load(strcat('insolation40sloped',num2str(n),'deq.mat'));
    %load(strcat('insolation40sloped',num2str(n),'dpole.mat'));
    %iDnet=[iDnetpole,iDnetflat,iDneteq];
    %load(strcat('insolation40sloped',num2str(n),'t.mat'));
    %load(strcat('insolation40sloped',num2str(n),'other.mat'));

    %%%%% Vary climate over time by adding pressure change
    Po=1220; %current pressure in pascals
    to=4.56; %current time since Sun formed in Gyr
    k=1274.224; %units of Pa/s, based on current loss rate from Lillis et al. (2017)
    timesunform=4.56-timetotoday; %for each step in track, time since the sun formed
    
    %this section is to make pretty plot with Maven values
    Pmat=zeros(size(timetotoday,1),6);
    Pmat(:,1)=Po;   %constant pressure
    counter=2;
    for aval=[2.20,2.71,3.2208,3.73,4.24]
       Pmat(:,counter)=Po-k.*to.^aval.*(timesunform.^(1-aval)-to^(1-aval))/(1-aval);
       counter=counter+1;
    end

    delTdiurnal=zeros(size(Pmat));
    roundP=round(Pmat);
    %for each pressure, calculate the diurnal temperature change
    for tstep=1:size(Pmat,1)
        for alpha=1:size(Pmat,2)
            if isnan(Pmat(tstep,alpha))
                delTdiurnal(tstep,alpha)=NaN;
            else
                delTdiurnal(tstep,alpha)=-25.73*log10(Pmat(tstep,alpha).*0.01)+104; %converted pressure from Pa to mbar for fit
            end
        end
    end
    
    delTdiurnal(delTdiurnal<0)=0; %remove temp differences below zero

    %%%%%% Incoming solar SW radiation
    load('rayleigh.mat')
    SWD_in=zeros(size(Pmat,1),size(Pmat,2),size(iDnet,2));
    alpha=0.3; %albedo
    for timestep=1:size(Pmat,1)
        for alphaval=1:size(Pmat,2)
            pdiffd=abs(pout(16,:)-Pmat(timestep,alphaval)); %difference between pressures in matrix and actual pressure at that timestep and alpha value
            [pmind,imind]=min(pdiffd);
            scatterd=scatter(16,imind);
            SWD_in(timestep,alphaval,:)=(1-scatterd).*(1-alpha).*iDnet(timestep,:);% Incoming SW from insolation (diurnal)
        end
    end
    
    %%%%%% Outgoing LW radiation
    sigma=5.670*10^-8; %stefan-boltzmann constant
    T=273.15; %assume surface temperature of 273.15 K
    emis=0.98; %emissivity
    LW_out=emis*sigma*T^4;

    %%%%%% Greenhouse effect
    green_in=zeros(size(Pmat)); %physically motivated using sqrt(P) for line broadening
    for timestep=1:size(Pmat,1)
        for alphaval=1:size(Pmat,2)
            green_in(timestep,alphaval)=9.981+9.662*10^-10*273^4*sqrt(Pmat(timestep,alphaval).*0.01);
        end
    end
    
    %%%%% Evaporative cooling following Kite et al 2013 eqn B4

    % Define Constants
    rh1=0.25; %Mars relative humidity
    rh2=0; %fake relative humidity
    Tsurf=273.15; %Surface temperature in K - held constant at freezing point
    %parameterizing Tbl based on fit to Kahre data
    %Tdiff=10.^(-0.2447.*log10(Pmat)+1.67);
    %Tbl=Tsurf-Tdiff;
    tempchanget=(1./0.75).*(10.^(-0.2507.*log10(Pmat)+2.151)); %fit to Mischna data instead of Kahre is more consistent
    Tbl=Tsurf-lumfraction.*tempchanget; %atmospheric temperature
    mc=0.044; %Molar mass of CO2 in kg/mol
    mw=0.018; %Molar mass of H2O in kg/mol
    R=8.314; %gas constant, units J/K*mol
    g=3.711; %Mars gravity in m/s^2
    Le=2830000; %latent heat of evaporation for water, in J/kg
    adjust=0.73/0.466; %adjustment to match results from Moore and Sears 2006, in W/m^2

    %Calculating intermediate variables
    rho_a=Pmat.*mc./R./Tbl; %density of air
    e_sat=3.69*10^12*exp(-6150/Tsurf); %saturation vapor pressure
    rho_sat=e_sat*mw/R/Tsurf; %saturation vapor density
    delta_nu=(rho_sat./rho_a).*(1-rh1); %difference between atm. and surface water mass fractions
    Da=(1.387.*10.^-5).*(Tbl./273.15).^(3/2).*(10.^5./Pmat); %diffusion coefficient of H2O in CO2
    nu_a=(1.48.*10.^-5).*(R.*Tbl./mc./Pmat).*((240+293.15)./(240+Tbl)).*(Tbl./293.15).^(3/2); %viscosity of air
    delrho=(mc-mw).*e_sat.*(1-rh1)./mc./Pmat; %delta rho/rho

    L_free_latent=Le.*0.14.*delta_nu.*rho_a.*Da.*((nu_a./Da).*(g./nu_a.^2).*delrho).^(1/3);
    L_Moore_adjust=L_free_latent.*adjust; %adjust to match experimental data

    %%%%%% Free Sensible Cooling
    
    %Calculate atmospheric temperature
    %tempchange=(1./0.75).*(10.^(-0.249.*log10(Pmat)+2.046)); %temperature change from fit to logT vs. logP curve from Kahre, at modern solar luminosity
    tempchange=(1./0.75).*(10.^(-0.2507.*log10(Pmat)+2.151)); %testing if fit to Mischna data makes a difference
    Tatm=Tsurf-lumfraction.*tempchange; %atmospheric temperature

    %constants
    cpair=770; %specific heat of atmosphere in J/K/kg, from Kite et al. 2013
    ka=0.0138; %Thermal conductivity of CO2 at approx. -16 Celsius in W/m/K (from http://ws680.nist.gov/publication/get_pdf.cfm?pub_id=907540)
    free_sensible=0.14.*(Tsurf-Tatm).*ka.*((cpair.*nu_a.*rho_a./ka).*(g./nu_a.^2).*delrho).^(1/3);

    %%%%% Forced Sensible Cooling
    
    %calculate atmospheric wind speed
    us=(1./0.925).*(-1.597.*log10(Pmat)+10.283).*(0.811.*lumfraction); % m/s wind from fit to wind vs. logP curve from Kahre, at modern solar luminosity
    Avonk=0.4; %von Karmen constant
    zanem=5.53; %anemometer height in m
    z0=10^-4; %roughness height in m
    A=Avonk^2/log(zanem/z0)^2;
    forced_sensible=rho_a.*cpair.*us.*A.*(Tsurf-Tatm);

    %%%%%% Forced Latent Cooling
    
    Mw=2.99*10^-26; %molecular mass of water in kg
    kb=1.381*10^-23; %Boltzmann constant in J/K
    forced_latent=A.*Le.*Mw./kb./Tbl.*us.*(e_sat).*(1-rh1);

    coolflux=forced_latent+forced_sensible+free_sensible+L_Moore_adjust;
        
    %%%% Conductive Cooling
    
    rho=350; %density in kg/m^3
    thermalcond=0.125; %thermal conductivity in W/m*K
    cp=1751; %specific heat in J/kg*K
    kappa=thermalcond/rho/cp; %thermal diffusivity in m^2/s
    tau=44100; %half of a Martian day in s
    dsd=2.32*sqrt(kappa*tau); %diurnal skin depth in m
    T=2*tau; %Martian day length in s
    omega=2*pi/T; %frequency of oscillations
    cond_cool=1./14400.*2.*(0.5.*delTdiurnal).*thermalcond.*sin(omega.*7200)./dsd./omega; %conductive cooling in W/m^2

    %%%%% Net forcing = incoming - outgoing
    netforceD=zeros(size(SWD_in,1),size(SWD_in,3),size(green_in,2)); %matrix of values for time x latitude x alpha
    for lat=1:size(SWD_in,3) %for each latitude
        for alph=1:size(green_in,2) %for each alpha value
            netforceD(:,lat,alph)=SWD_in(:,alph,lat)-LW_out+green_in(:,alph)-cond_cool(:,alph)-coolflux(:,alph);
        end
    end

    save(strcat('MavenT',num2str(n),'.mat'),'netforceD','timetotoday','timesunform','obl'); %flat surface at equator or 20 latitude
    %save(strcat('MavenslopedT',num2str(n),'.mat'),'netforceD','timetotoday','timesunform','obl'); %sloped surface at 20 latitude
    %save(strcat('Maven40slopedT',num2str(n),'.mat'),'netforceD','timetotoday','timesunform','obl'); %surface at 40 latitude
  
    
%    To create figure 6
%     windowsize=100000;
%     b=(1./windowsize)*ones(1,windowsize);
%     a=1;
%     meanD=zeros(size(netforceD));
%     for row=1:6
%         meanD(:,row)=filter(b,a,netforceD(:,1,row));
%     end
%     figure
%     f3=plot(timesunform,netforceD(:,1,1))
%     hold on
%     plot(timesunform,netforceD(:,1,4))
%     plot(timesunform,netforceD(:,1,5))
%     plot(timesunform(50000:end-50000),meanD(100000:end,1,1),'LineWidth',2)
%     plot(timesunform(50000:end-50000),meanD(100000:end,1,4),'LineWidth',2)
%     plot(timesunform(50000:end-50000),meanD(100000:end,1,5),'LineWidth',2)
%     hold off
%     xlabel('Time [Gyr]')
%     ylabel('Energy Available for Melting [W/m^2]')
%     legend('\alpha=0','\alpha=3.2208','\alpha=3.73')
end

%% Figure 4: Graph of Pressure vs. forcing at modern luminosity

Pmat=[600:1:2*10^6]';
load('PvsT.mat');

lumfraction=0.854; %for 2 Ga

%%%%% Parameterization for Delta T
delTdiurnal=zeros(size(Pmat));
roundP=round(Pmat);
for tstep=1:size(Pmat,1)
    delTdiurnal(tstep,1)=-25.73*log10(roundP(tstep).*0.01)+104; %converted pressure from Pa to mbar for fit
end

delTdiurnal(delTdiurnal<0)=0; %remove temperature differences that are below zero.

%%%%%% Outgoing LW radiation
sigma=5.670*10^-8; %stefan-boltzmann constant
T=273.15; %assume surface temperature of 273.15 K
emis=0.98; %emissivity
LW_out=emis*sigma*T^4;

%%% Greenhouse effect
green_in=zeros(size(Pmat));
for timestep=1:size(Pmat,1)
    for alphaval=1:size(Pmat,2)
        green_in(timestep,alphaval)=9.981+9.662*10^-10*273^4*sqrt(Pmat(timestep,alphaval).*0.01);
    end
end

%%%%% Evaporative cooling following Kite et al 2013 eqn B4

% Define Constants
rh1=0.25; %Mars relative humidity
rh2=0; %fake relative humidity
Tsurf=273.15; %Surface temperature in K - held constant at freezing point
%parameterizing Tbl based on fit to Kahre data
%Tdiff=10.^(-0.2447.*log10(Pmat)+1.67);
%Tbl=Tsurf-Tdiff;
tempchanget=(1./0.75).*(10.^(-0.2507.*log10(Pmat)+2.151)); %fit to Mischna data instead of Kahre is more consistent
Tbl=Tsurf-lumfraction.*tempchanget; %atmospheric temperature
mc=0.044; %Molar mass of CO2 in kg/mol
mw=0.018; %Molar mass of H2O in kg/mol
R=8.314; %gas constant, units J/K*mol
g=3.711; %Mars gravity in m/s^2
Le=2830000; %latent heat of evaporation for water, in J/kg
adjust=0.73/0.466; %adjustment to match results from Moore and Sears 2006, in W/m^2

%Calculating intermediate variables
rho_a=Pmat.*mc./R./Tbl; %density of air
e_sat=3.69.*10.^12.*exp(-6150./Tsurf); %saturation vapor pressure
rho_sat=e_sat*mw/R/Tsurf; %saturation vapor density
delta_nu=(rho_sat./rho_a).*(1-rh1); %difference between atm. and surface water mass fractions
Da=(1.387.*10.^-5).*(Tbl./273.15).^(3/2).*(10.^5./Pmat); %diffusion coefficient of H2O in CO2
nu_a=(1.48.*10.^-5).*(R.*Tbl./mc./Pmat).*((240+293.15)./(240+Tbl)).*(Tbl./293.15).^(3/2); %viscosity of air
delrho=(mc-mw).*e_sat.*(1-rh1)./mc./Pmat; %delta rho/rho

L_free_latent=adjust.*Le.*0.14.*delta_nu.*rho_a.*Da.*((nu_a./Da).*(g./nu_a.^2).*delrho).^(1/3);

%%%%%% Free Sensible Cooling
%Calculate atmospheric temperature
%tempchange=(1./0.75).*(10.^(-0.249.*log10(Pmat)+2.046)); %temperature change from fit to logT vs. logP curve from Kahre, at modern solar luminosity
tempchange=(1./0.75).*(10.^(-0.2507.*log10(Pmat)+2.151)); %testing if fit to Mischna data makes a difference
Tatm=Tsurf-lumfraction.*tempchange; %atmospheric temperature

%constants
cpair=770; %specific heat of atmosphere in J/K/kg, from Kite et al. 2013
ka=0.0138; %Thermal conductivity of CO2 at approx. -16 Celsius in W/m/K (from http://ws680.nist.gov/publication/get_pdf.cfm?pub_id=90754)
free_sensible=0.14.*(Tsurf-Tatm).*ka.*((cpair.*nu_a.*rho_a./ka).*(g./nu_a.^2).*delrho).^(1/3);

%%%%% Forced Sensible Cooling
%calculate atmospheric wind speed
us=(1./0.925).*(-1.597.*log10(Pmat)+10.283).*(0.811.*lumfraction); % m/s wind from fit to wind vs. logP curve from Kahre
Avonk=0.4; %von Karmen constant
zanem=5.53; %anemometer height in m
z0=10^-4; %roughness height in m
A=Avonk^2/log(zanem/z0)^2;
forced_sensible=rho_a.*cpair.*us.*A.*(Tsurf-Tatm);

%%%%%% Forced Latent Cooling
Mw=2.99*10^-26; %molecular mass of water in kg
kb=1.381*10^-23; %Boltzmann constant in J/K
forced_latent=A.*Le.*Mw./kb./Tbl.*us.*(e_sat).*(1-rh1);

%%%%%% Rayleigh Scattering
solarconstant=lumfraction.*1367.7; %in W/m^2
albedo=0.3; %albedo
load('rayleigh.mat')
rayleighd=zeros(size(Pmat));
for timestep=1:size(Pmat,1) %for each timestep
    integrald=0;
    for zenith=1:31
        pdiffd=abs(pout(zenith,:)-Pmat(timestep,1)); %difference between pressures in matrix and actual pressure at that timestep and alpha value
        [pmind,imind]=min(pdiffd);
        integrald=integrald+cos(SZAout(zenith,1)*pi/180)*(1-scatter(zenith,imind))+cos(SZAout(zenith+1,imind)*pi/180)*(1-scatter(zenith+1,imind));
    end
    rayleighd(timestep)=integrald*(1-albedo)*solarconstant*6/pi/30*4/24;
end

%%%% Conductive Cooling
rho=350; %density in kg/m^3
thermalcond=0.125; %thermal conductivity in W/m*K
cp=1751; %specific heat in J/kg*K
kappa=thermalcond/rho/cp; %thermal diffusivity in m^2/s
tau=44100; %half of a Martian day in s
dsd=2.32*sqrt(kappa*tau); %diurnal skin depth in m
T=2*tau; %Martian day length in s
omega=2*pi/T; %frequency of oscillations
cond_cool=1./14400.*2.*(0.5.*delTdiurnal).*thermalcond.*sin(omega.*7200)./dsd./omega; %conductive cooling in W/m^2

netdiurnal=rayleighd+green_in-L_free_latent-forced_latent-free_sensible-forced_sensible-cond_cool-LW_out;

%Figure 4
figure
hold on
plot(Pmat,rayleighd,'LineWidth',2)
plot(Pmat,LW_out*ones(size(Pmat)),'--','LineWidth',2)
plot(Pmat,green_in,'--','LineWidth',2)
plot(Pmat,cond_cool,'LineWidth',2)
plot(Pmat,L_free_latent,'LineWidth',2)
plot(Pmat,forced_latent,'LineWidth',2)
plot(Pmat,free_sensible,'LineWidth',2)
plot(Pmat,forced_sensible,'LineWidth',2)
plot(Pmat,abs(netdiurnal),'LineWidth',3)
title('Pressure vs. Forcing')
ylabel('Heating or Cooling [W/m^2]')
xlabel('Pressure [Pa]')
legend('F_{SW}','F_{LW}','F_{gh}','F_{cond}','F_{lfr}','F_{lfo}','F_{sfr}','F_{sfo}','F_{net}','Orientation','horizontal')%,'Net Mean','Net Diurnal')
hold off

%% Figure 5: Find minimum pressure to cause melting as a function of time
clear all, close all

Pmat=[600:1:2*10^6]';
load('PvsT.mat');

load('bahcall.mat');
smallt=[0:0.000001:8];
smallf=interp1(bahcall(:,1),bahcall(:,2),smallt)';

time=linspace(0.86,4.56,1000);
today=time(end);
timediff=abs(time-today);
timesincesunform=4.56-timediff;
timeindex=round(timesincesunform./0.000001+1);
%find fraction of sun's current solar constant that it was at that time
lumfraction=smallf(timeindex);

minvals=zeros(size(lumfraction));    %minimum pressure to allow melting

currentP=0; %pressure at last luminosity step

for luminosity=1:size(lumfraction,1)
    if currentP>0
        Pmat=[currentP-1500:1:currentP+100]';
    end
    %%%%% New Parameterization for Delta T
    delTdiurnal=zeros(size(Pmat));
    roundP=round(Pmat);
    for tstep=1:size(Pmat,1)
        delTdiurnal(tstep,1)=-25.73*log10(roundP(tstep).*0.01)+104; %converted pressure from Pa to mbar for fit
    end

    delTdiurnal(delTdiurnal<0)=0; %remove temperature differences that are below zero.

    %%%%%% Outgoing LW radiation
    sigma=5.670*10^-8; %stefan-boltzmann constant
    T=273.15; %assume surface temperature of 273.15 K
    emis=0.98; %emissivity
    LW_out=emis*sigma*T^4;

    %%%%%% Greenhouse effect
    green_in=zeros(size(Pmat));
    for timestep=1:size(Pmat,1)
        green_in(timestep,1)=9.981+9.662*10^-10*273^4*sqrt(Pmat(timestep,1).*0.01);
    end
    
    %%%%% Evaporative cooling following Kite et al 2013 eqn B4

    % Define Constants
    rh1=0.25; %Mars relative humidity
    rh2=0; %fake relative humidity
    Tsurf=273.15; %Surface temperature in K - held constant at freezing point
    %parameterizing Tbl based on fit to Kahre data
    %Tdiff=10.^(-0.2447.*log10(Pmat)+1.67);
    %Tbl=Tsurf-Tdiff;
    tempchanget=(1./0.75).*(10.^(-0.2507.*log10(Pmat)+2.151)); %fit to Mischna data instead of Kahre is more consistent
    Tbl=Tsurf-lumfraction(luminosity,1).*tempchanget; %atmospheric temperature
    mc=0.044; %Molar mass of CO2 in kg/mol
    mw=0.018; %Molar mass of H2O in kg/mol
    R=8.314; %gas constant, units J/K*mol
    g=3.711; %Mars gravity in m/s^2
    Le=2830000; %latent heat of evaporation for water, in J/kg
    adjust=0.73/0.466; %adjustment to match results from Moore and Sears 2006, in W/m^2

    %Calculating intermediate variables
    rho_a=Pmat.*mc./R./Tbl; %density of air
    e_sat=3.69.*10.^12.*exp(-6150./Tsurf); %saturation vapor pressure
    rho_sat=e_sat*mw/R/Tsurf; %saturation vapor density
    delta_nu=(rho_sat./rho_a).*(1-rh1); %difference between atm. and surface water mass fractions
    Da=(1.387.*10.^-5).*(Tbl./273.15).^(3/2).*(10.^5./Pmat); %diffusion coefficient of H2O in CO2
    nu_a=(1.48.*10.^-5).*(R.*Tbl./mc./Pmat).*((240+293.15)./(240+Tbl)).*(Tbl./293.15).^(3/2); %viscosity of air
    delrho=(mc-mw).*e_sat.*(1-rh1)./mc./Pmat; %delta rho/rho

    L_free_latent=adjust.*Le.*0.14.*delta_nu.*rho_a.*Da.*((nu_a./Da).*(g./nu_a.^2).*delrho).^(1/3);

    %%%%%% Free Sensible Cooling
    %Calculate atmospheric temperature
    %tempchange=(1./0.75).*(10.^(-0.249.*log10(Pmat)+2.046)); %temperature change from fit to logT vs. logP curve from Kahre, at modern solar luminosity
    tempchange=(1./0.75).*(10.^(-0.2507.*log10(Pmat)+2.151)); %testing if fit to Mischna data makes a difference
    Tatm=Tsurf-lumfraction(luminosity,1).*tempchange; %atmospheric temperature

    %constants
    cpair=770; %specific heat of atmosphere in J/K/kg, from Kite et al. 2013
    ka=0.0138; %Thermal conductivity of CO2 at approx. -16 Celsius in W/m/K (from http://ws680.nist.gov/publication/get_pdf.cfm?pub_id=90754)
    free_sensible=0.14.*(Tsurf-Tatm).*ka.*((cpair.*nu_a.*rho_a./ka).*(g./nu_a.^2).*delrho).^(1/3);

    %%%%% Forced Sensible Cooling
    %calculate atmospheric wind speed
    us=(1./0.925).*(-1.597.*log10(Pmat)+10.283).*(0.811.*lumfraction(luminosity,1)); % m/s wind from fit to wind vs. logP curve from Kahre
    Avonk=0.4; %von Karmen constant
    zanem=5.53; %anemometer height in m
    z0=10^-4; %roughness height in m
    A=Avonk^2/log(zanem/z0)^2;
    forced_sensible=rho_a.*cpair.*us.*A.*(Tsurf-Tatm);

    %%%%%% Forced Latent Cooling
    Mw=2.99*10^-26; %molecular mass of water in kg
    kb=1.381*10^-23; %Boltzmann constant in J/K
    forced_latent=A.*Le.*Mw./kb./Tbl.*us.*(e_sat).*(1-rh1);

    %%%%%% Rayleigh Scattering
    solarconstant=lumfraction(luminosity,1).*1367.7; %in W/m^2
    albedo=0.3; %albedo
    load('rayleigh.mat')
    rayleighd=zeros(size(Pmat));
    for timestep=1:size(Pmat,1) %for each timestep
        integrald=0;
        for zenith=1:31
            pdiffd=abs(pout(zenith,:)-Pmat(timestep,1)); %difference between pressures in matrix and actual pressure at that timestep and alpha value
            [pmind,imind]=min(pdiffd);
            integrald=integrald+cos(SZAout(zenith,1)*pi/180)*(1-scatter(zenith,imind))+cos(SZAout(zenith+1,imind)*pi/180)*(1-scatter(zenith+1,imind));
        end
        rayleighd(timestep)=integrald*(1-albedo)*solarconstant*6/pi/30*4/24;
    end

    %%%% Conductive Cooling
    rho=350; %density in kg/m^3
    thermalcond=0.125; %thermal conductivity in W/m*K
    cp=1751; %specific heat in J/kg*K
    kappa=thermalcond/rho/cp; %thermal diffusivity in m^2/s
    tau=44100; %half of a Martian day in s
    dsd=2.32*sqrt(kappa*tau); %diurnal skin depth in m
    T=2*tau; %Martian day length in s
    omega=2*pi/T; %frequency of oscillations
    cond_cool=1./14400.*2.*(0.5.*delTdiurnal).*thermalcond.*sin(omega.*7200)./dsd./omega; %conductive cooling in W/m^2

    netdiurnal=rayleighd+green_in-L_free_latent-forced_latent-free_sensible-forced_sensible-cond_cool-LW_out;
    
    temp=netdiurnal;
    goodpressures=Pmat(netdiurnal>0);
    minvals(luminosity,1)=min(goodpressures);
    currentP=min(goodpressures);
end

timetopresent=4.56-time(minvals>0);

% Figure 5
figure
plot(timetopresent,minvals(minvals>0))
%% Create Figure 7: Energy available for melting vs. obliquity for the warmest years in the model

warmcut=199;    %how many data points to show
for n=[1,2]
    load(strcat('Maven40slopedT',num2str(n),'.mat'));
    [maxvalsflat,indexflat]=sort(netforceD(:,2,6)); %40 lat, zero slope
    obl40flat=obl(indexflat(end-warmcut:end)); %obliquities for warmest 100 points
    force40flat=netforceD(indexflat(end-warmcut:end),2,6);
    [maxvalspole,indexpole]=sort(netforceD(:,1,6)); %40 lat, 20 poleward slope
    obl40pole=obl(indexpole(end-warmcut:end)); %obliquities for warmest 100 points
    force40pole=netforceD(indexpole(end-warmcut:end),1,6);
    [maxvalseq,indexeq]=sort(netforceD(:,3,6)); %40 lat, 20 equatorward slope
    obl40eq=obl(indexeq(end-warmcut:end)); %obliquities for warmest 100 points
    force40eq=netforceD(indexeq(end-warmcut:end),3,6);
    load(strcat('MavenslopedT',num2str(n),'.mat'))
    [max20s,index20s]=sort(netforceD(:,1,6)); %20 lat, 20 poleward slope
    obl20pole=obl(index20s(end-warmcut:end));
    force20pole=netforceD(index20s(end-warmcut:end),1,6);
    load(strcat('MavenT',num2str(n),'new.mat'))
    [maxeq,indexeq]=sort(netforceD(:,2,6)); %equator, no slope
    obleq=obl(indexeq(end-warmcut:end));
    forceeq=netforceD(indexeq(end-warmcut:end),2,6);
    [max20,index20]=sort(netforceD(:,1,6)); %20 lat, no slope
    obl20=obl(index20(end-warmcut:end));
    force20=netforceD(index20(end-warmcut:end),1,6);
    
    x=[force40pole;force40flat;forceeq;force20];
    y=[obl40pole;obl40flat;obleq;obl20];
    lat1=60*ones(200,1);
    lat2=40*ones(200,1);
    lat3=20*ones(200,1);
    lat4=0*ones(200,1);
    lats=[lat1;lat2;lat3;lat4];
    
    figure
    hold on
    scatterhist(x,y,'Group',lats,'MarkerSize',6,'Color','krbm','Marker','^oxs','Direction','out','Style','bar')
    hold off
    xlabel('Energy Available for Melting [W/m^2]')
    ylabel('Obliquity [degrees]')
end

%% Make graphs that show net energy vs. years exceeded for equator
for n=[1]%,2]
    %load(strcat('MavenT',num2str(n),'.mat'))   %flat surface at equator or 20 latitude
    load(strcat('MavenslopedT',num2str(n),'.mat'))  %sloped surface at 20 latitude
    %load(strcat('Maven40slopedT',num2str(n),'.mat'))   %surface at 40 latitude
    for alph=1:size(netforceD,3)
        netforceD(obl<40,1,alph)=NaN;   %for 0 or 20 latitude
        %netforceD(obl>40,1,alph)=NaN;  %for 40 latitude
        %netforceD(obl<30,1,alph)=NaN;
    end

    latnum=1;
    logforceD=log10(netforceD(:,latnum,:));
    
    nbins=100000;
    binsD=zeros(size(netforceD,3),nbins);
    xD=zeros(size(binsD));
    cumnumD=zeros(size(binsD));
    for alpha=1:size(netforceD,3)
        binsD(alpha,:)=hist(netforceD(:,latnum,alpha),100000);
        xD(alpha,:)=linspace(min(netforceD(:,latnum,alpha)),max(netforceD(:,latnum,alpha)),size(binsD,2));
        for step=1:size(binsD,2)
                cumnumD(alpha,step)=sum(binsD(alpha,step:size(binsD,2)))*1998.6; %1998.6=mean time difference between timesteps
        end
    end
    miny=min(min(log10(cumnumD)))-0.1;
    maxy=max(max(log10(cumnumD)))+0.1;
    
    %Figure 8
    figure
    hold on
    for alpha=2:size(netforceD,3)
        plot(xD(alpha,:),log10(cumnumD(alpha,:)),'LineWidth',2)
    end
    plot([15,15],[3 10],'r--','LineWidth',2)
    hold off
    title(['Energy Balance Histogram, 20 Latitude, 20 Degree Poleward Slope, Track ',num2str(n)])
    ylabel('log_{10}(Years in which Energy Exceeded)')
    xlabel('Diurnal Net Surface Energy [W/m^2]')
    box on
    grid minor
    legend('\alpha=2.20 (MAVEN -1\sigma)','\alpha=2.71 (MAVEN -0.5\sigma)','\alpha=3.2208 (MAVEN mean)','\alpha=3.73 (MAVEN +0.5\sigma)','\alpha=4.24 (MAVEN +1\sigma)')
    
end

%% Figure 9: Histogram of dry spells, marginalizing over all latitude and slope options
lengthlist=[];
drylengthlist=[];
longestdries=[];
alpha=5;
for n=[1]
% 0 Latitude, No Slope
    load(strcat('MavenT',num2str(n),'.mat'))
    netforceD(obl<40,2,alpha+1)=NaN;
    goodind=netforceD(:,2,alpha+1)>(15); %start from point >15
    goodind2=netforceD(:,2,alpha+1)>(0); %end when you drop below 0
    index=1:size(netforceD,1);
    start=index(goodind); %years in which E>15
    stop=index(goodind2); %years in which E>0
    timediffs=diff(timesunform); %actual timesteps from simulation
    d1=[1,diff(start)];
    startind=start(d1~=1); %years where E>15 that timestep and E<15 in the last timestep
    d2=diff(stop);
    stopind=stop(d2~=1); %years where E>0 that timestep and E<0 in the next timestep
    if stop(end)==size(timesunform,1) %if the most recent year has melting
        stopind=[stopind,size(timesunform,1)]; %add that last year to the list of stopping indices
    end
    length=zeros(size(startind,2),1);
    drylength=zeros(size(startind,2),1); %number of dry years between wet spells
    temp=0;
    lcount=1;
    for i=1:size(startind,2)
        if startind(i)>temp %if the starting index is later than the last stopping year
            goodstop=stopind(stopind>startind(i)); %possible stopping points
            if size(goodstop,2)>0 %if there are any stopping points after that starting point
                if goodstop(1)==size(timesunform,1)
                    length(lcount)=sum(timediffs(startind(i):end));
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                else
                    length(lcount)=sum(timediffs(startind(i):goodstop(1))); %number of timesteps in each melting period (stopping time - starting time + 1)
                    if i~=size(startind,2)
                        drylength(lcount)=sum(timediffs(goodstop(1):startind(i+1)));
                    end
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                end
                temp=goodstop(1); %the starting index is now the end of the last melting period
                lcount=lcount+1;
            end
        end
    end
    lengtheq=length(length~=0); %remove the empty placeholders in length
    drylength=drylength(drylength~=0);

    lengthlist=[lengthlist;lengtheq];
    drylengthlist=[drylengthlist;drylength];
    longestdries=[longestdries;max(drylength)];

% 20 Latitude, No Slope
    load(strcat('MavenT',num2str(n),'.mat'))
    netforceD(obl<40,1,alpha+1)=NaN;
    goodind=netforceD(:,1,alpha+1)>(15); %start from point >15
    goodind2=netforceD(:,1,alpha+1)>(0); %end when you drop below 0
    index=1:size(netforceD,1);
    start=index(goodind); %years in which E>15
    stop=index(goodind2); %years in which E>0
    timediffs=diff(timesunform); %actual timesteps from simulation
    d1=[1,diff(start)];
    startind=start(d1~=1); %years where E>15 that timestep and E<15 in the last timestep
    d2=diff(stop);
    stopind=stop(d2~=1); %years where E>0 that timestep and E<0 in the next timestep
    if stop(end)==size(timesunform,1) %if the most recent year has melting
        stopind=[stopind,size(timesunform,1)]; %add that last year to the list of stopping indices
    end
    length=zeros(size(startind,2),1);
    drylength=zeros(size(startind,2),1); %number of dry years between wet spells
    temp=0;
    lcount=1;
    for i=1:size(startind,2)
        if startind(i)>temp %if the starting index is later than the last stopping year
            goodstop=stopind(stopind>startind(i)); %possible stopping points
            if size(goodstop,2)>0 %if there are any stopping points after that starting point
                if goodstop(1)==size(timesunform,1)
                    length(lcount)=sum(timediffs(startind(i):end));
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                else
                    length(lcount)=sum(timediffs(startind(i):goodstop(1))); %number of timesteps in each melting period (stopping time - starting time + 1)
                    if i~=size(startind,2)
                        drylength(lcount)=sum(timediffs(goodstop(1):startind(i+1)));
                    end
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                end
                temp=goodstop(1); %the starting index is now the end of the last melting period
                lcount=lcount+1;
            end
        end
    end
    length20=length(length~=0); %remove the empty placeholders in length
    drylength=drylength(drylength~=0);

    lengthlist=[lengthlist;length20];
    drylengthlist=[drylengthlist;drylength];
    longestdries=[longestdries;max(drylength)];

% 20 Latitude, 20 Degree Poleward Slope
    load(strcat('MavenslopedT',num2str(n),'.mat'))
    netforceD(obl<40,1,alpha+1)=NaN;
    goodind=netforceD(:,1,alpha+1)>(15); %start from point >15
    goodind2=netforceD(:,1,alpha+1)>(0); %end when you drop below 0
    index=1:size(netforceD,1);
    start=index(goodind); %years in which E>15
    stop=index(goodind2); %years in which E>0
    timediffs=diff(timesunform); %actual timesteps from simulation
    d1=[1,diff(start)];
    startind=start(d1~=1); %years where E>15 that timestep and E<15 in the last timestep
    d2=diff(stop);
    stopind=stop(d2~=1); %years where E>0 that timestep and E<0 in the next timestep
    if stop(end)==size(timesunform,1) %if the most recent year has melting
        stopind=[stopind,size(timesunform,1)]; %add that last year to the list of stopping indices
    end
    length=zeros(size(startind,2),1);
    drylength=zeros(size(startind,2),1); %number of dry years between wet spells
    temp=0;
    lcount=1;
    for i=1:size(startind,2)
        if startind(i)>temp %if the starting index is later than the last stopping year
            goodstop=stopind(stopind>startind(i)); %possible stopping points
            if size(goodstop,2)>0 %if there are any stopping points after that starting point
                if goodstop(1)==size(timesunform,1)
                    length(lcount)=sum(timediffs(startind(i):end));
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                else
                    length(lcount)=sum(timediffs(startind(i):goodstop(1))); %number of timesteps in each melting period (stopping time - starting time + 1)
                    if i~=size(startind,2)
                        drylength(lcount)=sum(timediffs(goodstop(1):startind(i+1)));
                    end
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                end
                temp=goodstop(1); %the starting index is now the end of the last melting period
                lcount=lcount+1;
            end
        end
    end
    length20pole=length(length~=0); %remove the empty placeholders in length
    drylength=drylength(drylength~=0);

    lengthlist=[lengthlist;length20pole];
    drylengthlist=[drylengthlist;drylength];
    longestdries=[longestdries;max(drylength)];

% 40 Latitude, No Slope
    load(strcat('Maven40slopedT',num2str(n),'.mat'))
    netforceD(obl<30,2,alpha+1)=NaN;
    netforceD(obl>40,2,alpha+1)=NaN;
    goodind=netforceD(:,2,alpha+1)>(15); %start from point >15
    goodind2=netforceD(:,2,alpha+1)>(0); %end when you drop below 0
    index=1:size(netforceD,1);
    start=index(goodind); %years in which E>15
    stop=index(goodind2); %years in which E>0
    timediffs=diff(timesunform); %actual timesteps from simulation
    d1=[1,diff(start)];
    startind=start(d1~=1); %years where E>15 that timestep and E<15 in the last timestep
    d2=diff(stop);
    stopind=stop(d2~=1); %years where E>0 that timestep and E<0 in the next timestep
    if stop(end)==size(timesunform,1) %if the most recent year has melting
        stopind=[stopind,size(timesunform,1)]; %add that last year to the list of stopping indices
    end
    length=zeros(size(startind,2),1);
    drylength=zeros(size(startind,2),1); %number of dry years between wet spells
    temp=0;
    lcount=1;
    for i=1:size(startind,2)
        if startind(i)>temp %if the starting index is later than the last stopping year
            goodstop=stopind(stopind>startind(i)); %possible stopping points
            if size(goodstop,2)>0 %if there are any stopping points after that starting point
                if goodstop(1)==size(timesunform,1)
                    length(lcount)=sum(timediffs(startind(i):end));
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                else
                    length(lcount)=sum(timediffs(startind(i):goodstop(1))); %number of timesteps in each melting period (stopping time - starting time + 1)
                    if i~=size(startind,2)
                        drylength(lcount)=sum(timediffs(goodstop(1):startind(i+1)));
                    end
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                end
                temp=goodstop(1); %the starting index is now the end of the last melting period
                lcount=lcount+1;
            end
        end
    end
    length40=length(length~=0); %remove the empty placeholders in length
    drylength=drylength(drylength~=0);

    lengthlist=[lengthlist;length40];
    drylengthlist=[drylengthlist;drylength];
    longestdries=[longestdries;max(drylength)];

% 40 Latitude, 20 Degree Poleward Slope
    load(strcat('Maven40slopedT',num2str(n),'.mat'))
    netforceD(obl<30,1,alpha+1)=NaN;
    netforceD(obl>40,1,alpha+1)=NaN;
    goodind=netforceD(:,1,alpha+1)>(15); %start from point >15
    goodind2=netforceD(:,1,alpha+1)>(0); %end when you drop below 0
    index=1:size(netforceD,1);
    start=index(goodind); %years in which E>15
    stop=index(goodind2); %years in which E>0
    timediffs=diff(timesunform); %actual timesteps from simulation
    d1=[1,diff(start)];
    startind=start(d1~=1); %years where E>15 that timestep and E<15 in the last timestep
    d2=diff(stop);
    stopind=stop(d2~=1); %years where E>0 that timestep and E<0 in the next timestep
    if stop(end)==size(timesunform,1) %if the most recent year has melting
        stopind=[stopind,size(timesunform,1)]; %add that last year to the list of stopping indices
    end
    length=zeros(size(startind,2),1);
    drylength=zeros(size(startind,2),1); %number of dry years between wet spells
    temp=0;
    lcount=1;
    for i=1:size(startind,2)
        if startind(i)>temp %if the starting index is later than the last stopping year
            goodstop=stopind(stopind>startind(i)); %possible stopping points
            if size(goodstop,2)>0 %if there are any stopping points after that starting point
                if goodstop(1)==size(timesunform,1)
                    length(lcount)=sum(timediffs(startind(i):end));
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                else
                    length(lcount)=sum(timediffs(startind(i):goodstop(1))); %number of timesteps in each melting period (stopping time - starting time + 1)
                    if i~=size(startind,2)
                        drylength(lcount)=sum(timediffs(goodstop(1):startind(i+1)));
                    end
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                end
                temp=goodstop(1); %the starting index is now the end of the last melting period
                lcount=lcount+1;
            end
        end
    end
    length40pole=length(length~=0); %remove the empty placeholders in length
    drylength=drylength(drylength~=0);

    lengthlist=[lengthlist;length40pole];
    drylengthlist=[drylengthlist;drylength];
    longestdries=[longestdries;max(drylength)];

% 40 Latitude, 20 Degree Equatorward Slope
    load(strcat('Maven40slopedT',num2str(n),'.mat'))
    netforceD(obl<30,3,alpha+1)=NaN;
    netforceD(obl>40,3,alpha+1)=NaN;
    goodind=netforceD(:,3,alpha+1)>(15); %start from point >15
    goodind2=netforceD(:,3,alpha+1)>(0); %end when you drop below 0
    index=1:size(netforceD,1);
    start=index(goodind); %years in which E>15
    stop=index(goodind2); %years in which E>0
    timediffs=diff(timesunform); %actual timesteps from simulation
    d1=[1,diff(start)];
    startind=start(d1~=1); %years where E>15 that timestep and E<15 in the last timestep
    d2=diff(stop);
    stopind=stop(d2~=1); %years where E>0 that timestep and E<0 in the next timestep
    if stop(end)==size(timesunform,1) %if the most recent year has melting
        stopind=[stopind,size(timesunform,1)]; %add that last year to the list of stopping indices
    end
    length=zeros(size(startind,2),1);
    drylength=zeros(size(startind,2),1); %number of dry years between wet spells
    temp=0;
    lcount=1;
    for i=1:size(startind,2)
        if startind(i)>temp %if the starting index is later than the last stopping year
            goodstop=stopind(stopind>startind(i)); %possible stopping points
            if size(goodstop,2)>0 %if there are any stopping points after that starting point
                if goodstop(1)==size(timesunform,1)
                    length(lcount)=sum(timediffs(startind(i):end));
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                else
                    length(lcount)=sum(timediffs(startind(i):goodstop(1))); %number of timesteps in each melting period (stopping time - starting time + 1)
                    if i~=size(startind,2)
                        drylength(lcount)=sum(timediffs(goodstop(1):startind(i+1)));
                    end
                    counter=startind(i);
                    for in=startind(i):goodstop(1)
                        if any(start==in)
                            counter=in;
                        else
                        end
                    end
                end
                temp=goodstop(1); %the starting index is now the end of the last melting period
                lcount=lcount+1;
            end
        end
    end
    length40eq=length(length~=0); %remove the empty placeholders in length
    drylength=drylength(drylength~=0);

    lengthlist=[lengthlist;length40eq];
    drylengthlist=[drylengthlist;drylength];
    longestdries=[longestdries;max(drylength)];
    
    temp=lengthlist.*10^9;
    bincents=[3.7,4.0,4.3,4.6,4.9,5.2,5.5,5.8,6.1,6.4,6.7,7.0,7.3,7.6]; %centers of bins
    
    %Figure 9
    figure
    [n,xout]=hist(log10(temp),bincents);
    bar(xout, n, 'barwidth', 1, 'basevalue', 1);
    set(gca,'YScale','log')
    xlabel('log(Length of Wet Events [years])')
    ylabel('Number of Events')
end

%% Calculate geologic constraint values

longlake=-ones(3,8);    %Constraint 1: longlake>10^4 years (10 kyr)
tolivine=-ones(3,8);    %Constraint 2: tolivine<10^7 years (10 Myr)
cratercount=-ones(3,8); %Constraint 3: cratercount>10^7 years (10 Myr)
Ebigger=-ones(3,8);
Ebig=-ones(3,8);
runnum=1;
shift=0; %amount to move energy balance due to non-CO2 forcing.
for n=[1]%,2,3,4,5,6,7,8]
    for alpha=[0:6]
        %load(strcat('MavenT',num2str(n),'.mat'))   %flat surface at equator or 20 latitude
        load(strcat('MavenslopedT',num2str(n),'new.mat'))  %sloped surface at 20 latitude
        %load(strcat('Maven40slopedT',num2str(n),'.mat'))   %surface at 40 latitude
        netforceD(obl<40,3,alpha+1)=NaN;
        goodind=netforceD(:,1,alpha+1)>(15-shift); %start from point >15
        goodind2=netforceD(:,1,alpha+1)>(0-shift); %end when you drop below 0
        index=1:size(netforceD,1);
        start=index(goodind); %years in which E>15
        stop=index(goodind2); %years in which E>0
        timediffs=diff(timesunform); %actual timesteps from simulation
        if isempty(start)
            tolivine(alpha+1,runnum)=0;
            longlake(alpha+1,runnum)=0;
            cratercount(alpha+1,runnum)=0;
            Ebig(alpha+1,runnum)=sum(timediffs(stop))*10^9; %total amount of time with E>0 in years
            Ebigger(alpha+1,runnum)=0;
        else
            d1=[0,diff(start)];
            startind=start(d1~=1); %years where E>15 that timestep and E<15 in the last timestep
            d2=diff(stop);
            stopind=stop(d2~=1); %years where E>0 that timestep and E<0 in the next timestep
            if stop(end)==size(timesunform,1) %if the most recent year has melting
                stopind=[stopind,size(timesunform,1)]; %add that last year to the list of stopping indices
            end
            length=zeros(size(startind,2),1);
            olilength=zeros(size(startind,2),1); %number of years with 0<E<15 after the energy balance drops below E=15
            drylength=zeros(size(startind,2),1); %number of dry years between wet spells
            temp=0;
            lcount=1;
            test=0;
            for i=1:size(startind,2)
                if startind(i)>temp %if the starting index is later than the last stopping year
                    goodstop=stopind(stopind>startind(i)); %possible stopping points
                    if size(goodstop,2)>0 %if there are any stopping points after that starting point
                        if goodstop(1)==size(timesunform,1)
                            length(lcount)=sum(timediffs(startind(i):end));
                            counter=startind(i);
                            for index=startind(i):goodstop(1)
                                if any(start==index)
                                    counter=index;
                                else
                                end
                            end
                            olilength(lcount)=sum(timediffs(counter:end));
                        else
                            length(lcount)=sum(timediffs(startind(i):goodstop(1))); %number of timesteps in each melting period (stopping time - starting time + 1)
                            if i~=size(startind,2)
                                drylength(lcount)=sum(timediffs(goodstop(1):startind(i+1)));
                            end
                            counter=startind(i);
                            for index=startind(i):goodstop(1)
                                if any(start==index)
                                    counter=index;
                                else
                                end
                            end
                            olilength(lcount)=sum(timediffs(counter:goodstop(1)));
                        end
                        temp=goodstop(1); %the starting index is now the end of the last melting period
                        lcount=lcount+1;
                    end
                end
            end
            length=length(length~=0); %remove the empty placeholders in length
            olilength=olilength(olilength~=0);
            drylength=drylength(drylength~=0);
            yearlength=length.*10^9; %gives lengths in years

            tolivine(alpha+1,runnum)=olilength(end)*10^9/10^6; %length of last lake-forming episode must be less than 10^7 years in order to observe olivine (in Myr)
            longlake(alpha+1,runnum)=max(length)*10^9; %length of longest lake-forming episode
            if stopind(end)==size(timesunform,1)
                cratercount(alpha+1,runnum)=sum(timediffs(startind(1):end))*10^9/10^6; %in Myr
            else
                cratercount(alpha+1,runnum)=sum(timediffs(startind(1):stopind(end)))*10^9/10^6; %total amount of time with liquid water on the surface in Myr
            end
            Ebigger(alpha+1,runnum)=sum(timediffs(start))*10^9/10^3; %total amount of time with E>15 in kyr
            Ebig(alpha+1,runnum)=sum(timediffs(stop))*10^9; %total amount of time with E>0 in years
        end
    end
    runnum=runnum+1;
end

%% Fit to calculate greenhouse forcing

%read in data
clear all
load('wrf_tsk_range_ob15_p0.mat')
maxT15p0=maxT;
green15p0=maxgreenmaxT;
load('wrf_tsk_range_ob15_p6.mat')
maxT15p6=maxT;
green15p6=maxgreenmaxT;
Trange15p6=Trange;
load('wrf_tsk_range_ob15_p60.mat')
maxT15p60=maxT;
green15p60=maxgreenmaxT;
Trange15p60=Trange;
load('wrf_tsk_range_ob15_p600.mat')
maxT15p600=maxT;
green15p600=maxgreenmaxT;
Trange15p600=Trange;
load('wrf_tsk_range_ob15_p1200.mat')
maxT15p1200=maxT;
green15p1200=maxgreenmaxT;
Trange15p1200=Trange;
load('wrf_tsk_range_ob25_p0.mat')
maxT25p0=maxT;
green25p0=maxgreenmaxT;
load('wrf_tsk_range_ob25_p6.mat')
maxT25p6=maxT;
green25p6=maxgreenmaxT;
Trange25p6=Trange;
load('wrf_tsk_range_ob25_p60.mat')
maxT25p60=maxT;
green25p60=maxgreenmaxT;
Trange25p60=Trange;
load('wrf_tsk_range_ob25_p600.mat')
maxT25p600=maxT;
green25p600=maxgreenmaxT;
Trange25p600=Trange;
load('wrf_tsk_range_ob25_p1200.mat')
maxT25p1200=maxT;
green25p1200=maxgreenmaxT;
Trange25p1200=Trange;
load('wrf_tsk_range_ob35_p0.mat')
maxT35p0=maxT;
green35p0=maxgreenmaxT;
load('wrf_tsk_range_ob35_p6.mat')
maxT35p6=maxT;
green35p6=maxgreenmaxT;
Trange35p6=Trange;
load('wrf_tsk_range_ob35_p60.mat')
maxT35p60=maxT;
green35p60=maxgreenmaxT;
Trange35p60=Trange;
load('wrf_tsk_range_ob35_p600.mat')
maxT35p600=maxT;
green35p600=maxgreenmaxT;
Trange35p600=Trange;
load('wrf_tsk_range_ob35_p1200.mat')
maxT35p1200=maxT;
green35p1200=maxgreenmaxT;
Trange35p1200=Trange;
load('wrf_tsk_range_ob45_p0.mat')
maxT45p0=maxT;
green45p0=maxgreenmaxT;
load('wrf_tsk_range_ob45_p6.mat')
maxT45p6=maxT;
green45p6=maxgreenmaxT;
Trange45p6=Trange;
load('wrf_tsk_range_ob45_p60.mat')
maxT45p60=maxT;
green45p60=maxgreenmaxT;
Trange45p60=Trange;
load('wrf_tsk_range_ob45_p600.mat')
maxT45p600=maxT;
green45p600=maxgreenmaxT;
Trange45p600=Trange;
load('wrf_tsk_range_ob45_p1200.mat')
maxT45p1200=maxT;
green45p1200=maxgreenmaxT;
Trange45p1200=Trange;
load('wrf_tsk_range_ob60_p0.mat')
maxT60p0=maxT;
green60p0=maxgreenmaxT;
load('wrf_tsk_range_ob60_p6.mat')
maxT60p6=maxT;
green60p6=maxgreenmaxT;
Trange60p6=Trange;
load('wrf_tsk_range_ob60_p60.mat')
maxT60p60=maxT;
green60p60=maxgreenmaxT;
Trange60p60=Trange;
load('wrf_tsk_range_ob60_p600.mat')
maxT60p600=maxT;
green60p600=maxgreenmaxT;
Trange60p600=Trange;
load('wrf_tsk_range_ob60_p1200.mat')
maxT60p1200=maxT;
green60p1200=maxgreenmaxT;
Trange60p1200=Trange;

% make cuts to ignore poles, topo<-6000, topo>12000, albedo<0.15, ti<150
load('albedo.mat')
load('topo.mat')
load('ti.mat')
lat=[-87.5 -82.5 -77.5 -72.5 -67.5 -62.5 -57.5 -52.5 -47.5 -42.5 -37.5 -32.5 -27.5 -22.5 -17.5 -12.5 -7.5 -2.5  2.5  7.5 12.5 17.5 22.5 27.5 32.5 37.5 42.5 47.5 52.5 57.5 62.5 67.5 72.5 77.5 82.5 87.5];
long=5*[1:72];

albedo=reshape(albedo(14:23,:),[],1);
topo=reshape(topo(14:23,:),[],1);
ti=reshape(ti(14:23,:),[],1);

maxT15p0=reshape(maxT15p0(14:23,:),[],1);
cmaxT15p0=maxT15p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT15p6=reshape(maxT15p6(14:23,:),[],1);
cmaxT15p6=maxT15p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT15p60=reshape(maxT15p60(14:23,:),[],1);
cmaxT15p60=maxT15p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT15p600=reshape(maxT15p600(14:23,:),[],1);
cmaxT15p600=maxT15p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT15p1200=reshape(maxT15p1200(14:23,:),[],1);
cmaxT15p1200=maxT15p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT25p0=reshape(maxT25p0(14:23,:),[],1);
cmaxT25p0=maxT25p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT25p6=reshape(maxT25p6(14:23,:),[],1);
cmaxT25p6=maxT25p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT25p60=reshape(maxT25p60(14:23,:),[],1);
cmaxT25p60=maxT25p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT25p600=reshape(maxT25p600(14:23,:),[],1);
cmaxT25p600=maxT25p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT25p1200=reshape(maxT25p1200(14:23,:),[],1);
cmaxT25p1200=maxT25p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT35p0=reshape(maxT35p0(14:23,:),[],1);
cmaxT35p0=maxT35p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT35p6=reshape(maxT35p6(14:23,:),[],1);
cmaxT35p6=maxT35p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT35p60=reshape(maxT35p60(14:23,:),[],1);
cmaxT35p60=maxT35p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT35p600=reshape(maxT35p600(14:23,:),[],1);
cmaxT35p600=maxT35p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT35p1200=reshape(maxT35p1200(14:23,:),[],1);
cmaxT35p1200=maxT35p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT45p0=reshape(maxT45p0(14:23,:),[],1);
cmaxT45p0=maxT45p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT45p6=reshape(maxT45p6(14:23,:),[],1);
cmaxT45p6=maxT45p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT45p60=reshape(maxT45p60(14:23,:),[],1);
cmaxT45p60=maxT45p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT45p600=reshape(maxT45p600(14:23,:),[],1);
cmaxT45p600=maxT45p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT45p1200=reshape(maxT45p1200(14:23,:),[],1);
cmaxT45p1200=maxT45p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT60p0=reshape(maxT60p0(14:23,:),[],1);
cmaxT60p0=maxT60p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT60p6=reshape(maxT60p6(14:23,:),[],1);
cmaxT60p6=maxT60p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT60p60=reshape(maxT60p60(14:23,:),[],1);
cmaxT60p60=maxT60p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT60p600=reshape(maxT60p600(14:23,:),[],1);
cmaxT60p600=maxT60p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
maxT60p1200=reshape(maxT60p1200(14:23,:),[],1);
cmaxT60p1200=maxT60p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);

green15p0=reshape(green15p0(14:23,:),[],1);
cgreen15p0=green15p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green15p6=reshape(green15p6(14:23,:),[],1);
cgreen15p6=green15p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green15p60=reshape(green15p60(14:23,:),[],1);
cgreen15p60=green15p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green15p600=reshape(green15p600(14:23,:),[],1);
cgreen15p600=green15p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green15p1200=reshape(green15p1200(14:23,:),[],1);
cgreen15p1200=green15p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green25p0=reshape(green25p0(14:23,:),[],1);
cgreen25p0=green25p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green25p6=reshape(green25p6(14:23,:),[],1);
cgreen25p6=green25p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green25p60=reshape(green25p60(14:23,:),[],1);
cgreen25p60=green25p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green25p600=reshape(green25p600(14:23,:),[],1);
cgreen25p600=green25p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green25p1200=reshape(green25p1200(14:23,:),[],1);
cgreen25p1200=green25p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green35p0=reshape(green35p0(14:23,:),[],1);
cgreen35p0=green35p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green35p6=reshape(green35p6(14:23,:),[],1);
cgreen35p6=green35p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green35p60=reshape(green35p60(14:23,:),[],1);
cgreen35p60=green35p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green35p600=reshape(green35p600(14:23,:),[],1);
cgreen35p600=green35p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green35p1200=reshape(green35p1200(14:23,:),[],1);
cgreen35p1200=green35p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green45p0=reshape(green45p0(14:23,:),[],1);
cgreen45p0=green45p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green45p6=reshape(green45p6(14:23,:),[],1);
cgreen45p6=green45p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green45p60=reshape(green45p60(14:23,:),[],1);
cgreen45p60=green45p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green45p600=reshape(green45p600(14:23,:),[],1);
cgreen45p600=green45p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green45p1200=reshape(green45p1200(14:23,:),[],1);
cgreen45p1200=green45p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green60p0=reshape(green60p0(14:23,:),[],1);
cgreen60p0=green60p0(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green60p6=reshape(green60p6(14:23,:),[],1);
cgreen60p6=green60p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green60p60=reshape(green60p60(14:23,:),[],1);
cgreen60p60=green60p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green60p600=reshape(green60p600(14:23,:),[],1);
cgreen60p600=green60p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
green60p1200=reshape(green60p1200(14:23,:),[],1);
cgreen60p1200=green60p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);

Trange15p6=reshape(Trange15p6(14:23,:),[],1);
cTrange15p6=Trange15p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange15p60=reshape(Trange15p60(14:23,:),[],1);
cTrange15p60=Trange15p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange15p600=reshape(Trange15p600(14:23,:),[],1);
cTrange15p600=Trange15p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange15p1200=reshape(Trange15p1200(14:23,:),[],1);
cTrange15p1200=Trange15p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange25p6=reshape(Trange25p6(14:23,:),[],1);
cTrange25p6=Trange25p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange25p60=reshape(Trange25p60(14:23,:),[],1);
cTrange25p60=Trange25p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange25p600=reshape(Trange25p600(14:23,:),[],1);
cTrange25p600=Trange25p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange25p1200=reshape(Trange25p1200(14:23,:),[],1);
cTrange25p1200=Trange25p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange35p6=reshape(Trange35p6(14:23,:),[],1);
cTrange35p6=Trange35p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange35p60=reshape(Trange35p60(14:23,:),[],1);
cTrange35p60=Trange35p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange35p600=reshape(Trange35p600(14:23,:),[],1);
cTrange35p600=Trange35p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange35p1200=reshape(Trange35p1200(14:23,:),[],1);
cTrange35p1200=Trange35p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange45p6=reshape(Trange45p6(14:23,:),[],1);
cTrange45p6=Trange45p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange45p60=reshape(Trange45p60(14:23,:),[],1);
cTrange45p60=Trange45p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange45p600=reshape(Trange45p600(14:23,:),[],1);
cTrange45p600=Trange45p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange45p1200=reshape(Trange45p1200(14:23,:),[],1);
cTrange45p1200=Trange45p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange60p6=reshape(Trange60p6(14:23,:),[],1);
cTrange60p6=Trange60p6(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange60p60=reshape(Trange60p60(14:23,:),[],1);
cTrange60p60=Trange60p60(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange60p600=reshape(Trange60p600(14:23,:),[],1);
cTrange60p600=Trange60p600(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);
Trange60p1200=reshape(Trange60p1200(14:23,:),[],1);
cTrange60p1200=Trange60p1200(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);

%find pressure as a function of topography
H=10000; %scale height on mars in meters

ctopo=topo(albedo>0.15 & topo>-6000 & topo<12000 & ti>150);

p0=zeros(size(ctopo));
p6=6*exp(-ctopo./H);
p60=60*exp(-ctopo./H);
p600=600*exp(-ctopo./H);
p1200=1200*exp(-ctopo./H);

temps=[cmaxT15p6;cmaxT15p60;cmaxT15p600;cmaxT15p1200;cmaxT25p6;cmaxT25p60;cmaxT25p600;cmaxT25p1200;cmaxT35p6;cmaxT35p60;cmaxT35p600;cmaxT35p1200;cmaxT45p6;cmaxT45p60;cmaxT45p600;cmaxT45p1200;cmaxT60p6;cmaxT60p60;cmaxT60p600;cmaxT60p1200];
green=[cgreen15p6;cgreen15p60;cgreen15p600;cgreen15p1200;cgreen25p6;cgreen25p60;cgreen25p600;cgreen25p1200;cgreen35p6;cgreen35p60;cgreen35p600;cgreen35p1200;cgreen45p6;cgreen45p60;cgreen45p600;cgreen45p1200;cgreen60p6;cgreen60p60;cgreen60p600;cgreen60p1200];
pressures=[p6;p60;p600;p1200;p6;p60;p600;p1200;p6;p60;p600;p1200;p6;p60;p600;p1200;p6;p60;p600;p1200];
delts=[cTrange15p6;cTrange15p60;cTrange15p600;cTrange15p1200;cTrange25p6;cTrange25p60;cTrange25p600;cTrange25p1200;cTrange35p6;cTrange35p60;cTrange35p600;cTrange35p1200;cTrange45p6;cTrange45p60;cTrange45p600;cTrange45p1200;cTrange60p6;cTrange60p60;cTrange60p600;cTrange60p1200];

psmall=[p6;p60;p600;p1200];
dsmall=[Trange35p6;Trange35p60;Trange35p600;Trange35p1200];
obl=[15*ones(1200,1);25*ones(1200,1);35*ones(1200,1);45*ones(1200,1);60*ones(1200,1)];

t4=temps.^4;
sqrtP=sqrt(pressures);

%from cftool using a+b*T^4*sqrt(P), rsquared=0.98
calgreen3=9.981+9.662.*10.^-10.*t4.*sqrtP;
res3=calgreen3-green;

%figure
%scatter(green,calgreen3,10,pressures)
%hold on
%plot(green,green,'r')
%hold off