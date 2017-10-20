function [ netiM,netiD,netiP,timediff,lumfraction,indexD ] = insolation_sloped( obl,ecc,Lp,time,iM,iP,iD,smallf )
%Use illumination lookup table to calculate illumination history for an obliquity track, then use
%that and solar history to get insolation history.

%INUPTS:
%   obl = obliquity history
%   ecc = eccentricity history
%   Lp = longitude of perihelion history
%   time = time of each point in the obliquity history
%   iM = mean illumination lookup table
%   iP = peak illumination loookup table
%   iD = warmest 4-hour period illumination lookup table
%   smallf = lookup table for sun's luminosity as a function of time

%OUTPUTS:
%   netiM = mean insolation for each point on the obliquity track
%   netiP = peak insolation for each point on the obliquity track
%   netiD = insolation during warmest 4-hour period for each point on the
%   obliquity track
%   timediff = time in Gyr before today for each point on the track

%round values to search lookup table
oblindex=round(obl./5+1); %obliquity
eccindex=round(ecc./0.01+1); %eccentricity
Lpindex=round(Lp./10+1); %longitude of perihelion
toobig=Lpindex>36;
Lpindex(toobig)=ones(size(Lpindex(toobig)));

%want to investigate 20 degree sloping surface at 20S.
iMfrac=zeros(size(ecc,1),1); %each point, 1 latitude
iPfrac=zeros(size(ecc,1),1);
iDfrac=zeros(size(ecc,1),1);
indexM=zeros(size(ecc,1),1);
indexP=zeros(size(ecc,1),1);
indexD=zeros(size(ecc,1),1);
for ind1=1:size(ecc,1)
        iMrow1=iM(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,1);
        [iMfrac(ind1,1),indexM(ind1,1)]=max(iMrow1);
        iProw1=iP(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,1);
        [iPfrac(ind1,1),indexP(ind1,1)]=max(iProw1);
        iDrow1=iD(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,1);
        [iDfrac(ind1,1),indexD(ind1,1)]=max(iDrow1);
end

%now that we have illuminations, convert to insolations
%Sun is currently 4.5 Gyr old
timeingyr=time./(10^9);
today=timeingyr(end);
timediff=abs(timeingyr-today);
timesincesunform=4.56-timediff;
timeindex=round(timesincesunform./0.000001+1);

%find fraction of sun's current solar constant that it was at that time
lumfraction=smallf(timeindex);

%calculate insolation
solarconstant=1367.7; %in W/m^2
netiM=solarconstant.*lumfraction.*iMfrac;
netiD=solarconstant.*lumfraction.*iDfrac;
netiP=solarconstant.*lumfraction.*iPfrac;

end

