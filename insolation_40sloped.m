function [ netiM,netiD,netiP,timediff,lumfraction,iDfrac ] = insolation_40sloped( obl,ecc,Lp,time,iM,iP,iD,smallf )
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

%want to investigate surface at 40 degrees south with different slopes
iMfrac=zeros(size(ecc,1),3); %each point, 3 slopes
iPfrac=zeros(size(ecc,1),3);
iDfrac=zeros(size(ecc,1),3);
for ind1=1:size(ecc,1)
        iMrow1=iM(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,1);
        iMrow2=iM(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,2);
        iMrow3=iM(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,3);
        iMfrac(ind1,1)=max(iMrow1);
        iMfrac(ind1,2)=max(iMrow2);
        iMfrac(ind1,3)=max(iMrow3);
        iProw1=iP(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,1);
        iProw2=iP(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,2);
        iProw3=iP(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,3);
        iPfrac(ind1,1)=max(iProw1);
        iPfrac(ind1,2)=max(iProw2);
        iPfrac(ind1,3)=max(iProw3);
        iDrow1=iD(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,1);
        iDrow2=iD(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,2);
        iDrow3=iD(oblindex(ind1),eccindex(ind1),Lpindex(ind1),:,3);
        iDfrac(ind1,1)=max(iDrow1);
        iDfrac(ind1,2)=max(iDrow2);
        iDfrac(ind1,3)=max(iDrow3);
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
netiM=zeros(size(iMfrac));
netiD=zeros(size(iDfrac));
netiP=zeros(size(iPfrac));
netiM(:,1)=solarconstant.*lumfraction.*iMfrac(:,1);
netiM(:,2)=solarconstant.*lumfraction.*iMfrac(:,2);
netiM(:,3)=solarconstant.*lumfraction.*iMfrac(:,3);
netiD(:,1)=solarconstant.*lumfraction.*iDfrac(:,1);
netiD(:,2)=solarconstant.*lumfraction.*iDfrac(:,2);
netiD(:,3)=solarconstant.*lumfraction.*iDfrac(:,3);
netiP(:,1)=solarconstant.*lumfraction.*iPfrac(:,1);
netiP(:,2)=solarconstant.*lumfraction.*iPfrac(:,2);
netiP(:,3)=solarconstant.*lumfraction.*iPfrac(:,3);

end

