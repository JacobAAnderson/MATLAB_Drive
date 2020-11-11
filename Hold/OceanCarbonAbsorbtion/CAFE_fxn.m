function [NPP,NPP_z,AP,AP_z2,zseq,Eksurf,KPURsurf,zeu]=CAFE_fxn(mld,sst,adg_443,chl,bbp_443,bbp_s,PAR,lat,yd,sal,aph_443)
%
%          mld,            #Mixed layer Depth [m]
%          sst,            #Sea surface temperature [Deg C]
%          adg_443,        #Absorption coefficient of detritus and gelbstoff [m-1]
%          chl,            #Cholorphyll a concentration [mg m-3]
%          bbp_443,        #Particulate backscaterring coefficient [m-1]
%          bbp_s,          #Slope of backscattering coefficient [dimensionless]
%          PAR,            #Daily integrated PAR [mol photons m-2 day-1]
%          lat,            #Latitude, nothern hemisphere is postive, decimal degrees
%          yd,             #Day of year
%          sal = 32.5,     #Salinity of water [PSU]
%          aph_443){       #Phytoplankton absorption coefficient [m-1]
%
%  #Notes on Input
%  # 1) If Mixed layer depth is not known, can pass deep value (>250 m) to assume
%  #    vertically homogenous water column
%  # 2) IOP data from GIOP-DC (Werdell et al. 2013)
%  #    - GIOP-DC assumes slope of adg = 0.018 [m-1]
%  #    - GIOP-DC assumes spectral shape of phyto absorption coefficient
%  #      is a function of Chl (Bricuad et al. 1998)
%  # 3) Have provision that if aph_443 is not passed, then calculate absorption coefficient
%  #    from Chl-a and Bricuad et al (1995)


%  #########################################################################################
%  # Step 1: Derive IOPs at 10 nm increments from 400 to 700 nm

%  #Variables
%  #wv: wavelength [nm]
%  #aw: Absorption of Pure Water (Pope and Fry 1997)
%  #aphi: Phytoplantkon absorption coefficient
%  #A.Bricaud, E.Bricaud: Spectral shape of phyto. abs. coeff (Bricaud et al. 1998)
%  #bbw:  Pure water backscattering coefficients derived from Zhang et al. (2009)

%  #Comments
%  #  GIOP-DC assumes slope of adg = 0.018 [m-1]
%  #  GIOP-DC assumes spectral shape of aphi is a function of Chl (Bricuad et al. 1998)
%  #########################################################################################

  wv=400:10:700;

  aw=[0.00663, 0.00473, 0.00454, 0.00495, 0.00635, 0.00922, 0.00979, 0.0106, 0.0127, ...
            0.015, 0.0204, 0.0325, 0.0409, 0.0434, 0.0474, 0.0565, 0.0619, 0.0695, 0.0896, ...
            0.1351, 0.2224,0.2644, 0.2755, 0.2916, 0.3108, 0.34, 0.41, 0.439, 0.465, 0.516, ...
            0.624];

  A_Bricaud=[0.0241, 0.0287, 0.0328, 0.0359, 0.0378, 0.0350, 0.0328, 0.0309, 0.0281, ...
                 0.0254, 0.0210, 0.0162, 0.0126, 0.0103, 0.0085, 0.0070, 0.0057, 0.0050, ...
                 0.0051, 0.0054, 0.0052, 0.0055, 0.0061, 0.0066, 0.0071, 0.0078, 0.0108, ...
                 0.0174, 0.0161, 0.0069, 0.0025];

  E_Bricaud=[0.6877, 0.6834, 0.6664, 0.6478, 0.6266, 0.5993, 0.5961, 0.5970, 0.5890, ...
                 0.6074, 0.6529, 0.7212, 0.7939, 0.8500, 0.9036, 0.9312, 0.9345, 0.9298, ...
                 0.8933, 0.8589, 0.8410, 0.8548, 0.8704, 0.8638, 0.8524, 0.8155, 0.8233, ...
                 0.8138, 0.8284, 0.9255, 1.0286];

  aphi = aph_443 .*  A_Bricaud .* chl.^(E_Bricaud) / (0.03711 .* chl.^(0.61479));

  a = aw + aphi + adg_443 .* exp( -0.018 .* (wv - 443));

  bbw=NaN.*ones(1,31);
  for w=1:31
    [t1, t2, t3, t4]=betasw_ZHH2009(wv(w), sal, sst,0.039); bbw(w)=t3./2;
  end
  
  bb=bbw + bbp_443 .* (443 ./ wv).^(bbp_s);

%  #########################################################################################
%  # Step 2: Calculate the fraction of photons absorbed by phytoplankton

%  # Assume %5 is reflected/upwelled from surface

%  # PARfraction derived from ASTMG173 reference spectrum
%  # http://rredc.nrel.gov/solar/spectra/am1.5/astmg173/astmg173.html
%  #########################################################################################

  PARfraction=[0.00227, 0.00218, 0.00239, 0.00189, 0.00297, 0.00348, 0.00345, 0.00344, ...
                   0.00373, 0.00377, 0.00362, 0.00364, 0.00360, 0.00367, 0.00354, 0.00368, ...
                   0.00354, 0.00357, 0.00363, 0.00332, 0.00358, 0.00357, 0.00359, 0.00340, ...
                   0.00350, 0.00332, 0.00342, 0.00347, 0.00342, 0.00290, 0.00314];

  absorbed_photons = PAR .* 0.95 * intgrt(wv',PARfraction' .* aphi'./a');
  
%  #########################################################################################
%  # Step 3: Derive Kd following Lee et al 2005, Eq. 11

%  #dec: declination of sun
%  #DL: Daylength calculated following Westberry et al 2008
%  #solzen: Solar zenith angle
%  #m0, m1, m2, m3: Coefficients from Lee et al. 2005
%  #########################################################################################

  lat=lat .* pi./180;
  dec= 23.5 .* cos((2 .* pi .* (yd - 172))./365) *pi./180;

  DL= -1 .* tan(lat) * tan(dec);
  DL(DL>1)= 1;                         %#Check for daylengths less than 0 hours
  DL(DL<(-1))= -1 ;                     %  #Check for daylengths greater than 24 hours
  DL = acos(DL) ./ pi;             %#Daylength in days
  
  solzen = 90 - asin (sin(lat) * sin(dec) - cos(lat) .* cos(dec) .*  cos(pi)) .*180./pi ;
  
  m0 = abs(1+0.005*solzen);
  m1 = 4.18; 
  m2 = 0.52; 
  m3 = 10.8;
  
  kd = m0 .* a + m1 .* (1 - m2 .* exp(-m3 .* a)) .* bb;

%  #########################################################################################
%  # Step 4: Construct water column irradiance through time (t) and depth (z)
%
%  #tseq:   Divides diurnal period into 101 increments [dimensionless]
%  #kdpar:  Attenuation coefficient of downwelling PAR [m-1]
%  #        Derived from Kd(490) following Morel et al. 2007.
%  #zeu:    Euphotic depth [m]. Taken as 0.1 mol photons m-2 day-1 isolume
%  #zseq:   Divides euphotic depth into 101 increments [m]
%  #delz:   Incremental change  in depth [m]
%  #E.tz:   Array of irradiance through time and depth [mol photons m-2 d-1]
%  #AP.tz:  Array of absorbed photons through time and depth [mol photons m-3 d-1]
%  #########################################################################################

  %tseq=0:0.01:1;
  tseq=0:0.02:1;
  kdpar = 0.0665 + 0.874 .* kd(10)  - 0.00121 ./ kd(10);
  zeu = -log (0.1 ./(PAR * 0.95)) ./ kdpar;
  zseq =linspace(0,ceil(zeu),101);
  delz = zseq(2) - zseq(1);

%  #Setup arrays
  %E_tz = zeros(101, 101); % #Irradiance
  %AP_tz = zeros(101, 101); % #Absorbed photons
  E_tz = zeros(101, 51); % #Irradiance
  AP_tz = zeros(101, 51); % #Absorbed photons

%  #Max surface PAR at solar noon
  PAR_noon = pi .* PAR .* 0.95 ./ 2 .* PARfraction;

  for t=1:51
  %for t=1:101
    for z=1:101
      E_tz(z,t) = intgrt(wv',PAR_noon' .* sin(pi .* tseq(t)) .* exp(-kd' .* zseq(z)));
      AP_tz(z,t) = intgrt(wv',PAR_noon' .* sin(pi .* tseq(t)) .* exp(-kd' .* zseq(z)) .* aphi');
    end
  end

%  #Integrate AP through time and depth
  AP_z = trapz(tseq,AP_tz,2);
  PAR_z = trapz(tseq',E_tz')';
  AP = intgrt(zseq',AP_z);
  
%  #Derive Upwelling Irradiance, absorbed energy is from Section 2
  Eu = absorbed_photons / AP ;

%  #Modify Irradiance and absorbed energy to account for upwelled irradiance
  E_tz = E_tz .* Eu;
  AP_tz = AP_tz .* Eu;

%  #########################################################################################
%  # Step 5: CALCULATE EK through depth

%  #Surface Ek is calculated following:
%  #I extended this algorithm to propagate Ek beneath the MLD

%  #IML: Median Mixed Layer Irradiance [mol photons m-2 hour-1]

%  #########################################################################################

  IML = (PAR * 0.95 ./ (DL.* 24)) .* exp(-0.5 .* kdpar * mld);
  Ek = 19 .* exp(0.038 .* (PAR .* 0.95 / (DL .* 24)) .^ 0.45 ./kdpar) .*ones(101,1);

  if mld < zeu

    Ek = Ek .* (1 + exp(-0.15 .* (PAR .* 0.95 ./ (DL .* 24)))) ./ (1 + exp(-3 .* IML));

%    #Find indices of zseq deeper than MLD
    deep = find(zseq > mld);
    Eg = (PAR ./ DL)  .* exp(-1.*kdpar.*zseq');
    Eg_mld = (PAR ./ DL)  .* exp(-1.*kdpar.*mld);
    Ek(deep) = 10 + (Ek(deep) - 10) ./ (Eg_mld - 0.1) .* (Eg(deep) -0.1);
  end

  Ek(Ek < 10) = 10;          %#Ensure Ek is no smaller than 110 umol m-2 s-1
  Ek = Ek * 0.0864;  %#Convert to mol photons/m2/day

%  #########################################################################################
%  # Step 6: CALCULATE SCF (Spectral Correction Factor)
%  # KPUR is spectrally scaled EK
%  #########################################################################################
  SCF = NaN.*ones(101,1);

  for z=1:101
    E_tzw = PAR_noon .* sin(pi .* tseq(50)) .* exp(-kd .* zseq(z));
    AQ_tzw = E_tzw .* aphi;
    SCF(z) = intgrt(wv',AQ_tzw') ./ (intgrt(wv',E_tzw') * mean(aphi)) /1.3;
  end

  KPUR = Ek./SCF;
  KPURsurf=KPUR(1);
%  #########################################################################################
%  # Step 7: Tie PHIMax to Ek
%  #########################################################################################

  phirange = [0.018, 0.030];
  ekrange = [150*86400/1e6, 10*86400/1e6];
  slope = (phirange(2) - phirange(1)) ./ (ekrange(2) - ekrange(1));
  phimax = phirange(2) + (Ek - ekrange(2)) .* slope;
  phimax(phimax < phirange(1)) = phirange(1);
  phimax(phimax > phirange(2)) = phirange(2);

%  #########################################################################################
%  # Step 8: Deep Chlorophyll Maxima is scaled to Ek
%  #########################################################################################

  if mld < zeu

    aphi_fact = ones(101,1);
    aphi_fact(deep) = 1 + Ek(1)./Ek(deep) .* 0.15;

%   Recalculate irradiance and absorbed energy over time depth array
    for i=1:length(wv)
        aphi_z(:,:,i)=aphi(i).*(aphi_fact * ones(1,51));
        adg_z(:,:,i)=adg_443.*exp(-0.018.*(wv(i)-443)) .* ones(101,51);
        aw_z(:,:,i)=aw(i) .*ones(101,51);
    end
    a_z=aphi_z + adg_z + aw_z;

    for i=1:length(wv), bb_z(:,:,i)=bb(i).*ones(101,51); end
    kd_z = m0.*a_z + m1 .*(1- m2.*exp(-m3.*a_z)) .*bb_z;
    
    E0_wv=PAR_noon' * sin(pi.*tseq); 
    par_t=trapz(wv',E0_wv);
    E_wv=E0_wv;
    for i=1:length(zseq)-1
        E_wv = E_wv .* exp(-squeeze(kd_z(i,:,:))' .* (zseq(i+1) - zseq(i)));
        AP_tz2(i,:) = trapz(wv',E_wv .* squeeze(aphi_z(i,:,:))');
        E_tz2(i,:) = trapz(wv',E_wv);
    end
    E_tz2=[par_t;E_tz2];
    AP_tz2=[trapz(wv',E0_wv .* squeeze(aphi_z(1,:,:))');AP_tz2];
    %#Modify absorbed energy to account for upwelled irradiance
    AP_tz2 =  AP_tz2 .* Eu;
  else
    for i=1:length(wv), aphi_z(:,:,i) = ones(101,51).* aphi(i); end
    AP_tz2=AP_tz; 
    E_tz2=E_tz;
  end


%  #########################################################################################
%  # Step 9: Final Calculations
%
%  #Calculate carbon specific growth Rate
%
%  #########################################################################################
  
  %Net Primary Production
  NPP_tz = AP_tz2 .* (phimax * ones(1,51)) .* tanh( (KPUR*ones(1,51)) ./E_tz) .* 12000;
  NPP_z = trapz(tseq, NPP_tz, 2);
  NPP   = trapz(zseq,NPP_z);

  %Absorbed photons
  AP_z2 = trapz(tseq,AP_tz2, 2);
  AP = trapz(zseq,AP_z2);
  
  Eksurf=Ek(1);
  %Fraction of absorbed photons dissipated
  %NPQ =AP_tz .* (phimax * ones(1,51)) .* 12000;
  %NPQ_z = trapz(tseq,NPQ, 2);
  %NPQ = 1 - NPP/trapz(zseq,NPQ_z);
  
  %Calculate NPP, AQ, and NPQ and mu in the mixed layer

  
%  if (mld > max(zseq))   %MLD is deeper than euphotic depth
%    NPP_mld = NPP;
%    AP_mld  = AP;
%    NPQ_mld = NPQ;
%  else                 %MLD is shallower than euphotic depth (i.e. vertical structure)
%    [foo,ind] =sort(abs(mld-zseq));
%    AP_mld = trapz(zseq(1:ind),AP_z(1:ind));
%    NPP_mld = trapz(zseq(1:ind),NPP_z(1:ind));
%  end

  %bbp_470 = bbp_443.*(443/470).^(bbp_s);
  %Cphyto = 12128 * bbp.470 + 0.59 ;
%  Cphyto = 13000 * (bbp_443 - 0.00035);
%  mu = NPP_mld ./ (Cphyto .* mld);

%  #####################################################
%  # Return Data
%  #####################################################
%
%  return(list(AP=AP,              #photons absorbed by phytoplankton [mol photons m-2 d-1]
%              AP.mld=AP.mld,          #photons absorbed by phytoplankton in the Mixed layer [mol photons m-2 d-1]
%              AP.sat = absorbed_photons,
%              NPP=NPP,            #Net Primary Production [mg C m-2 d-1]
%              NPP.mld=NPP.mld,        #Net Primary Production in the mixed layer [mg C m-2 d-1]
%              NPQ=NPQ,            #From of absorbed energy lost to non-photochemical quenching
%              Ek=Ek[1],           #Ek in the Surface
%              KPUR=KPUR,          #KPur (i.e. spectrally scaled Ek)
%              aphi.z=aphi.z[,1],  #Profile of phytoplantkon absorption coefficienct
%              zeu=zeu,            #Eupthotic Depth [m]
%              z=zseq,             #Depths over which data are evaluated
%              NPP.z=NPP.z,        #Daily volumetric net primary production [mg C m-3 d-1]  
%              Eu=Eu))             #Fraction of upwelled irradiance



