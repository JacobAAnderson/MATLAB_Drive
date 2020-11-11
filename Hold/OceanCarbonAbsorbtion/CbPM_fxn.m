function [ppz, mu, chlz, Cz, PARz, Kz, Zeu,Ig]=CbPM_fxn(par,K490,mld,zno3,dl,chl,bbp)
%%%       this will be a first crack at calculating growth rate profiles
%%%       products needed: Chl, bbp, PAR, K490, MLD, Zno3,
%%%       daylength
%%%
%%%       Dependencies:  austinPetzold_1986.m, intgrt.m (can use trapz.m from MATLAB)
%%%                      morel1988_Kd.m,
%%% -- HARDWIRED STUFF -- %%%
ppz=[]; mu=[]; chlz=[]; Cz=[]; PARz=[]; Kz=[]; Zeu=[];
y0 = 0.0003;
z = [1:300]';
lambda = [400 412 443 490 510 555 625 670 700];
parFraction = [0.0029 0.0032 0.0035 0.0037 0.0037 0.0036 0.0032 0.0030 0.0024];        %%% in energy fractions
%%%

%%% -- MAKE SPECTRAL ED BASED ON PROSCRIBED FRACTIONS ABOVE -- %%%

%%% -- calculate median PAR in MLD at each pixel
Ed0 = par .* parFraction;
Klambda = austinPetzold_1986(lambda, K490);
Ez_mld = Ed0(:).*0.975.*exp(-Klambda(:).*mld./2);
PAR_mld = intgrt(lambda',Ez_mld)./dl;
%%% -- ESTIMATE CHL:C FROM SATELLITE AND CORRESPONDING NUTRIENT STRESS (i.e., difference from max) -- %%%
C = (bbp-0.00035) .* 13000;          %%% from B et al. (2005)
chlC_sat = chl ./ C;
delCHLC = [0.022 + (0.045 - 0.022).*exp(-3.*PAR_mld)] - chlC_sat;

%%% -- CALCULATE Kd OFFSET WHICH WILL BE CARRIED THROUGH TO DEPTH AND PICKS UP CONTRIBUTION FROM NON-CHL ATTEN. -- %%%
[Kd,wave] = morel1988_Kd(chl); 
Kd = interp1(wave,Kd,lambda');
Kdif = Klambda(:) - Kd;

%%% -- stuff for estimating Cz profile -- %%%
R = 0.1; 
scalar = 1./R;

%%% -- NOW ITERATE THROUGH WATER COLUMN AND RECONSTRUCT LIGHT FIELD, CALCULATE Kd, CHL:C, AND MU -- %%%
  for m=1:length(z)
      if z(m)<mld || m == 1
         [Kd,wave] = morel1988_Kd(chl); 
         Kd=interp1(wave,Kd,lambda');
         chl_C(m) = chlC_sat; 
         chlz(m) = chl_C(m) .* C;
         mu(m) = 2 .* (chlC_sat-y0)./([0.022 + (0.045 - 0.022).*exp(-3.*PAR_mld)]-y0) .* [1 - exp(-5.*PAR_mld)];
         Ezlambda(:,m) = Ed0.*0.975.*exp(-Klambda.*z(m));
         PARz(m) = intgrt(lambda',Ezlambda(:,m));
         Kz(m,:) = Klambda;
         prcnt(m) = real(PARz(m)) ./ (par.*0.975);
         Cz(m)=C;
      else
         %%% -- IF BELOW MIXED LAYER MUST TREAT PROPERTIES DIFFERENTLY -- %%%
         [Kd,wave] = morel1988_Kd( chlz(m-1) ); 
         Kd=interp1(wave,Kd,lambda');
         Kd = Kd + Kdif;
         Kz(m,:)=Kd;
         Ezlambda(:,m) = Ezlambda(:,m-1).*exp(-Kd.*1);    %%% -- basically, it calculates the light level based on that
                                                          %%% -- just above it less 1 meter with the new Kd(lambda)
         PARz(m) = intgrt(lambda',Ezlambda(:,m));
         prcnt(m) = real(PARz(m)) ./ (par.*0.975);
           if mu(m-1) >= R
              Cz(m) = C;
           else
              Cz(m) = C.*[1-scalar.*(R-mu(m-1))];
           end
         deltaZ = zno3-z(m);
         deltaZ(deltaZ<0)=0;
         chl_C(m) = [0.022 + (0.045-0.022).*exp(-3.*PARz(m)./dl)] - (delCHLC .* (1-exp(-0.075.*deltaZ)'));
         chlz(m) = chl_C(m) .* Cz(m);
         mu(m) = 2 .* (chlC_sat-y0)./([0.022 + (0.045 - 0.022).*exp(-3.*PAR_mld)]-y0) .* [1 - exp(-5.*PARz(m)./dl)];
      end
  end %%% -- end depth loop

%%% -- AND FINALLY COMPUTE NPP -- %%%
ind=~isnan(prcnt);
if ind==1, Zeu = interp1(prcnt(ind),z(ind),.01); end
ppz = mu .* Cz;
Ig=PAR_mld;
return
