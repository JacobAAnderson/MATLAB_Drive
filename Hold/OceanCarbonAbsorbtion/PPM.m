function NPP = PPM(bbp, chl, irr, k490, MLD, dl)
    
  %________________________________________________________________________
  % 2) Declare Spectral Values
  lambda      = [ 400, 412, 443, 490, 510, 555, 625, 670, 700];
  parFraction = [ 0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024];
  X           = [ .11748, .122858, .107212, .07242, .05943, .03996, .04000, .05150, .03000];
  e           = [.64358, .653270, .673358, .68955, .68567, .64204, .64700, .69500, .60000]; 
  Kw          = [.01042, .007932, .009480, .01660, .03385, .06053, .28400, .43946, .62438];
  
  
  %________________________________________________________________________
  % 3) Initalize Values
  
  y0   = 0.0003;
  z    = 0:200;
  r    = 0.1;
  uMax = 2.0;
  
  Klambda = NaN(length(chl), 9);
  E       = NaN(length(chl), 9);
  Kdif    = NaN(length(chl), 9);
  PAR     = zeros(length(chl));
  
  
  % Compute median PAR in MLD
  
%   for i = 1:9
%     Klambda(:,i) = austinPetzold_1986(lambda(i),k490);
%     E(:,i)       = irr * parFraction(i) * 0.975 * exp(-Klambda(:,i) * MLD / 2.0);
%   end
%   
%   for i = 1:8, PAR = PAR + (lambda(i+1) - lambda(i)) * (E(:,i) + E(:,i+1)) / 2; end
  
  % PAR = PAR / dl  # Changed from above
  irr    = irr / dl;
  PAR_kd = irr^0.45/k490;
  
  % PPM model coefficients
  a = -16.8;
  b = 1.55;
  c = 47.027;
  d = 0.0125;
  
  % Photoacclimation terms 
  dm            = 19 *exp(0.038*PAR_kd);
  sm            = (1+exp(-0.15*irr))/(1+exp(-3*PAR));
  sm(sm < 1)    = 1;
  dm(dm > 1000) = 1000;
  theta         = dm*sm;  
  
  % use either bbp or model C
  carbon = 13000 * (bbp - 0.00035);
  % carbon = dm * sm * chl; 
  
  % Light limited function for depth 
  IgFunc = 1- exp(-5 * irr);
  
  % Light limitation function  
  lightmu = I(1 / dm * a) + b; 
  
  % Nutrient limitation function  
  nutmu = I(1 / (dm*sm) * c) + d; 
  
  % Estimate phytoplankton growth 
  mu = lightmu * nutmu *  IgFunc;
  mu(mu > uMax) = uMax;
  
  %Need to find out E just above MLD for later on
%   for i = 1:9 
%     Kdif(:,i) = Klambda(:,i) - (Kw(i) +  X(i) * chl^e(i));
%     E(:,i)     = irr * parFraction(i) * 0.975 * exp(-Klambda(:,i) * floor(MLD));
%   end
  
  %________________________________________________________________________
  % 3) Compute NPP through depth
  NPP = NaN(length(irr), 201);
  
  % Parameters dependent on level z-1 ?????????????????????
%   mu.z  = mu;
%   chl.z = chl;
%   E.z   = E;
  
  % save(mu, file = paste(getwd(), "/Model_output/mu/", month,".RData", sep=""))  
  % save(carbon, file = paste(getwd(), "/Model_output/carbon/", month,".RData", sep=""))  
  % save(MLD, file = paste(getwd(), "/Model_output/MLD/", month,".RData", sep=""))  
  
  % Fill in NPP where data is present. Will overwrite NPP for  data beneath the MLD
  for k = 1:201, NPP(:,k) = mu * carbon; end 
  
  % Start at shallowest MLD  
  surface = ceiling(min(MLD))+1;
  
  % rm(Klambda);rm(Kdif);
  % rm(E); rm(E.z)
  % rm(bbp);rm(chl);rm(dl);rm(theta)
  % gc()
  
  
  for k = 1:201
   
    % k= 10 
    % deep <- which(MLD<=z[k], arr.ind = T)
    
%     irr.z  = irr * exp(-k490*z(k));
    z  = irr * exp(-k490*z(k));
    
    IgFunc = 1 - exp(-5 * irr.z);
    
    lightmu(k) = I(1 / dm(k) * a) + b; 
    nutmu(k)   = I(1 / (dm(k)*sm(k)) * c) + d; 
    
    mu(k) = lightmu(k) * nutmu(k) *  IgFunc(k);
    
    mu(k, mu(k)>uMax) = uMax;
    
    ind <- which(mu.z <  r & MLD<=z[k], arr.ind=T)

    % chl[k]    = chl_c[k] * carbon[k]
    NPP(:,k) = round(mu(k) * carbon,3);
    
    % Parameters dependent on level z-1
    mu.z(k) = mu(k);
    % chl.z[k]  <- chl[k]
    % NPP.z[,k]  = round(mu.z[k] * carbon,3)
    % for l = 1:9, E.z[k,l] = E[k,l]; end
    
    print(k)
    
  end
  
  
  % Need to parse between above and below the mixed layer 

  % Remove some data if necessary 
  % rm(E.z); rm(deltaZ); rm(IgFunc);
  % rm(chl.z); rm(mu.z); rm(nutTempFunc);
  
  NPP.surface = NaN(2160,4320);   
  NPP.all     = NaN(2160,4320); 
  NPP.deep    = NaN(2160,4320);
  
  % For integration, truncate MLD to 200 m
  MLD(MLD>200) = 200;
  
  for i = 1:length(irr)
    % Integrate
    fn <- splinefun(z,NPP[i,])       
    try
      NPP.surface(ocean(i,1),ocean(i,2)) = round(integrate(fn,0,MLD(i)));
    catch
        NPP.all(ocean(i,1),ocean(i,2)) =  round(integrate(fn,0,200));
    end
  end
  
  
  % NPP.all   <- NPP.all[ocean]  
  % NPP.surface   <- NPP.surface[ocean]  
  NPP.deep = NPP.all - NPP.surface;
  
  % Save Rdata
  % save(NPP.surface, file = paste(getwd(), "/Model_output/NPP Surface/", month,".RData", sep="")) 
  % save(NPP.deep, file = paste(getwd(), "/Model_output/NPP Deep/", month,".RData", sep=""))  
  % save(NPP.all, file = paste(getwd(), "/Model_output/NPP/", month,".RData", sep=""))  
 
end  