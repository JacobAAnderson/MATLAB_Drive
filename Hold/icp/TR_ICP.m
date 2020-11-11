


function [TT, err] = TR_ICP(fixed, moving)

iter = size(moving,2);
err = zeros(iter+1, 1);
TT = zeros(iter+1,2);

[match, mindist] = match_bruteForce(fixed, moving);

TT(1,:) = [0,0];
err(1) = sqrt(sum(mindist.^2)/length(mindist));   

for ii = 1: iter 
    
    tt =  moving(:,ii) - fixed(:,match(ii));
    tt(3) = 0;
   
    test = moving - tt;
    
    [~, mindist] = match_bruteForce(fixed, test);
    
    err(ii+1) = sqrt(sum(mindist.^2)/length(mindist));    
    TT(ii+1,:) = tt(1:2)';
end

end


function [match, mindist] = match_bruteForce(fixed, moving)
    m = size(moving,2);
    n = size(fixed,2);    
    match = zeros(1,m);
    mindist = zeros(1,m);
    for ki=1:m
        d=zeros(1,n);
        for ti=1:3
            d=d+(fixed(ti,:)-moving(ti,ki)).^2;
        end
        [mindist(ki),match(ki)]=min(d);
    end
    
    mindist = sqrt(mindist);
end


%% Junk!!!!!!!!!!!!!!!!!!!!!!!!!!!
%             new_xyz(3,:) = 0;
%             global_xyz(3,:) = 0;
%             moving = pointCloud(new_xyz',    'Color', repmat([0  0  1],size(new_xyz',1),1));
%             fixed =  pointCloud(global_xyz', 'Color', repmat([0  1  0],size(global_xyz',1),1));
%             [tform,movingReg, rms] = pcregistericp(moving,fixed, 'Tolerance',[0.01, 0.00] );
%             [TR, TT, ER] = icp(new_xyz, global_xyz);
%             eul = tform2eul(tform.T);            
%             fprintf( "RMS: %f\n", rms)
%             fprintf('\nTransform\n')
%             disp(tform.T)
%             disp(eul)

%             icp_xyz = movingReg.Location;
             
%             figure
%             plot3(global_xyz(1,:), global_xyz(2,:), global_xyz(3,:),'og')
%             hold on 
%             plot3(new_xyz(1,:), new_xyz(2,:), new_xyz(3,:),'*r')
%             plot3(icp_xyz(:,1), icp_xyz(:,2), icp_xyz(:,3),'*b')
%             hold off
%             view(0,90)

