

function [X, Y, R, MODEL] = RANSAC(points, maxRadius, limit)

model = [];
X = [];
Y = [];
R = [];

if size(points,1) < 3, return, end

for ii = 1: max(size(points))
    
    indx = randi(max(size(points)),1, 3);
    
    [x1,y1,r1] = FitCircle(points(indx,:));
    
    if r1 > maxRadius || isnan(r1) || r1 < 1, continue, end
    
    test_model = NaN(size(points));
    cnt = 1;
    for jj = 1: max(size(points))
        
        test = [points(indx(1:2),:); points(jj,:)];
        [x2,y2,r2] = FitCircle(test);
        
        
        if isnan(r2) || abs(x1-x2) > limit || abs(y1-y2) > limit || abs(r1-r2) > limit
            continue
            
        else
            test_model(cnt,:) = points(jj,:);
            cnt = cnt + 1;
        end
        
    end
    
    if max(size(model)) < max(size(test_model))
        model = test_model(~isnan(test_model));
        X = x1;
        Y = y1;
        R = r1;
    end
    
end

if nargout > 3, MODEL = model; end

end


function [x,y,r] = FitCircle(points)

x = [];
y = [];
r = [];

if size(points,1) < 3, return, end

x1 = points(1, 1);
y1 = points(1, 2);

x2 = points(2, 1);
y2 = points(2, 2);

x3 = points(3, 1);
y3 = points(3, 2);

A = x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2;
B = (x1*x1 + y1*y1)*(y3-y2)       + (x2*x2 + y2*y2)*(y1-y3)       + (x3*x3 + y3*y3)*(y2-y1);
C = (x1*x1 + y1*y1)*(x2-x3)       + (x2*x2 + y2*y2)*(x3-x1)       + (x3*x3 + y3*y3)*(x1-x2);
D = (x1*x1 + y1*y1)*(x3*y2-x2*y3) + (x2*x2 + y2*y2)*(x1*y3-x3*y1) + (x3*x3 + y3*y3)*(x2*y1-x1*y2);

x = -B/(2*A);
y = -C/(2*A);
r = sqrt((B*B + C*C - 4*A*D)/(4*A*A));

end
