clc


distances = [norm([1,4]-[2,1])
             norm([1,4]-[4,2])
             norm([2,1]-[4,2])];


f = [ -4;  3;  -1];

A = [];

b =  [];

Aeq = [];
beq = [];

intcon = [1,2,3];

lb = zeros(3,1);
ub = ones(3,1); 

x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);

disp(x)




