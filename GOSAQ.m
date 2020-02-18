Ra = [1, 2];
Rd = [1, 2];
Pa = [-2, -3];
Pd = [-3, -2];
lam = 1.5;
K = 1;

options = optimoptions('fmincon','SpecifyObjectiveGradient',true);

L = -500;
U = 500;

while U - L > .001
    r = (L + U)/2;
    objf = @(y)obj(y, Ra, Rd, Pa, Pd, lam, r);
    boundsf = @(y)bounds(y, Ra, Pa, lam, K);
    [x,fval] = fmincon(objf, [.9,.9], [],[],[],[],[],[], boundsf, options)
    if fval > 0
        L = r;
    else 
        U = r;
    end
end


function [fx, gradx] = obj(y, Ra,Rd,Pa,Pd,lam,r)
    fx = 0;
    targets = length(y);
    for i = 1:targets
        fx = fx + (y(i) * exp(lam * Ra(i)) * (r - Pd(i)));
        fx = fx + ((((Rd(i) - Pd(i)) * exp(lam * Ra(i))) / (lam * (Ra(i) - Pa(i)))) * (y(i) * log(y(i))));
    end
    for i = 1:targets
        gradx(i) = 0;
        gradx(i) = gradx(i) + (exp(lam * Ra(i)) * (r - Pd(i)));
        gradx(i) = gradx(i) + ((((Rd(i) - Pd(i)) * exp(lam * Ra(i))) / (lam * (Ra(i) - Pa(i)))) * (log(y(i))));
        gradx(i) = gradx(i) + (((Rd(i) - Pd(i)) * exp(lam * Ra(i))) / (lam * (Ra(i) - Pa(i))));
    end
end

function [c,ceq] = bounds(y, Ra, Pa, lam, K)   
    targets = length(y);
    c(1) = -K;
    for i = 1:targets
        c(1) = c(1) + (log(y(i))/(lam * (Pa(i) - Ra(i))));
    end
    for i = 1:targets
        c(i + 1) = y(i) - 1;
    end
    for i = 1:targets
        c(1 + targets + i) = -y(i) + exp(lam * (Pa(i) - Ra(i)));
    end
    ceq = [];
end
    
    


