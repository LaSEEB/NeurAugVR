function [logans] = cmp(a,b)

   logans = abs(a-b)<1e-9;
   
end