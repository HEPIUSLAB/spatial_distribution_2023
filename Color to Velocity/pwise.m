% Function: pwise
%
% Purpose: Failed piecewise function for optimization (internal use)
%
% Input parameters:
%   b: double (parameters to optimize)
%   x: double (input color data)
%
% Output parameters:
%   v: double (output velocity)
%
% Created by: Denis Routkeitch (droutke1@jhmi.edu)

function v = pwise(b,x)
    
    v = (b(2)+b(3)*x(:,1)+b(4)*x(:,2)+b(5)*x(:,3))+max(b(1)-x(:,3), 0).*(b(6)+b(7)*x(:,1) + b(8)*x(:,2));

%     v = abs(min(b(10)-x(:,1), 0).*min(b(11)-x(:,2), 0).*min(b(1)-x(:,3), 0)).*(b(2)+b(3)*x(:,1)+b(4)*x(:,2)+b(5)*x(:,3))+...
%         (b(6)+b(7)*x(:,1) + b(8)*x(:,2) + b(9)*x(:,3));
    

%     v = (b(2)+b(3)*x(:,1)+b(4)*x(:,2)+b(5)*x(:,3))+ ...
%         max(b(1)-x(:,3), 0).*(b(6)+b(7)*x(:,1) + b(8)*x(:,2) + b(9)*x(:,2)) + ...
%         max(b(10)-x(:,3), 0).*(b(11)+b(12)*x(:,1) + b(13)*x(:,2) + b(14)*x(:,2));
end