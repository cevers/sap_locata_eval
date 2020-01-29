function delta = angdiff(x, y)
%This function is for internal use only. It may be removed in the future.

%ANGDIFF Find difference between two angles
%
%   DTHETA = angdiff(THETA1, THETA2) computes difference between angles 
%   THETA1 and THETA2 and returns the difference DTHETA within the interval 
%   [-pi pi]. Positive, odd multiples of pi map to pi and negative, odd 
%   multiples of pi map to -pi.

% Copyright 2014-2015 The MathWorks, Inc.

%#codegen


if nargin == 1
    % Single input means that x should be a numeric vector
    % Syntax: ANGDIFF(X)
    
    % Calculate difference
    d = diff(x);
else
    % There are two inputs
    % Syntax: ANGDIFF(X,Y)    
    
    d = y - x;
end

% Make sure the output is in the [-pi,pi) range
delta = wrapToPi(d);

end
