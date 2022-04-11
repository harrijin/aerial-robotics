function [eak] = voltageConverter(Fk,NBk,P)
% voltageConverter : Generates output voltages appropriate for desired
%                    torque and thrust.
%
%
% INPUTS
%
% Fk --------- Commanded total thrust at time tk, in Newtons.
%
% NBk -------- Commanded 3x1 torque expressed in the body frame at time tk, in
%              N-m.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%
% OUTPUTS
%
% eak -------- Commanded 4x1 voltage vector to be applied at time tk, in
%              volts. eak(i) is the voltage for the ith motor.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Harrison Jin
%+==============================================================================+  
% Get quadParams
r = P.quadParams.rotor_loc;
cm = P.quadParams.cm(1);
eamax = P.quadParams.eamax;
kF = P.quadParams.kF;
kN = P.quadParams.kN;
kT = kN(1)/kF(1);
% Set Fmax
Fmax = kF(1)*(cm*eamax)^2;
% Initialize conversion matrix
conversionMat = [1 1 1 1; ...
                 r(2,1) r(2,2) r(2,3) r(2,4); ...
                 -r(1,1) -r(1,2) -r(1,3) -r(1,4); ...
                 -kT kT -kT kT];
% Loop until Fi is zero or positive and less than or equal to Fmax
alpha = 1;
beta = 0.9;
F = (Fmax+1) * ones(4,1);
while (max(F) > Fmax && alpha > 0)
    F = conversionMat\[min(Fk, 4*beta*Fmax);alpha*NBk];
    alpha = alpha - 0.01;
end
F(F<0) = 0;
% Convert F to voltages
eak = sqrt(F./kF)/cm;
end