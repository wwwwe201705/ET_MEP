function [ e_s ] = e_sat( temperature_meas, varargin )
%UNTITLED4 Summary of this function goes here
% Version 02: Temperature input is in Celsius
%   Calculates saturation vapor pressure (mb) for a given temperature (Celsius) (version 2 edit)
%	Allows for vector input
%	Allows for specification of water vaporization or ice sublimation

% optional input
% varargin('model type', default = water, if varargin{1} == 'ice' then change Lv to Lv of ice, string input, must say 'water' or 'ice' if used
% if not specified, then use water, Lv = 2.5E06

if nargin > 1 % if varargin is given (meaning the number of input arguments is greater than 1)
    
    if isequal(varargin{1},'water') == 1
        Lv = 2.5E06; %latent heat of vaporization (water)
    end
    
    if isequal(varargin{1},'ice') == 1
        Lv = 2.83E06; % latent heat of sublimation (ice)
    end
    
else % if varargin not specified
    Lv = 2.5E06; %latent heat of vaporization (water)
end

Rv = 461;
To = 273;
T = temperature_meas; % make sure temperature measurement is in Celsius!
T = T + 273.15; % convert to Kelvin

e_s = 6.11.*exp((Lv/Rv)*(1/To - 1./T));


end

