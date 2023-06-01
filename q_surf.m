function [qsurf] = q_surf( e_vapor_press, varargin)
%UNTITLED3 Summary of this function goes here
%   given vapor pressure (at surface), e (mb)
%   given atmospheric pressure (at surface), P (mb)
%   calculate surface specific humidity

% optional inputs
% varargin{1}: atmopsheric pressure in mb

if nargin > 1 % if atmospheric pressure specified, then use it

atmos_press = varargin{1}; % specified atmospheric pressure, (mb)

else % if pressure not specified, then use 1013 mb as default value
atmos_press = 1013; % if atmospheric pressure not specified, then use default value of 1013 mb
end


qsurf = 0.622.*e_vapor_press./atmos_press;


end

