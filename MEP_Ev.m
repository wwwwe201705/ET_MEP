function [ GMEP, EMEP, HMEP, B,varargout ] = MEP_Ev( Rn, Ts, varargin )
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ GMEP, EMEP, HMEP, varargout ] = MEP_Module( Model_Type, Rn, Ts, varargin )
MEP_Module Maximum Entropy Production model of surface fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Variable definitions:

    

Rn: Net radiation (W/m^2), positive downward, column vector of time series observations
qs: specific humidity at the evaporating surface (kg/kg), column vector of time series observations
Ts: temperature of evaporating surface (Celsius), column vector of time series observations
z: theoretical height above surface at which Monin Obukhov Similarity Theory applies. Constant (scalar) value (meters)
I: media thermal inertial (TIU units), column vector of time series observations
stability: atmospheric stability characterization (binary variable), 1 = unstable, 0 = stable

soil_moisture: volumetric surface layer soil water content (m^3/m^3)
porosity: soil porosity (m^3/m^3)
b_param: parameter (exponent) for the saturation ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Required Input (depending on Model Type)

Rn, Ts, qs, I, z
[ GMEP, EMEP, HMEP, varargout ] = MEP_Module( Rn, Ts, qs, I, z, optional_inputs )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Optional Input
stability: atmospheric stability characterization (binary variable), 1 = unstable, 0 = stable

RnL: Net longwave radiation (W/m^2), positive downward, column vector of time series observations
        Note: Net Longwave is required for water-snow-ice case
air_density: air density (kg/m^3), column vector of time series observations or can be single value
air_pressure: air pressure (mb), column vector of time series observations, default is 1013 mb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
GMEP: Ground heat flux, positive downward (W/m^2)
EMEP: Latent heat flux, positive upward (W/m^2)
HMEP: Sensible heat flux, positive upward (W/m^2)

% Optional Outputs
B: Inverse Bowen ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


%% Convert Ts to Kelvin units
Ts_Celsius = Ts; % not converted (useful for q_surf and e_sat calculations)
Ts = Ts + 273.15;

%% Parameters

Model_Type = 1;

T0 = 273.15;
Tr = 300;
Cp = 1006;
Rv = 461;
Rd = 287.06;
a = 1;
g = 9.81;
k = 0.4;
Lv = 2.5E06;
Liv = 2.83E06;
alpha =1;
beta = 4.7;
r1 = 15;
r2 = 9;
C1 = [sqrt(3)/alpha, 2/(1+2*alpha)];
C2 = [r2/2, 2*beta];
e16 = 1/6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE INPUT PARSER

p = inputParser;

% default inputs
default_RnL = NaN;
default_air_density = 1.2; % kg/m^3
default_air_pressure = 1013; % mb
default_stability = NaN;
default_B_formula_type = 1; % B formula, if = 2, then use the alternative "negative" formula

% optional inputs
addParameter(p, 'stability',default_stability);
addParameter(p, 'RnL',default_RnL);
addParameter(p, 'air_density',default_air_density);
addParameter(p, 'air_pressure',default_air_pressure);
addParameter(p, 'B_formula_type',default_B_formula_type);

% default_target_folder = 'home_directory';
% default_Figure = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch Model_Type
    %%
    case 1 % Soil Surface Model #1, surface specific humidity given
        qs = varargin{1};
        I = varargin{2};
        z = varargin{3};
        
        addRequired(p,'Rn',@isnumeric);
        addRequired(p,'Ts',@isnumeric);
        addRequired(p,'qs',@isnumeric);
        addRequired(p,'I',@isnumeric);
        addRequired(p,'z',@isnumeric);
        
        % parse inputs (p, required inputs, optional inputs as varargin{:})
        parse(p, Rn, Ts, varargin{:});
        
        Rn = p.Results.Rn;
        Ts = p.Results.Ts;
        qs = p.Results.qs;
        I = p.Results.I;
        z = p.Results.z;
        stability = p.Results.stability;
        B_formula_type = p.Results.B_formula_type;
        
        % if air density is given as a scalar, then convert it to a vector
        air_density = p.Results.air_density;
        if length(air_density) == 1
            air_density = air_density.*(ones(length(Rn(:,1)), 1)); % convert air density to column vector
        end
        
        % if air pressure is given as a scalar, then convert it to a vector
        air_pressure = p.Results.air_pressure;
        if length(air_pressure) == 1
            air_pressure = air_pressure.*(ones(length(Rn(:,1)), 1)); % convert air pressure to column vector
        end
        
        % if soil thermal inertial is given as a scalar, then convert it to a vector
        if exist('I', 'var') == 1
            if length(I) == 1
                I = I .* (ones(length(Rn(:,1)), 1)); % convert I to column vector
            end
        end
        
        % if stability is given as a scalar, then convert it to a vector
        if exist('stability', 'var') == 1
            if length(stability) == 1 && isnan(stability) == 1 % if stability vector not provided,
                % then determine stability based on net radiation
                % if Rn +, assume daytime, unstable, therefore stability = 1
                % otherwise stable, stability = 0
                for i = 1:length(Rn)
                    if Rn(i) > 0
                        stability(i,1) = 1;
                    else
                        stability(i,1) = 0;
                    end
                end% next time step
            end % end stability vector length check
        end
        
        % Net Longwave Radiation
        RnL = p.Results.RnL;
        
        % Calculate Variables
        I0 = zeros(size(Rn));
        unsta = find(stability ==1);
        sta = find(stability ==0);
        I0(unsta) = air_density(unsta).*Cp.*sqrt(C1(1).*k.*z).*(C2(1).*k.*z.*g./(air_density(unsta).*Cp.*Tr)).^(e16);
        I0(sta) = air_density(sta).*Cp.*sqrt(C1(2).*k.*z).*(C2(2).*k.*z.*g./(air_density(sta).*Cp.*Tr)).^(e16);
        ice_num = find(Ts < T0); % ice surface
        Lv_vec = ones(size(Ts))*Lv;
        Lv_vec(ice_num) = Liv;
        
        % Calculate sigma and B
        sg = sqrt(a).*Lv_vec.^2.*qs./(Cp.*Rv.*Ts.^2);
        if B_formula_type == 1
            B = 6 .* (sqrt(1+11./36.*sg)-1);
        else
            B = -6 .* (sqrt(1+11./36.*sg)+1);
            
        end
        Bsr = B./sg;
        IdI0 = I./I0;
        
        % MEP Calculations
        for n = 1:1:length(Rn)
            if isnan(B(n)) == 0 && isnan(Rn(n)) == 0
                try % attempt to run fzero function
                    HMEP(n) = fzero(@(H) abs(H)^(e16)*(B(n)+1)*H + Bsr(n)*IdI0(n)*H - abs(H)^(e16)*Rn(n),0.5*Rn(n));
                    EMEP(n) = B(n)*HMEP(n);
                catch % if fzero fails, then set MEP fluxes to NaN
                    HMEP(n) = NaN;
                    EMEP(n) = NaN;
                end % end error catching
            else
                HMEP(n) = NaN;
                EMEP(n) = NaN;
            end
        end % next time step
        HMEP = HMEP';
        EMEP = EMEP';
        GMEP = Rn - HMEP - EMEP;
        
        
        
        
end % end switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

