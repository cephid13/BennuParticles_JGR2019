function a_th = thermal_accel(r, A2M, phi, sun_dir, num_facets, nhat, cvec, AreaF, Tmap, Rsun)

% Inputs
%   r - particle position vector in Bennu-fixed frame [km]
%   A2M - particle area-to-mass ratio [km^2/kg]
%   phi - angle from midnight to Bennu prime meridian (xb) [deg]
%   sun_dir - unit vector AWAY from sun in body-fixed frame
%   num_facets - duh
%   nhat - facet normals in body-fixed frame [unit vectors]
%   cvec - facet centers in body-fixed frame [km]
%   AreaF - facet areas [km^2]
%   Tmap - Bennu temperature model from Rozitis (described below)
%   Rsun - current asteroid distance from sun [km]
%
% Outputs
%   a_th - thermal acceleration on particle in body-fixed frame [km/s^2]



% The other has the temperature per deg of local solar time on the x-axis, 360 columns starting with sun over the x-axis, with one row per facet. 

% Universal constants
c = 2.99792458e8; % [m/s]
Boltz = 5.670367e-8; % [W/m^2/K^4]
GR = 1368; % [W/m^2 at 1AU]
B = 2/3; 

% Bennu model parameters
emiss = 0.9; 
albedo = .044; % mean albedo
R_bennu = 134594309.932292; % distance from Bennu to sun at temperature file epoch of Jan 20th 00:00:00 [km]




T_ind = ceil(phi);
T0 = Tmap(:,T_ind);
T4 = (R_bennu/Rsun)^2.*T0.^4;


trp_all = zeros(3, num_facets);
for ii = 1:num_facets
    
    uvec_sc = cvec(:,ii) - r;
    u_sc = norm(uvec_sc);
    uhat_sc = uvec_sc./u_sc;
    
    % Simple facet to spacecraft visibility (no self shadowing)
    cosAlpha = dot(nhat(:,ii), -uhat_sc);
    if cosAlpha < 0
        continue; % H = 0, so doesn't contribute to TRP force
    end
    
    % Simple sun visibility (no self shadowing)
    cosTheta = dot(sun_dir, nhat(:,ii));
    if cosTheta > 0
        tau = 1;
    else
        tau = 0;
    end
    
    % Compute pressure term from this facet
    A_u = min(AreaF(ii)/u_sc^2, 1);
    Pj = (tau*albedo*GR*cosTheta + emiss*Boltz*T4(ii)) * (cosAlpha*A_u/c/pi);
    
    % compute f term
    fj = -uhat_sc;
    
    trp_all(:, ii) = Pj.*fj;
    
end

a_th = sum(trp_all,2).*(1+B)*A2M*1e6;

a_th = 1e-3.*a_th; % convert from [m/s^2] to [km/s^2] for output
    

end