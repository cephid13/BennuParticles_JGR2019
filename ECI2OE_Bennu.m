function [oe_out, radius, hhat, edotArgX] = ECI2OE_Bennu(state_hist, T, et0str)
% ECI2OE 
%
% Usages:   oe_out = ECI2OE(pos,vel,varargin)
%
% Description:
%   This function follows the outline from Chobotov's "Orbital Mechanics"
%   page 66 to compute the orbital elements starting with position and 
%   velocity.
%
% Inputs:
%   pos - 3x1 ECI position vector
%   vel - 3x1 ECI velocity vector
%   varargin - (optional) a flag, if set to a non-zero value uses ft,
%              otherwise defaults to km
% Outputs:
%   oe_out - 6x1 vector of orbit elements = [a, e, i, w, Om, nu]. Units of
%            'a' match the varargin choice (default km). All angles are in
%            rads. i in [0 pi], others in [0 2pi]
%
% Modification history:
%   7/25/2012 - Jay - created
%

rho = 1.19e12; % [kg/km^3]
Vol = 0.06165181123859803;                   % [km^3]     Asteroid's Volume from shape model header
G           =  6.67384e-20;             % [km^3/(kg*s^2)] Gravitational constant
M = rho*Vol;        % Mass [kg]
mu_ast = G*M;

oe_out = zeros(size(state_hist));
radius = zeros(size(state_hist,1),1);
hhat = zeros(size(state_hist,1),3);
edotArgX = zeros(size(state_hist,1),1);

cspice_furnsh( { 'naif0012.tls', 'pck00010.tpc', 'de424.bsp', 'bennu_draft_v12.tpc', 'orx_v13.tf', 'orx_181231_190305_190211_od098-N-M0D-P-M1D_v1.bsp'} );
et0 = cspice_str2et( {et0str} );

SN_2_J2000_epoch = cspice_pxform('ORX_BENNU_SUN_NORTH','J2000',et0);
uhat_epoch = SN_2_J2000_epoch(:,1);
uhat_epoch = uhat_epoch./norm(uhat_epoch);

for ii = 1:size(state_hist,1)

    SN_2_J2000 = cspice_pxform('ORX_BENNU_SUN_NORTH','J2000',et0+T(ii));
    uhat = SN_2_J2000(:,1);
    uhat = uhat./norm(uhat);
    ang_sun_move = acos(dot(uhat_epoch,uhat));
    if ~isreal(ang_sun_move)
        ang_sun_move = 0;
    end
    
    epochF_2_currF = [cos(ang_sun_move) -sin(ang_sun_move) 0; sin(ang_sun_move) cos(ang_sun_move) 0; 0 0 1];
    
    
    r = epochF_2_currF*state_hist(ii,1:3)';
    rmag = norm(r);
    radius(ii) = rmag;
%     v = epochF_2_currF*(state_hist(ii,4:6)' - cross([0;0;-pi/180/24/3600],state_hist(ii,1:3)')); % tecnically right, but unnecessary for slow heliocentric orbit
    v = epochF_2_currF*state_hist(ii,4:6)';
    vmag = norm(v);
    
    % Semi-major axis [ft]
    a = 1/((2/rmag)-(vmag^2/mu_ast));
    
    What = cross(r,v)/norm(cross(r,v));
    
    hhat(ii,:) = What; % Xhat component of hhat
    
    % Inclination [rad]
    i = acos(dot(What,[0;0;1]));
    
    evec = (1/mu_ast)*((vmag^2-(mu_ast/rmag)).*r - dot(r,v).*v);
    
    % Eccentricity
    e = norm(evec);
    
    tempVec = cross(What,evec)./e;
    edotArgX(ii) = tempVec(1);
    
    Nhat = cross([0;0;1],What)/norm(cross([0;0;1],What));
    
    % Right Ascension of the Ascending Node [rad]
    craan = dot([1;0;0],Nhat);
    sraan = dot(cross([1;0;0],Nhat),[0;0;1]);
%     if ~isreal(craan) || ~isreal(sraan)
%         disp('what?');
%     end
    raan = atan2(sraan,craan);
    if(raan < 0)
        raan = raan + 2*pi;
    end
    
    % Argument of Perigee [rad]
    cw = dot(Nhat,evec)/e;
    sw = dot(cross(Nhat,evec)./e,What);
    w = atan2(sw,cw);
    if(w < 0)
        w = w + 2*pi;
    end
    
    % True Anomaly [rad]
    ctheta = dot(evec,r)/(e*rmag);
    stheta = dot((cross(evec,r)/(e*rmag)),What);
    theta = atan2(stheta,ctheta);
    if(theta < 0)
        theta = theta + 2*pi;
    end
    
    oe_out(ii,:) = [a e i w raan theta];

end

cspice_kclear;

end % end of orb_elements