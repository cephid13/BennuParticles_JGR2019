function state_hist_Sun_rot = ACI2SunRot_Bennu(state_hist, T, et0str)
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


cspice_furnsh( { 'naif0012.tls', 'pck00010.tpc', 'de424.bsp', 'bennu_draft_v12.tpc', 'orx_v13.tf', 'orx_181231_190305_190211_od098-N-M0D-P-M1D_v1.bsp'} );
et0 = cspice_str2et( {et0str} );

SN_2_J2000_epoch = cspice_pxform('ORX_BENNU_SUN_NORTH','J2000',et0);
uhat_epoch = SN_2_J2000_epoch(:,1);

state_hist_Sun_rot = zeros(size(state_hist));

for ii = 1:size(state_hist,1)

    SN_2_J2000 = cspice_pxform('ORX_BENNU_SUN_NORTH','J2000',et0+T(ii));
    uhat = SN_2_J2000(:,1);
    ang_sun_move = acos(dot(uhat_epoch,uhat));
    
    epochF_2_currF = [cos(ang_sun_move) -sin(ang_sun_move) 0; sin(ang_sun_move) cos(ang_sun_move) 0; 0 0 1];
    
%     if abs(det(epochF_2_currF)-1) > 1e-12
%         disp('hi')
%     end
    
    state_hist_Sun_rot(ii,1:3) = (epochF_2_currF*state_hist(ii,1:3)')';
    state_hist_Sun_rot(ii,4:6) = (epochF_2_currF*state_hist(ii,4:6)')';

end

cspice_kclear;

end % end of orb_elements