function [T_all, X_all] = particles_thermal_extend(inFile)
% This function creates a simiple simulation that shows what happens to a
% population of spherical particles that are lofted off the surface of a
% spun-up asteroid. We use a simple oblate spheroid gravity field (only
% second order - needs upgrading) and spherical constant density particles
% that don't interact gravitationally with one another. SRP is present to
% see the effect of size sorting in the lofted material.

% Due to SRP, probably easiest to propagate this in a Sun-Asteroid rotating
% frame to start; X-away from Sun, Z-along the asteroid spin vector

% Need to add:
%   - Shadowing for SRP
%   - polyhedral gravity model

load(inFile); % containts state_init, srpInfo, caseNumber, caseStr

% Parent body parameters
% based roughly on Bennu
% oblate ellipsoid, a=b at equator
rho = 1.19e12; % [kg/km^3]
R_a = 0.272; % [km]
R_c = 0.250; % [km]
Vol = 0.06165181123859803;                   % [km^3]     Asteroid's Volume from shape model header
G           =  6.67384e-20;             % [km^3/(kg*s^2)] Gravitational constant
M = rho*Vol;        % Mass [kg]
mu_ast = G*M;           % gravitational constant [km^3/s^2]
I_s = (1/5)*M*(R_a^2 + R_c^2); % [kg km^2]
I_z = (1/5)*M*(R_a^2 + R_a^2);
w_spin = (2*pi/(4.297461*3600)); % Bennu spin rate
% w_spin = (2*pi/(2.55*3600)); % increased spin rate to loft particles
pole_vec = [0; 0; 1];
% a_helio = 0.9*149597870;    % 0.9 AU in [km]
G1 = 1e14;  % [kg km / s^2]

% Bennu orbit info
cspice_furnsh( { 'naif0012.tls', 'pck00010.tpc', 'de424.bsp', 'bennu_draft_v12.tpc', 'orx_v13.tf', 'orx_181231_190305_190211_od098-N-M0D-P-M1D_v1.bsp'} );
Bennu_id = '2101955';
ets = cspice_str2et( {'2019 JAN 06 20:58:28.000 UTC', '2019 Jan 19 00:53:40.900 UTC', '2019 Feb 11 23:27:28.480 UTC'} );

% start time is 40 days after epoch for this case
ets = ets + 40*24*3600.*ones(size(ets));

if caseNumber == 1 % Jan 6, Near
    et_t0 = ets(1);
elseif caseNumber == 2 % Jan 6, Far
    et_t0 = ets(1);
elseif caseNumber == 3 % Old Jan 19
    et_t0 = ets(2);
elseif caseNumber == 4 % Jan 19, High res
    et_t0 = ets(2);
elseif caseNumber == 5 % Feb 11, High res
    et_t0 = ets(3);
end

SN_2_J2000_epoch = cspice_pxform('ORX_BENNU_SUN_NORTH','J2000',et_t0);
uhat_epoch = SN_2_J2000_epoch(:,1);

[posBennu,~]= cspice_spkpos( num2str(Bennu_id), et_t0, 'ORX_BENNU_ORBIT_FIXED', 'NONE', 'SUN' );
a_helio = posBennu(1);

% Want theta0 from SPICE - just take as the angle between Bennu's body
% frame and the Sun-North frame at t0
% then we have negative signs in computing theta0 because Sun-North is 180
% degrees around Z-axis off my frame used here
BennuBF_2_SN = cspice_pxform( 'IAU_BENNU', 'ORX_BENNU_SUN_NORTH', et_t0 );
theta0 = atan2(-BennuBF_2_SN(2,1),-BennuBF_2_SN(1,1));

% Locations from Science paper
% Jan 6 Near - body-fixed
pos1 = [0.0510, -0.0353, -0.2307]';

% Jan 6 Far- body-fixed
pos2 = [0.1227, -0.0359, -0.1991]';

% Old not used
pos3 = [0.0986, -0.1090, 0.1841]'; % Near radiant point, 1/19

% Jan 19 High Res 
r4 = 0.247507105834902;
latd4 = 20.63;
longd4 = 335.40;
pos4 = r4.*[cosd(longd4)*cosd(latd4); cosd(latd4)*sind(longd4); sind(latd4)]; 

% Feb 11 High Res 
r5 = 0.246413180488004;
latd5 = 20.68;
longd5 = 60.17;
pos5 = r5.*[cosd(longd5)*cosd(latd5); cosd(latd5)*sind(longd5); sind(latd5)]; 

Tmap = load('bennu_12288_v20_temperatures.txt');

ShapeModelFile = '../Shapes - 75cm delivery/Palmer/g_12580mm_spc_obj_0000n00000_v020_noComments.obj';
Polygon = GenerateVertices(ShapeModelFile);
Num_vec = [length(Polygon.Vertex), length(Polygon.FacetNormal), length(Polygon.Edges_1)];
load('Palmer_v20_facet_areas.mat','AreaF');

info = load('../Shapes - 75cm delivery/Palmer/LimbInfo.mat','vlimbs','flimbs','maxLimbR');
vlimbs = info.vlimbs;
flimbs = info.flimbs;
maxLimbR = info.maxLimbR;
rmax = 0.287955928119145;
rmin = 0.214679035301331;

% Normalizing values
norm_R = rmin; %min(maxLimbR); % [km]
norm_t = sqrt(norm_R^3/mu_ast);
norm_M = M;
norm_V = sqrt(mu_ast/norm_R);


% Solar radiation pressure
P = G1./a_helio^2;

% SRP cannonball coefficients [km/s^2]
SRP_coeffs_part = srpInfo;


N_part = length(srpInfo);

% Run simulations
t0 = 0; %40*24*3600;
tf = (437-40)*24*3600; % Bennu period is 436.65 days (wikipedia) 

options = odeset('RelTol',1e-3,'AbsTol',1e-6,'Events',@stop_conditions);

for jj = 1:N_part

    x0_temp = state_init((jj-1)*6+1:jj*6);
    x0_temp(1:3) = x0_temp(1:3)./norm_R;
    x0_temp(4:6) = x0_temp(4:6)./norm_V;
    SRP_coeffs_part_temp = SRP_coeffs_part(jj);
    A2M = SRP_coeffs_part_temp/(1 + (4/9).*.04); % Assumes used 4% reflectivity!!
    
    [T,X,te,xe,ie] = ode45(@derivs_3D, [t0 tf]./norm_t, x0_temp, options,rho,Num_vec,Polygon,length(SRP_coeffs_part_temp),P.*SRP_coeffs_part_temp./norm_R.*norm_t^2,vlimbs./norm_R,flimbs,norm_R,norm_t,w_spin,theta0,uhat_epoch,et_t0, a_helio, AreaF, Tmap, A2M);
    
    T_all{jj} = T.*norm_t;
    X(:,1:3) = X(:,1:3).*norm_R;
    X(:,4:6) = X(:,4:6).*norm_V;
    X_all{jj} = X;
    if ~isempty(te)
        Tevent_all{jj} = te.*norm_t;
        xe(:,1:3) = xe(:,1:3).*norm_R;
        xe(:,4:6) = xe(:,4:6).*norm_V;
        Xevent_all{jj} = xe;
        Index_event_all{jj} = ie;
    else
        Tevent_all{jj} = [];
        Xevent_all{jj} = [];
        Index_event_all{jj} = [];
    end
    
    disp(['Particle #' num2str(jj) ' of ' num2str(N_part) ' done!']);
    
    % save
    save(['Extend_40days_pos' num2str(caseNumber) caseStr '_redo'],'T_all','X_all','Tevent_all','Xevent_all','Index_event_all','-v7.3');
    
end




cspice_kclear;

end

function dx = derivs_3D(t,x,density,Num_vec,Polygon,N,SRP_coeffs_part,vlimbs,flimbs,norm_R,norm_t,w,theta0,uhat_epoch,et, a_epoch, AreaF, Tmap, A2M)

% Iz and Is are input as mass-normalized moments of inertia
% N is the number of particles
% SRP_coeffs_part is the 

dx = zeros(6,N);

xloop = reshape(x,6,N);

SN_2_J2000 = cspice_pxform('ORX_BENNU_SUN_NORTH','J2000',et+t*norm_t);
uhat = SN_2_J2000(:,1);
ang_sun_move = real(acos(dot(uhat_epoch,uhat)));

epochF_2_currF = [cos(ang_sun_move) -sin(ang_sun_move) 0; sin(ang_sun_move) cos(ang_sun_move) 0; 0 0 1];

sun_dir = -[cos(ang_sun_move); -sin(ang_sun_move); 0]; % nearly true (doesn't account for non-perfect retrograde pole) direction vector TOWARD sun in my working frame

Bennu_id = '2101955';
[posBennu,~]= cspice_spkpos( num2str(Bennu_id), et+t*norm_t, 'ORX_BENNU_ORBIT_FIXED', 'NONE', 'SUN' );
a_helio = posBennu(1);

% Want theta from SPICE - here need rotation of body since t0 which can
% still just do by formula below - theta0 should come from PCK
theta = w*t*norm_t + theta0;
B_ACI = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];


mu_Sun = 132712440040.944600; % Sun GM from DE424


for ii = 1:N
    
    xvec = xloop(:,ii);
    r = norm(xvec(1:3)).*norm_R;
%     rdot = dot(xvec(1:3),xvec(4:6))/r;
    
    shadowed = 0;
    inside = 0; 
    for jj = 1:size(flimbs,1)
        
        [shadowed,~] = in_triangle_C_mex(vlimbs(:,flimbs(jj,1)),vlimbs(:,flimbs(jj,2)),vlimbs(:,flimbs(jj,3)),epochF_2_currF*xvec(1:3),[-1;0;0]);
        
%         if shadowed
%             if t*norm_t > 900
%                 for kk = 1:size(flimbs,1)
%                     [inside,~] = in_triangle_C_mex(vlimbs(:,flimbs(kk,1)),vlimbs(:,flimbs(kk,2)),vlimbs(:,flimbs(kk,3)),epochF_2_currF*[0; norm(xvec(1:2)); xvec(3)],[-1;0;0]);
%                     if inside
%                         break;
%                     end
%                 end
%             end
%             break;
%         end
    end
    
    
%     if inside %|| r > 10 %r < 1 %&& rdot < 0 %Omega < -1
%         
% %         dx(6*(ii-1)+1:6*(ii-1)+6) = 0;
%         dx(:,ii) = 0;
%     
%     else
        
        dx_temp = zeros(6,1);
        dx_temp(1:3) = xvec(4:6);
        
        rvec_unN = xvec(1:3).*norm_R;
        rvec_unN_Body = B_ACI*rvec_unN;
        
        [~, acc, ~, ~, Omega] = GetShapeModelGravity_mex(density, Num_vec, rvec_unN_Body',...
        Polygon.Vertex, ...
        Polygon.FacetVertex_1, Polygon.FacetVertex_2, Polygon.FacetVertex_3,...
        Polygon.EdgesIndex, Polygon.EdgeFacetNumberDirection_1, Polygon.EdgeFacetNumberDirection_2,...
        Polygon.Edges_1, Polygon.Edges_2);
        
        dx_temp(4:6) = (B_ACI'*acc)./norm_R.*norm_t^2;
        
        
        
        if ~shadowed || xvec(1) <= 0
%             dx_temp(4) = dx_temp(4) + SRP_coeffs_part(ii);
            
            % not spinning
            dx_temp(4:6) = dx_temp(4:6) - SRP_coeffs_part(ii)*(a_epoch/a_helio)^2.*sun_dir;
            
            % spinning
%             dx_temp(4:6) = dx_temp(4:6) - SRP_coeffs_part(ii)*cos(wSpin*t*norm_t)^2*(a_epoch/a_helio)^2.*sun_dir;
        end
        
        d_pos          = -a_helio.*sun_dir;  % [km] Asteroid's position in the SCI frame
        S_vec    = d_pos + xvec(1:3).*norm_R;
        Acce_tidal = mu_Sun * (-1/norm(S_vec)^3*S_vec + 1/norm(d_pos)^3*d_pos);
        
        dx_temp(4:6) = dx_temp(4:6) + Acce_tidal./norm_R.*norm_t^2;
        
        % Thermal forces from Bennu
        mid_dir_body = -B_ACI*sun_dir;
        phi = 360 - mod(atan2d(mid_dir_body(2),mid_dir_body(1)),360);
        %a_th = thermal_accel_half(rvec_unN_Body, A2M, phi, -mid_dir_body, Num_vec(2), Polygon.FacetNormal', Polygon.FacetCenter', AreaF, Tmap, a_helio);
        a_th = thermal_accel(rvec_unN_Body, A2M, phi, -mid_dir_body, Num_vec(2), Polygon.FacetNormal', Polygon.FacetCenter', AreaF, Tmap, a_helio);
        dx_temp(4:6) = dx_temp(4:6) + (3/5).*(B_ACI'*a_th)./norm_R.*norm_t^2; % 3/5 takes out the Lambertian for re-radiation - assuming spinning and reradiating isothermally
        
        
        dx(:,ii) = dx_temp;
    
%     end

end

dx = reshape(dx,6*N,1);

end

function [value, isterminal, direction] = stop_conditions(t, y,density,Num_vec,Polygon,N,SRP_coeffs_part,vlimbs,flimbs,norm_R,norm_t,w,theta0,uhat_epoch,et, a_epoch, AreaF, Tmap, A2M)

r = norm(y(1:3));

rdot = dot(y(1:3), y(4:6));

% 4 cases to watch for
% 1: inside body when radius is less than 210 m (STOP)
% 2: outside hill radius when radius is greater than 35 km (STOP)
% 3: apoapse when rdot goes through zero decreasing
% 4: periapse when rdot goes through zero increasing
value = [r-(.210/norm_R); r-(35/norm_R); rdot; rdot];

isterminal = [1; 1; 0; 0];

direction = [0; 0; -1; 1];

end
