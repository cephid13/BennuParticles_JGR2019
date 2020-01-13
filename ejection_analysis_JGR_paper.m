% ejection analysis
caseNum = 23;
if caseNum == 20
    caseStr = '_pos2_thermal_35ths_newPos';
elseif caseNum == 21
    caseStr = '_pos1_thermal_35ths_newPos';
elseif caseNum == 22
    caseStr = '_pos5_thermal_35ths_newPos';
elseif caseNum == 23
    caseStr = '_pos4_thermal_35ths_newPos_redo';
else
    error('Don''t use any other cases in paper');    
end

ORExLineColors = [...
    0 0 0; ...          % 1 black
    210 17 71; ...      % 2 red
    255 135 65; ...     % 3 orange
    255 223 90; ...     % 4 yellow
    113 210 162; ...    % 5 light green
    64 150 181; ...     % 6 light blue
    100 85 165; ...     % 7 light purple
    12 51 114]./256;    % 8 indigo

ORExcolorMap= [ 0.3843 , 0.2941 , 0.6588; 0.3137, 0.3725 , 0.6902; 0.1961 , 0.4510 , 0.7216; 0.1451 , 0.5255 , 0.7294; 0.0471 , ...
0.6000 , 0.7333; 0.1451 ,0.6667 , 0.6980; 0.1922 ,0.7451 , 0.6549; 0.3333 , 0.7843 , 0.6431; 0.4431 , 0.8235 , 0.6353; 0.5490 , ...
0.8588 , 0.6235; 0.6471 , 0.8902 , 0.6078; 0.7373 , 0.9176 , 0.5882; 0.8196 , 0.9490 , 0.5647; 0.8745 , 0.9647 , 0.5882; 0.9294 , ...
0.9843 , 0.6118; 0.9647 , 0.9961 , 0.6588; 1.0000 , 1.0000 , 0.7098; 1.0000 , 0.9647 , 0.6510; 1.0000 , 0.9294 , 0.5765; 1.0000 , ...
0.8784 , 0.5098; 1.0000 , 0.8275 , 0.4471; 1.0000 , 0.7647 , 0.3882; 1.0000 , 0.6941 , 0.3294; 1.0000 , 0.6157 , 0.2902; 1.0000, ...
0.5294 , 0.2549; 1.0000 , 0.4510 , 0.2353; 1.0000 , 0.3647 , 0.2118; 0.9765 , 0.2941 , 0.2431; 0.9412 , 0.2235 , 0.2706; 0.8824 , ...
0.1412 , 0.2745; 0.8235 , 0.0667 , 0.2784; 0.7490 , 0.0353 , 0.2667];

%% Find impact/escape time
if caseNum <= 8
    vmags = [10:2:30 40:10:100];
else
    vmags = [10:2:30];
end
Tend_all = zeros(1606,length(vmags));
Inds_short = cell(length(vmags),1);
Inds_1_5days = cell(length(vmags),1);
Inds_5_40days = cell(length(vmags),1);
Inds_max40days = cell(length(vmags),1);
for jj = 1:length(vmags) %2:25
    load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(jj)) 'cms_lowTol' caseStr],'T','X');
    for ii = 1:1606
        ind = find(any(diff(X(:,(ii-1)*6 +1:ii*6)),2)==0,1);
        if isempty(ind)
            ind = size(X,1);
            Inds_max40days{jj} = [Inds_max40days{jj}; ii];
        end
        Tend_all(ii,jj) = T(ind)/(3600*24);
        if Tend_all(ii,jj) > 1
            if Tend_all(ii,jj) < 5
                Inds_1_5days{jj} = [Inds_1_5days{jj}; ii];
            else
                Inds_5_40days{jj} = [Inds_5_40days{jj}; ii];
            end
        end
    end
    Inds_short{jj} = [1:1606]';
    Inds_short{jj}(unique([Inds_1_5days{jj}; Inds_5_40days{jj}; Inds_max40days{jj}])) = [];
end

figure; plot(Tend_all,'o'); ylabel('Time of impact/escapse'); 

figure; hold on; histogram(Tend_all,50); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log')

figure; htest = histogram(reshape(Tend_all,1606*length(vmags),1),[0 .25 .5 .75 1 5 40],'Normalization','probability');

figure; hold on; bar(htest.Values.*100); ylabel('\% of Objects');
set(gca,'XTick',[1:6]);
set(gca,'XTickLabels',{'0 - 6 hrs','6 - 12 hrs','12 - 18 hrs','18-24 hrs','1 - 5 days','5 - 40+ days'})
set(gca,'XTickLabelRotation',75)


%% Get launch conditions and A/m ratios for sets of interest
part_radius_short = [];
part_radius_1_5days = [];
part_radius_5_40days = [];
part_radius_max40days = [];
Am_short = [];
Am_1_5days = [];
Am_5_40days = [];
Am_max40days = [];
velmag_short = [];
velmag_1_5days = [];
velmag_5_40days = [];
velmag_max40days = [];
velaz_short = [];
velaz_1_5days = [];
velaz_5_40days = [];
velaz_max40days = [];
velel_short = [];
velel_1_5days = [];
velel_5_40days = [];
velel_max40days = [];
times_short = [];
times_1_5days = [];
times_5_40days = [];
times_max40days = [];

% recreate radius set here since didn't save it
num_el = 7; % including 90 degrees
num_az = 12; 
num_rad = 22;
N_part = ((num_el-1)*num_az + 1)*num_rad;
radius_part = zeros(N_part,1);
radius_part(1:num_rad) = [.1 .5 1:20];
radius_part(num_rad+1:end) = repmat(radius_part(1:num_rad),num_az*(num_el-1),1);

% load az, el, SRP from any file
load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(1)) 'cms_lowTol' caseStr],'SRP_coeffs_part','vel_az','vel_el');

% Gather data for each set
% <1 day
for ii = 1:length(vmags)
    Inds = Inds_short{ii};
    part_radius_short = [part_radius_short; radius_part(Inds)]; % cm
    Am_short = [Am_short; SRP_coeffs_part(Inds)'./(1+.04*4/9).*1e6]; % m^2/kg 
    velmag_short = [velmag_short; vmags(ii).*ones(length(Inds),1)]; % cm/s
    velaz_short = [velaz_short; vel_az(Inds)]; % rad, measured from East
    velel_short = [velel_short; vel_el(Inds)]; % rad, measured from tangent plane   
    times_short = [times_short; Tend_all(Inds,ii)];
end

% 1-5 days
for ii = 1:length(vmags)
    Inds = Inds_1_5days{ii};
    part_radius_1_5days = [part_radius_1_5days; radius_part(Inds)]; % cm
    Am_1_5days = [Am_1_5days; SRP_coeffs_part(Inds)'./(1+.04*4/9).*1e6]; % m^2/kg 
    velmag_1_5days = [velmag_1_5days; vmags(ii).*ones(length(Inds),1)]; % cm/s
    velaz_1_5days = [velaz_1_5days; vel_az(Inds)]; % rad, measured from East
    velel_1_5days = [velel_1_5days; vel_el(Inds)]; % rad, measured from tangent plane   
    times_1_5days = [times_1_5days; Tend_all(Inds,ii)];
end

% 5-40 days
for ii = 1:length(vmags)
    Inds = Inds_5_40days{ii};
    part_radius_5_40days = [part_radius_5_40days; radius_part(Inds)]; % cm
    Am_5_40days = [Am_5_40days; SRP_coeffs_part(Inds)'./(1+.04*4/9).*1e6]; % m^2/kg 
    velmag_5_40days = [velmag_5_40days; vmags(ii).*ones(length(Inds),1)]; % cm/s
    velaz_5_40days = [velaz_5_40days; vel_az(Inds)]; % rad, measured from East
    velel_5_40days = [velel_5_40days; vel_el(Inds)]; % rad, measured from tangent plane
    times_5_40days = [times_5_40days; Tend_all(Inds,ii)];
end

% 40 days only (max length)
for ii = 1:length(vmags)
    Inds = Inds_max40days{ii};
    part_radius_max40days = [part_radius_max40days; radius_part(Inds)]; % cm
    Am_max40days = [Am_max40days; SRP_coeffs_part(Inds)'./(1+.04*4/9).*1e6]; % m^2/kg 
    velmag_max40days = [velmag_max40days; vmags(ii).*ones(length(Inds),1)]; % cm/s
    velaz_max40days = [velaz_max40days; vel_az(Inds)]; % rad, measured from East
    velel_max40days = [velel_max40days; vel_el(Inds)]; % rad, measured from tangent plane    
    times_max40days = [times_max40days; Tend_all(Inds,ii)];
end

% Plot data
% Vel Az vs El vs Mag
figure; hold on; plot3(velaz_1_5days, velel_1_5days, velmag_1_5days, 'Marker','o','MarkerSize',15,'MarkerFaceColor',[0 .45 .74],'LineStyle','none'); plot3(velaz_5_40days, velel_5_40days, velmag_5_40days, 'Marker','^','MarkerSize',15,'MarkerFaceColor',[.85 .32 .01],'LineWidth',0.5,'LineStyle','none'); plot3(velaz_max40days, velel_max40days, velmag_max40days, 'Marker','+','MarkerSize',15,'Color',[0 0 0],'LineStyle','none'); 
xlabel('Launch Az [rad]'); ylabel('Launch El [rad]'); zlabel('Launch Vel [cm/s]');
legend('1-5 days','5-40 days','40 days')

% A/M ratio vs vel Mag
figure; hold on; plot(Am_1_5days, velmag_1_5days, 'Marker','o','MarkerSize',15,'MarkerFaceColor',[0 .45 .74],'LineStyle','none'); plot(Am_5_40days, velmag_5_40days, 'Marker','^','MarkerSize',15,'MarkerFaceColor',[.85 .32 .01],'LineWidth',0.5,'LineStyle','none'); plot(Am_max40days, velmag_max40days, 'Marker','+','MarkerSize',15,'Color',[0 0 0],'LineStyle','none'); 
set(gca,'XScale','log')
xlabel('A/m [m$^2$/kg]'); ylabel('Launch Vel [cm/s]');
legend('1-5 days','5-40 days','40 days')

%% Classify based on periapse passages for all cases
Inds_direct_escape = cell(length(vmags),1);
Inds_suborbital = cell(length(vmags),1);
Inds_suborbital_bad = cell(length(vmags),1);
Inds_escape = cell(length(vmags),1);
Inds_orbital = cell(length(vmags),1);
Inds_multi_orbital = cell(length(vmags),1);
Inds_direct_escape_long = cell(length(vmags),1);
Inds_escape_long = cell(length(vmags),1);
Tend_all_direct_escape = [];
Tend_all_suborbital = [];
Tend_all_escape = [];
Tend_all_orbital = [];
rpop_low = zeros(length(T),1);
rpop_med = zeros(length(T),1);
rpop_high = zeros(length(T),1);
rpop_surf = zeros(length(T),1);
for jj = 1:length(vmags) 
    
    load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(jj)) 'cms_lowTol' caseStr],'T','X');
    for ii = 1:1606
        rmag_all{ii} = sqrt(dot_fast_sparse(X(:,(ii-1)*6 + 1:(ii-1)*6 + 3)', X(:,(ii-1)*6 + 1:(ii-1)*6 + 3)'));
        
        rpop_ind_surf = find(T./3600/24 >= Tend_all(ii,jj));
        rpop_tend_ind = min(rpop_ind_surf);
        if isempty(rpop_tend_ind)
            rpop_tend_ind = length(T);
        else
            rpop_tend_ind = rpop_tend_ind - 1;
        end
        rpop_ind_low = find(rmag_all{ii}(1:rpop_tend_ind) <= 1);
        rpop_ind_med = find(rmag_all{ii}(1:rpop_tend_ind) <= 5 & rmag_all{ii}(1:rpop_tend_ind) > 1);
        rpop_ind_high = find(rmag_all{ii}(1:rpop_tend_ind) > 5);
        rpop_low(rpop_ind_low) = rpop_low(rpop_ind_low) + 1;
        rpop_med(rpop_ind_med) = rpop_med(rpop_ind_med) + 1;
        rpop_high(rpop_ind_high) = rpop_high(rpop_ind_high) + 1;
        rpop_surf(rpop_ind_surf) = rpop_surf(rpop_ind_surf) + 1;
        
        drmag = diff(rmag_all{ii});
        sdrmag = sign(drmag);
        dsdrmag = diff(sdrmag);
        temp_indp = find(dsdrmag==2)+1;
        if ~isempty(temp_indp) && length(temp_indp) == 1 && dsdrmag(temp_indp) == -1
            disp(['Fake periapse in V = ' num2str(jj) ' case ' num2str(ii)]);
            indp_all{ii} = [];
        else
            indp_all{ii} = temp_indp;
        end
        inda_all{ii} = find(dsdrmag==-2)+1;
        
        T_periapse{ii} = T(indp_all{ii});
        T_apoapse{ii} = T(inda_all{ii});
        Periods_periapse{ii} = diff(T_periapse{ii});
        Periods_apoapse{ii} = diff(T_apoapse{ii});
        
        if isempty(indp_all{ii}) && isempty(inda_all{ii})
            if rmag_all{ii}(end) > 1
                % Direct escape
                Inds_direct_escape{jj} = [Inds_direct_escape{jj}; ii];
                if Tend_all(Inds_direct_escape{jj}(end),jj) > 10
%                     figure(fig_long_escape); plot(T./3600./24,rmag_all{ii})
                    Inds_direct_escape_long{jj} = [Inds_direct_escape_long{jj}; ii Tend_all(Inds_direct_escape{jj}(end),jj)];
                end
            else
                % suborbital
                Inds_suborbital{jj} = [Inds_suborbital{jj}; ii];
                Inds_suborbital_bad{jj} = [Inds_suborbital_bad{jj}; ii];
            end
        elseif rmag_all{ii}(end) >= 35
            % Later escape
            Inds_escape{jj} = [Inds_escape{jj}; ii];
            if Tend_all(Inds_escape{jj}(end),jj) > 10
%                 figure(fig_long_escape); plot(T./3600./24,rmag_all{ii})
                Inds_escape_long{jj} = [Inds_escape_long{jj}; ii Tend_all(Inds_escape{jj}(end),jj)];
            end
        else
            if isempty(indp_all{ii})
                % suborbital
                Inds_suborbital{jj} = [Inds_suborbital{jj}; ii];
            else
                if indp_all{ii}(1) > find(T./3600/24 == Tend_all(ii,jj),1)
                    % bad identification of periapse, should be suborbital
                    Inds_suborbital{jj} = [Inds_suborbital{jj}; ii];
                    indp_all{ii} = [];
                    T_periapse{ii} = [];
                else
                    % Orbital
                    Inds_orbital{jj} = [Inds_orbital{jj}; ii];
                    if length(indp_all{ii}) > 1
                        Inds_multi_orbital{jj} = [Inds_multi_orbital{jj}; ii];
                    end
                end
            end
        end
        
    end
    rmag_really_all{jj} = rmag_all;
    indp_really_all{jj} = indp_all;
    inda_really_all{jj} = inda_all;
    T_periapse_all{jj} = T_periapse;
    T_apoapse_all{jj} = T_apoapse;
    Periods_periapse_all{jj} = Periods_periapse;
    Periods_apoapse_all{jj} = Periods_apoapse;
    
    Tend_all_direct_escape = [Tend_all_direct_escape; Tend_all(Inds_direct_escape{jj},jj)];
    Tend_all_suborbital = [Tend_all_suborbital; Tend_all(Inds_suborbital{jj},jj)];
    Tend_all_escape = [Tend_all_escape; Tend_all(Inds_escape{jj},jj)];
    Tend_all_orbital = [Tend_all_orbital; Tend_all(Inds_orbital{jj},jj)];
end

figure; hold on; 
plot(T./3600./24,rpop_low./17666.*100,'Color',ORExLineColors(5,:));
plot(T./3600./24,rpop_med./17666.*100,'Color',ORExLineColors(3,:));
plot(T./3600./24,rpop_high./17666.*100,'Color',ORExLineColors(2,:));
plot(T./3600./24,rpop_surf./17666.*100,'Color',ORExLineColors(7,:));
xlabel('Time [days]'); ylabel('$\%$ particles'); legend('$<$ 1 km','1-5 km','5+ km','crashed/escaped','Location','best');
xlim([0 10]);

figure; hold on; fig_long_escape = gcf; 
for jj = 1:length(vmags)
    if ~isempty(Inds_direct_escape_long{jj})
        for ii = 1:size(Inds_direct_escape_long{jj},1)
            plot(T./3600./24,rmag_really_all{jj}{Inds_direct_escape_long{jj}(ii,1)},'r')
        end
    end
    if ~isempty(Inds_escape_long{jj})
        for ii = 1:size(Inds_escape_long{jj},1)
            plot(T./3600./24,rmag_really_all{jj}{Inds_escape_long{jj}(ii,1)},'b')
        end
    end
end
ylabel('r [km]'); xlabel('t [days]'); title('Red = direct escape, Blue = escape post periapse');
% In Pos1, velCase=5, numbers 834 and 1054 are the direct escape cases that take > 40 days
% Also in Pos1, long escape plotted case #481 (which is not direct) is interesting

figure; hold on; htest1 = histogram(Tend_all_direct_escape,0:41); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Direct Escape');
figure; hold on; htest2 = histogram(Tend_all_suborbital,0:41); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Suborbital');
figure; hold on; htest3 = histogram(Tend_all_escape,0:41); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Escape');
figure; hold on; htest4 = histogram(Tend_all_orbital,0:41); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Orbital');

Tend_all_stacked = [htest1.Values; htest2.Values; htest4.Values; htest3.Values]';
figure; bar(.5:1:40.5,Tend_all_stacked,'stacked'); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); 
legend('Direct Escape','Suborbital','Orbital','Escape')

% figure; htest = histogram(reshape(Tend_all,1606*length(vmags),1),[0 .25 .5 .75 1 5 40],'Normalization','probability');
% 
% figure; hold on; bar(htest.Values.*100); ylabel('\% of Objects');
% set(gca,'XTick',[1:6]);
% set(gca,'XTickLabels',{'0 - 6 hrs','6 - 12 hrs','12 - 18 hrs','18-24 hrs','1 - 5 days','5 - 40+ days'})
% set(gca,'XTickLabelRotation',75)

% Revise to add extended info for long lived bodies
% Need to load extended data
% figure out the final event index code to see their fate (1 crash, 2 escape)
% can get their final time from the last event time, plus 40 days from
% initial simulation
count_to_v_case = zeros(0,2);
for ii = 1:length(vmags)
    count_to_v_case = [count_to_v_case; ii.*ones(length(Inds_max40days{ii}),1) Inds_max40days{ii}];
end
if caseNum == 21
    extend = load('Extend_40days_pos1_pos1_thermal_35ths_newPos_redo.mat');
elseif caseNum == 20
    extend = load('Extend_40days_pos2_pos2_thermal_35ths_newPos_redo.mat');
elseif caseNum == 22
    extend = load('Extend_40days_pos5_pos5_thermal_35ths_newPos_redo.mat');
elseif caseNum == 23
    extend = load('Extend_40days_pos4_pos4_thermal_35ths_newPos_redo_redo.mat');
else
    error('fix loading extended data');
end

for ii = 1:size(count_to_v_case,1)
    if ~isempty(find(Inds_direct_escape{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1))
        fateNum = 1;
        modIndex = find(Inds_direct_escape{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1);
        timeIndex = modIndex;
        for jj = 1:count_to_v_case(ii,1)-1
            timeIndex = timeIndex + length(Inds_direct_escape{jj});
        end
        if length(extend.Index_event_all{ii})>1
            disp('Shouldn''t be a direct escape any more! But need to fix if you want to update that');
        end
        Tend_all_direct_escape(timeIndex) = Tend_all_direct_escape(timeIndex) + extend.Tevent_all{ii}(end)/3600/24;
        
    elseif ~isempty(find(Inds_suborbital{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1))
        fateNum = 2;
        disp(['I didn''t expect any in suborbital... weird long lofting! Case #' num2str(ii)]);
        
        modIndex = find(Inds_suborbital{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1);
        timeIndex = modIndex;
        for jj = 1:count_to_v_case(ii,1)-1
            timeIndex = timeIndex + length(Inds_suborbital{jj});
        end
        if length(extend.Index_event_all{ii}) > 1
            if extend.Index_event_all{ii}(end) == 1
                disp(['Case #' num2str(ii) 'Changed to an orbital case!']);
                modIndexUpdate = find(Inds_orbital{count_to_v_case(ii,1)}>count_to_v_case(ii,2),1);
                if isempty(modIndexUpdate)
                    modIndexUpdate = 1;
                end
                timeIndexUpdate = modIndexUpdate;
                for jj = 1:count_to_v_case(ii,1)-1
                    timeIndexUpdate = timeIndexUpdate + length(Inds_orbital{jj});
                end
                Tend_all_orbital = [Tend_all_orbital(1:timeIndexUpdate-1); Tend_all_suborbital(timeIndex) + extend.Tevent_all{ii}(end)/3600/24; Tend_all_orbital(timeIndexUpdate:end)];
                Inds_orbital{count_to_v_case(ii,1)} = sort([Inds_orbital{count_to_v_case(ii,1)}; count_to_v_case(ii,2)]);
                Tend_all_suborbital(timeIndex) = [];
                Inds_suborbital{count_to_v_case(ii,1)}(modIndex) = [];
            elseif extend.Index_event_all{ii}(end) == 2
                disp(['Case #' num2str(ii) 'Changed to an escape case!']);
                modIndexUpdate = find(Inds_escape{count_to_v_case(ii,1)}>count_to_v_case(ii,2),1);
                if isempty(modIndexUpdate)
                    modIndexUpdate = 1;
                end
                timeIndexUpdate = modIndexUpdate;
                for jj = 1:count_to_v_case(ii,1)-1
                    timeIndexUpdate = timeIndexUpdate + length(Inds_escape{jj});
                end
                Tend_all_escape = [Tend_all_escape(1:timeIndexUpdate-1); Tend_all_suborbital(timeIndex) + extend.Tevent_all{ii}(end)/3600/24; Tend_all_escape(timeIndexUpdate:end)];
                Inds_escape{count_to_v_case(ii,1)} = sort([Inds_escape{count_to_v_case(ii,1)}; count_to_v_case(ii,2)]);
                Tend_all_suborbital(timeIndex) = [];
                Inds_suborbital{count_to_v_case(ii,1)}(modIndex) = [];
            end
        elseif extend.Index_event_all{ii}(end) == 2
            disp(['Case #' num2str(ii) 'Changed to a direct escape case!']);
            modIndexUpdate = find(Inds_direct_escape{count_to_v_case(ii,1)}>count_to_v_case(ii,2),1);
            if isempty(modIndexUpdate)
                modIndexUpdate = 1;
            end
            timeIndexUpdate = modIndexUpdate;
            for jj = 1:count_to_v_case(ii,1)-1
                timeIndexUpdate = timeIndexUpdate + length(Inds_direct_escape{jj});
            end
            Tend_all_direct_escape = [Tend_all_direct_escape(1:timeIndexUpdate-1); Tend_all_suborbital(timeIndex) + extend.Tevent_all{ii}(end)/3600/24; Tend_all_direct_escape(timeIndexUpdate:end)];
            Inds_direct_escape{count_to_v_case(ii,1)} = sort([Inds_direct_escape{count_to_v_case(ii,1)}; count_to_v_case(ii,2)]);
            Tend_all_suborbital(timeIndex) = [];
            Inds_suborbital{count_to_v_case(ii,1)}(modIndex) = [];
        else
            % just a long trajectory suborbital case
            Tend_all_suborbital(timeIndex) = Tend_all_suborbital(timeIndex) + extend.Tevent_all{ii}(end)/3600/24;
        end
        
    elseif ~isempty(find(Inds_escape{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1))
        fateNum = 3;
        
        modIndex = find(Inds_escape{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1);
        timeIndex = modIndex;
        for jj = 1:count_to_v_case(ii,1)-1
            timeIndex = timeIndex + length(Inds_escape{jj});
        end
        if extend.Index_event_all{ii}(end) == 1
            disp(['Case #' num2str(ii) ' Needs to be changed to an orbital case - not escaping']);
        end
        Tend_all_escape(timeIndex) = Tend_all_escape(timeIndex) + extend.Tevent_all{ii}(end)/3600/24;
        
        
    elseif ~isempty(find(Inds_orbital{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1))
        fateNum = 4;
        
        modIndex = find(Inds_orbital{count_to_v_case(ii,1)}==count_to_v_case(ii,2),1);
        timeIndex = modIndex;
        for jj = 1:count_to_v_case(ii,1)-1
            timeIndex = timeIndex + length(Inds_orbital{jj});
        end
        if extend.Index_event_all{ii}(end) > 2
            disp(['Case #' num2str(ii) ' still orbiting!']);
        end
        if extend.Index_event_all{ii}(end) == 2
            disp(['Case #' num2str(ii) ' changed to an escape case']);
            modIndexUpdate = find(Inds_escape{count_to_v_case(ii,1)}>count_to_v_case(ii,2),1);
            if isempty(modIndexUpdate)
                modIndexUpdate = 1;
            end
            timeIndexUpdate = modIndexUpdate;
            for jj = 1:count_to_v_case(ii,1)-1
                timeIndexUpdate = timeIndexUpdate + length(Inds_escape{jj});
            end
            Tend_all_escape = [Tend_all_escape(1:timeIndexUpdate-1); Tend_all_orbital(timeIndex) + extend.Tevent_all{ii}(end)/3600/24; Tend_all_escape(timeIndexUpdate:end)];
            Inds_escape{count_to_v_case(ii,1)} = sort([Inds_escape{count_to_v_case(ii,1)}; count_to_v_case(ii,2)]);
            Tend_all_orbital(timeIndex) = [];
            Inds_orbital{count_to_v_case(ii,1)}(modIndex) = [];
        elseif extend.Index_event_all{ii}(end) == 1 % still classified as orbital, but it crashed
            Tend_all_orbital(timeIndex) = Tend_all_orbital(timeIndex) + extend.Tevent_all{ii}(end)/3600/24;
        else
            % still orbiting after a year!
            Tend_all_orbital(timeIndex) = Tend_all_orbital(timeIndex) + extend.T_all{ii}(end)/3600/24;
            if Tend_all_orbital(timeIndex) < 437
                Tend_all_orbital(timeIndex) = 437.1;
            end
        end
    end
end

figure; hold on; htest5 = histogram(Tend_all_direct_escape,0:438); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Direct Escape');
figure; hold on; htest6 = histogram(Tend_all_suborbital,0:438); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Suborbital');
figure; hold on; htest7 = histogram(Tend_all_escape,0:438); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Escape');
figure; hold on; htest8 = histogram(Tend_all_orbital,0:438); ylabel('\# of Objects'); xlabel('Time of impact/escape'); set(gca,'YScale','log'); title('Orbital');

Tend_all_stacked = [htest8.Values; htest5.Values; htest6.Values; htest7.Values]';
figure; bar(.5:1:437.5,Tend_all_stacked,1,'stacked'); ylabel('\# of Objects'); xlabel('Time of impact/escape'); %set(gca,'YScale','log'); 
legend('Orbital','Direct Escape','Suborbital','Escape')
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(4,:); % Escape
barAxis.Children(2).FaceColor = ORExLineColors(7,:); % Suborbital
barAxis.Children(3).FaceColor = ORExLineColors(2,:); % Direct escape
barAxis.Children(4).FaceColor = ORExLineColors(5,:); % Orbital
% NEED TO PLAY A COUPLE PLOTTING GAMES - MAKE THE FIRST DAY BAR WIDER, AND
% GET RID OF ALL THE OPEN SPACE BETWEEN ~150 - 425 DAYS
% maybe use 'Normalization','cumcount' to do this with the histograms..

Tend_all_stacked_combined = [Tend_all_stacked(1:50,:); sum(Tend_all_stacked(51:100,:)); sum(Tend_all_stacked(101:150,:)); sum(Tend_all_stacked(151:200,:)); sum(Tend_all_stacked(201:250,:)); sum(Tend_all_stacked(251:300,:));  sum(Tend_all_stacked(301:350,:)); sum(Tend_all_stacked(351:400,:));  sum(Tend_all_stacked(401:437,:)); Tend_all_stacked(438,:)];
figure; bar([.5:1:58.5],Tend_all_stacked_combined,1,'stacked'); ylabel('\# of Objects (out of 17666)'); xlabel('Time of impact/escape [days]'); %set(gca,'YScale','log'); 
legend('Orbital','Direct Escape','Suborbital','Escape')
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(4,:); % Escape
barAxis.Children(2).FaceColor = ORExLineColors(7,:); % Suborbital
barAxis.Children(3).FaceColor = ORExLineColors(2,:); % Direct escape
barAxis.Children(4).FaceColor = ORExLineColors(5,:); % Orbital
xlim([0 50])

%% update final plot for next figures
% zoom settings
ylim([0 100]);
legend('off');
xlabel('')
ylabel('')

% long results settings
xlim([51 59])
set(gca,'XTick',[51.5:1:58.5])
set(gca,'XTickLabels',{'50-100','100-150','150-200','200-250','250-300','350-400','400-437','Surviving'})
legend('Orbital','Direct Escape','Suborbital','Escape');
ylabel('\# of Objects (out of 17666)'); xlabel('Time of impact/escape [days]');
ylim('auto')

%% Mapping
vmags = [10:2:30];

omega_vec = [0;0;(2*pi/(4.297461*3600))];

w = (2*pi/(4.297461*3600));
if caseNum == 19 || caseNum == 23
    theta0 = -1.503584070566941;
elseif caseNum == 20 || caseNum == 21
    theta0 = -1.468496493225085;
elseif caseNum == 22
    theta0 = -2.602938699489971;
else
    error('Check case to set theta0');
end


Inds_reimpact = cell(length(vmags),1);
Lat0_all = zeros(1,length(vmags)); 
Long0_all = zeros(1,length(vmags)); 
Latf_all = zeros(1606,length(vmags)); 
Longf_all = zeros(1606,length(vmags)); 
tot_impact = 0;
rf_all = zeros(3,1606*length(vmags));
rf_ind = 1;
% vertical_speeds = zeros(1606,length(vmags)); 
% impact_angle = zeros(1606,length(vmags)); 
for jj = 1:length(vmags) 
    load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(jj)) 'cms_lowTol' caseStr],'T','X');
    B_ACI = [cos(theta0) sin(theta0) 0; -sin(theta0) cos(theta0) 0; 0 0 1];
    r0 = B_ACI*X(1,1:3)';
    lat0 = asin(r0(3)/norm(r0));
    long0 = atan2(r0(2),r0(1));
    if long0 < 0
        long0 = long0 + 2*pi;
    end
    longEastWest = mod(long0 + pi, 2*pi);
    
    Lat0_all(jj) = lat0;
    Long0_all(jj) = long0;
    
    for ii = 1:1606
        
        % if it survives 40 days, don't bother continuing analysis
        if ~isempty(find(Inds_max40days{jj} == ii,1))
             Latf_all(ii,jj) = NaN;
             Longf_all(ii,jj) = NaN;
%              vertical_speeds(ii,jj) = NaN;
%              impact_angle(ii,jj) = NaN;
            continue
        end
        
        theta = w*T(end) + theta0;
        B_ACI = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
        
        rf = B_ACI*real(X(end,(ii-1)*6 +1:(ii-1)*6 +3)');
        
        % if it escapes, ignore
        if norm(rf) > 1
             Latf_all(ii,jj) = NaN;
             Longf_all(ii,jj) = NaN;
%              vertical_speeds(ii,jj) = NaN;
%              impact_angle(ii,jj) = NaN;
            continue
        else
            Inds_reimpact{jj} = [Inds_reimpact{jj}; ii];
        end
        
        
        rf_all(:,rf_ind) = rf;
        rf_ind = rf_ind + 1;
        
        latf = asin(rf(3)/norm(rf));
        longf = atan2(rf(2),rf(1));
        if longf < 0
            longf = longf + 2*pi;
        end
        
        Latf_all(ii,jj) = latf;
        Longf_all(ii,jj) = longf;
        
%         % landing speeds
%         vi = real(X(end,(ii-1)*6 +4:(ii-1)*6 +6)');
%         ri = real(X(end,(ii-1)*6 +1:(ii-1)*6 +3)');
%         
%         vb = vi - cross(omega_vec, ri);
%         
%         up = ri./norm(ri);
%         east = cross(omega_vec, ri)./norm(cross(omega_vec, ri));
%         north = cross(up,east)./norm(cross(up,east));
%         
%         vb_enu = [east north up]'*vb;
% %         body_vels{jj} = [body_vels{jj} vb_enu];
%         
% %         impact_speeds(counter) = norm(vb);
%         vertical_speeds(ii,jj) = vb_enu(3);
%         impact_angle(ii,jj) = acosd(dot(vb,ri)/norm(vb)/norm(ri));
        
    end
    
    tot_impact = tot_impact + length(Inds_reimpact{jj});
    
end

if rf_ind == tot_impact+1
    rf_all(:,rf_ind:end) = [];
else
    error('Bad rf_ind indexing?')
end

% all data
figure; hold on; 
plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','s','Color',ORExLineColors(2,:),'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',ORExLineColors(2,:),'MarkerSize',18);
scatter(180/pi.*reshape(Longf_all(~isnan(Longf_all)),tot_impact,1),180/pi.*reshape(Latf_all(~isnan(Latf_all)),tot_impact,1),'filled','MarkerFaceColor',ORExLineColors(5,:)); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]);
set(gca,'XTick',0:90:360,'YTick',-90:45:90);
legend('Launch Site','Reimpacts','Orientation','horizontal','Location','best')
grid on;

% figure; plot(180/pi.*reshape(Latf_all(~isnan(Latf_all)),tot_impact,1), reshape(impact_angle(~isnan(impact_angle)),tot_impact,1),'x'); xlabel('Latitude [deg]'); ylabel('Impact Angle (180 = radial in) [deg]');
% figure; plot(reshape(impact_angle(~isnan(impact_angle)),tot_impact,1), reshape(vertical_speeds(~isnan(vertical_speeds)),tot_impact,1),'x'); xlabel('Impact Angle (180 = radial in) [deg]'); ylabel('Vertical speed [km/s]');

figure; hist_lat = histogram(180/pi.*reshape(Latf_all(~isnan(Latf_all)),tot_impact,1), -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);
figure; hist_long = histogram(180/pi.*reshape(Longf_all(~isnan(Longf_all)),tot_impact,1), 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);

figure; barh(-85:10:85,hist_lat.Values,1,'stacked'); xlabel(['\# of Objects (out of ' num2str(tot_impact) ')']); ylabel('Latitude [deg]'); ylim([-90 90]);
set(gca,'YTick',-90:45:90);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(5,:); % Orbital
xvals = xlim; 
hold on; line([0 .99*xvals(2)],[180/pi*Lat0_all(1) 180/pi*Lat0_all(1)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('Reimpacts','Launch Site')

figure; bar(5:10:355,hist_long.Values,1,'stacked'); ylabel(['\# of Objects (out of ' num2str(tot_impact) ')']); xlabel('East Longitude [deg]'); xlim([0 360]);
set(gca,'XTick',0:90:360);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(5,:); % Orbital
yvals = ylim; 
hold on; line([180/pi*Long0_all(1) 180/pi*Long0_all(1)],[0 .99*yvals(2)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('Reimpacts','Launch Site')

figure;
hold on; 
hall = histogram2(180/pi.*reshape(Longf_all(~isnan(Longf_all)),tot_impact,1),180/pi.*reshape(Latf_all(~isnan(Latf_all)),tot_impact,1),0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);



% suborbital vs other
Longf_sub_all = [];
Latf_sub_all = [];
Longf_orb_all = [];
Latf_orb_all = [];
for jj = 1:length(vmags)
    sub = Longf_all(Inds_suborbital{jj},jj);
    Longf_sub_all = [Longf_sub_all; sub(~isnan(sub))];

    sub = Latf_all(Inds_suborbital{jj},jj);
    Latf_sub_all = [Latf_sub_all; sub(~isnan(sub))];    

    
    sub = Longf_all(Inds_orbital{jj},jj);
    Longf_orb_all = [Longf_orb_all; sub(~isnan(sub))];

    sub = Latf_all(Inds_orbital{jj},jj);
    Latf_orb_all = [Latf_orb_all; sub(~isnan(sub))];
    
end

figure; hold on; 
plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','s','Color',ORExLineColors(2,:),'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',ORExLineColors(2,:),'MarkerSize',18);
scatter(180/pi.*Longf_orb_all,180/pi.*Latf_orb_all,'filled','MarkerFaceColor',ORExLineColors(5,:));
scatter(180/pi.*Longf_sub_all,180/pi.*Latf_sub_all,'filled','MarkerFaceColor',ORExLineColors(7,:)); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]);
plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','o','MarkerFaceColor',ORExLineColors(2,:));
legend('Launch Site','Orbital','Suborbital','Orientation','horizontal','Location','best')
set(gca,'XTick',0:90:360,'YTick',-90:45:90); 
grid on;

figure; hist_sub_long = histogram(180/pi.*Longf_sub_all, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_sub_lat = histogram(180/pi.*Latf_sub_all, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);
figure; hist_orb_long = histogram(180/pi.*Longf_orb_all, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_orb_lat = histogram(180/pi.*Latf_orb_all, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);

Latf_sub_orb_stacked = [hist_sub_lat.Values; hist_orb_lat.Values]';
figure; barh(-85:10:85,Latf_sub_orb_stacked,1,'stacked'); xlabel(['\# of Objects (out of ' num2str(tot_impact) ')']); ylabel('Latitude [deg]'); ylim([-90 90]);
set(gca,'YTick',-90:45:90);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(5,:); % Orbital
barAxis.Children(2).FaceColor = ORExLineColors(7,:); % Suborbital
xvals = xlim; 
hold on; line([0 .99*xvals(2)],[180/pi*Lat0_all(1) 180/pi*Lat0_all(1)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('Suborbital','Orbital','Launch Site')

Longf_sub_orb_stacked = [hist_sub_long.Values; hist_orb_long.Values]';
figure; bar(5:10:355,Longf_sub_orb_stacked,1,'stacked'); ylabel(['\# of Objects (out of ' num2str(tot_impact) ')']); xlabel('East Longitude [deg]'); xlim([0 360]);
set(gca,'XTick',0:90:360);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(5,:); % Orbital
barAxis.Children(2).FaceColor = ORExLineColors(7,:); % Suborbital
yvals = ylim; 
hold on; line([180/pi*Long0_all(1) 180/pi*Long0_all(1)],[0 .99*yvals(2)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('Suborbital','Orbital','Launch Site')

figure;
subplot(2,1,1); hsub = histogram2(Longf_sub_all.*180/pi,Latf_sub_all.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); title('Suborbital'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
subplot(2,1,2); horb = histogram2(Longf_orb_all.*180/pi,Latf_orb_all.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); title('Orbital'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
cmax = max([max(max(hsub.Values)) max(max(horb.Values))]);
subplot(2,1,1); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);
subplot(2,1,2); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);

% elevation bands
Inds_el_high = find(round(vel_el*180/pi)>=60);
Inds_el_med = find(round(vel_el*180/pi)>=30 & round(vel_el*180/pi)<60);
Inds_el_low = find(round(vel_el*180/pi)<30);
Longf_el_high = [];
Latf_el_high = [];
Longf_el_med = [];
Latf_el_med = [];
Longf_el_low = [];
Latf_el_low = [];
for jj = 1:length(vmags)
    sub = Longf_all(Inds_el_high,jj);
    Longf_el_high = [Longf_el_high; sub(~isnan(sub))];

    sub = Latf_all(Inds_el_high,jj);
    Latf_el_high = [Latf_el_high; sub(~isnan(sub))];    
    
    sub = Longf_all(Inds_el_med,jj);
    Longf_el_med = [Longf_el_med; sub(~isnan(sub))];

    sub = Latf_all(Inds_el_med,jj);
    Latf_el_med = [Latf_el_med; sub(~isnan(sub))]; 
    
    sub = Longf_all(Inds_el_low,jj);
    Longf_el_low = [Longf_el_low; sub(~isnan(sub))];

    sub = Latf_all(Inds_el_low,jj);
    Latf_el_low = [Latf_el_low; sub(~isnan(sub))]; 
    
end

figure; hold on; 
plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','s','Color',ORExLineColors(2,:),'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',ORExLineColors(2,:),'MarkerSize',18);
scatter(180/pi.*Longf_el_high,180/pi.*Latf_el_high,'filled','MarkerFaceColor',ORExLineColors(8,:));
scatter(180/pi.*Longf_el_med,180/pi.*Latf_el_med,'filled','MarkerFaceColor',ORExLineColors(3,:));
scatter(180/pi.*Longf_el_low,180/pi.*Latf_el_low,'filled','MarkerFaceColor',ORExLineColors(5,:)); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]);
legend('Launch Site','$\theta \geq 60^\circ$', '$30^\circ \leq \theta < 60^\circ$','$0^\circ \leq \theta < 30^\circ$','Orientation','horizontal','Location','best')
set(gca,'XTick',0:90:360,'YTick',-90:45:90); 
grid on;

figure; hist_long_el_high = histogram(180/pi.*Longf_el_high, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_lat_el_high = histogram(180/pi.*Latf_el_high, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);
figure; hist_long_el_med = histogram(180/pi.*Longf_el_med, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_lat_el_med = histogram(180/pi.*Latf_el_med, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);
figure; hist_long_el_low = histogram(180/pi.*Longf_el_low, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_lat_el_low = histogram(180/pi.*Latf_el_low, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);

Latf_el_stacked = [hist_lat_el_high.Values; hist_lat_el_med.Values; hist_lat_el_low.Values]';
figure; barh(-85:10:85,Latf_el_stacked,1,'stacked'); xlabel(['\# of Objects (out of ' num2str(tot_impact) ')']); ylabel('Latitude [deg]'); ylim([-90 90]);
set(gca,'YTick',-90:45:90);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(5,:); % low
barAxis.Children(2).FaceColor = ORExLineColors(3,:); % med
barAxis.Children(3).FaceColor = ORExLineColors(8,:); % high
xvals = xlim; 
hold on; line([0 .99*xvals(2)],[180/pi*Lat0_all(1) 180/pi*Lat0_all(1)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('$\theta \geq 60^\circ$', '$30^\circ \leq \theta < 60^\circ$','$0^\circ \leq \theta < 30^\circ$','Launch Site');

Longf_el_stacked = [hist_long_el_high.Values; hist_long_el_med.Values; hist_long_el_low.Values]';
figure; bar(5:10:355,Longf_el_stacked,1,'stacked'); ylabel(['\# of Objects (out of ' num2str(tot_impact) ')']); xlabel('East Longitude [deg]'); xlim([0 360]);
set(gca,'XTick',0:90:360);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(5,:); % low
barAxis.Children(2).FaceColor = ORExLineColors(3,:); % med
barAxis.Children(3).FaceColor = ORExLineColors(8,:); % high
yvals = ylim; 
hold on; line([180/pi*Long0_all(1) 180/pi*Long0_all(1)],[0 .99*yvals(2)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('$\theta \geq 60^\circ$', '$30^\circ \leq \theta < 60^\circ$','$0^\circ \leq \theta < 30^\circ$','Launch Site');

figure;
subplot(3,1,1); hh = histogram2(Longf_el_high.*180/pi,Latf_el_high.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); title('$\theta \geq 60^\circ$'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
subplot(3,1,2); hm = histogram2(Longf_el_med.*180/pi,Latf_el_med.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); title('$30^\circ \leq \theta < 60^\circ$'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
subplot(3,1,3); hl = histogram2(Longf_el_low.*180/pi,Latf_el_low.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); title('$0^\circ \leq \theta < 30^\circ$'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
cmax = max([max(max(hh.Values)) max(max(hm.Values)) max(max(hl.Values))]);
subplot(3,1,1); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);
subplot(3,1,2); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);
subplot(3,1,3); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);

% azimuth bands
Inds_az_east = find(round(vel_az*180/pi)==0 | round(vel_az*180/pi)==30 | round(vel_az*180/pi)==330);
Inds_az_north = find(round(vel_az*180/pi)==60 | round(vel_az*180/pi)==90 | round(vel_az*180/pi)==120);
Inds_az_west = find(round(vel_az*180/pi)==150 | round(vel_az*180/pi)==180 | round(vel_az*180/pi)==210);
Inds_az_south = find(round(vel_az*180/pi)==240 | round(vel_az*180/pi)==270 | round(vel_az*180/pi)==300);
Longf_az_east = [];
Latf_az_east = [];
Longf_az_north = [];
Latf_az_north = [];
Longf_az_west = [];
Latf_az_west = [];
Longf_az_south = [];
Latf_az_south = [];
for jj = 1:length(vmags)
    sub = Longf_all(Inds_az_east,jj);
    Longf_az_east = [Longf_az_east; sub(~isnan(sub))];

    sub = Latf_all(Inds_az_east,jj);
    Latf_az_east = [Latf_az_east; sub(~isnan(sub))];    
    
    sub = Longf_all(Inds_az_north,jj);
    Longf_az_north = [Longf_az_north; sub(~isnan(sub))];

    sub = Latf_all(Inds_az_north,jj);
    Latf_az_north = [Latf_az_north; sub(~isnan(sub))]; 
    
    sub = Longf_all(Inds_az_west,jj);
    Longf_az_west = [Longf_az_west; sub(~isnan(sub))];

    sub = Latf_all(Inds_az_west,jj);
    Latf_az_west = [Latf_az_west; sub(~isnan(sub))]; 
    
    sub = Longf_all(Inds_az_south,jj);
    Longf_az_south = [Longf_az_south; sub(~isnan(sub))];

    sub = Latf_all(Inds_az_south,jj);
    Latf_az_south = [Latf_az_south; sub(~isnan(sub))];
    
end

figure; hold on; 
plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','s','Color',ORExLineColors(2,:),'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',ORExLineColors(2,:),'MarkerSize',18);
scatter(180/pi.*Longf_az_east,180/pi.*Latf_az_east,'filled','MarkerFaceColor',ORExLineColors(8,:));
scatter(180/pi.*Longf_az_north,180/pi.*Latf_az_north,'filled','MarkerFaceColor',ORExLineColors(3,:));
scatter(180/pi.*Longf_az_west,180/pi.*Latf_az_west,'filled','MarkerFaceColor',ORExLineColors(5,:));
scatter(180/pi.*Longf_az_south,180/pi.*Latf_az_south,'filled','MarkerFaceColor',ORExLineColors(7,:)); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]);
legend('Launch Site','East', 'North','West','South','Orientation','horizontal','Location','best')
set(gca,'XTick',0:90:360,'YTick',-90:45:90); 
grid on;

figure; hist_long_az_east = histogram(180/pi.*Longf_az_east, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_lat_az_east = histogram(180/pi.*Latf_az_east, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);
figure; hist_long_az_north = histogram(180/pi.*Longf_az_north, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_lat_az_north = histogram(180/pi.*Latf_az_north, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);
figure; hist_long_az_west = histogram(180/pi.*Longf_az_west, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_lat_az_west = histogram(180/pi.*Latf_az_west, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);
figure; hist_long_az_south = histogram(180/pi.*Longf_az_south, 0:10:360); xlabel('East Longitude [deg]'); xlim([0 360]);
figure; hist_lat_az_south = histogram(180/pi.*Latf_az_south, -90:10:90); xlabel('Latitude [deg]'); xlim([-90 90]);

Latf_az_stacked = [hist_lat_az_east.Values; hist_lat_az_north.Values; hist_lat_az_west.Values; hist_lat_az_south.Values]';
figure; barh(-85:10:85,Latf_az_stacked,1,'stacked'); xlabel(['\# of Objects (out of ' num2str(tot_impact) ')']); ylabel('Latitude [deg]'); ylim([-90 90]);
set(gca,'YTick',-90:45:90);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(7,:); % south
barAxis.Children(2).FaceColor = ORExLineColors(5,:); % west
barAxis.Children(3).FaceColor = ORExLineColors(3,:); % north
barAxis.Children(4).FaceColor = ORExLineColors(8,:); % east
xvals = xlim; 
hold on; line([0 .99*xvals(2)],[180/pi*Lat0_all(1) 180/pi*Lat0_all(1)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('East', 'North','West','South','Launch Site');

Longf_az_stacked = [hist_long_az_east.Values; hist_long_az_north.Values; hist_long_az_west.Values; hist_long_az_south.Values]';
figure; bar(5:10:355,Longf_az_stacked,1,'stacked'); ylabel(['\# of Objects (out of ' num2str(tot_impact) ')']); xlabel('East Longitude [deg]'); xlim([0 360]);
set(gca,'XTick',0:90:360);
barAxis = gca;
barAxis.Children(1).FaceColor = ORExLineColors(7,:); % south
barAxis.Children(2).FaceColor = ORExLineColors(5,:); % west
barAxis.Children(3).FaceColor = ORExLineColors(3,:); % north
barAxis.Children(4).FaceColor = ORExLineColors(8,:); % east
yvals = ylim; 
hold on; line([180/pi*Long0_all(1) 180/pi*Long0_all(1)],[0 .99*yvals(2)],'LineStyle','--','Color',ORExLineColors(1,:));
legend('East', 'North','West','South','Launch Site');

figure;
subplot(2,2,1); he = histogram2(Longf_az_east.*180/pi,Latf_az_east.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); title('East'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
subplot(2,2,2); hn = histogram2(Longf_az_north.*180/pi,Latf_az_north.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); xlim([0 360]); ylim([-90 90]); title('North'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
subplot(2,2,3); hw = histogram2(Longf_az_west.*180/pi,Latf_az_west.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); xlabel('East Longitude [deg]'); ylabel('Latitude [deg]'); xlim([0 360]); ylim([-90 90]); title('West'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
subplot(2,2,4); hs = histogram2(Longf_az_south.*180/pi,Latf_az_south.*180/pi,0:10:360,-90:10:90,'DisplayStyle','Tile'); colormap(ORExcolorMap); xlabel('East Longitude [deg]'); xlim([0 360]); ylim([-90 90]); title('South'); colorbar; set(gca,'XTick',0:90:360,'YTick',-90:45:90);
cmax = max([max(max(he.Values)) max(max(hn.Values)) max(max(hs.Values)) max(max(hw.Values))]);
subplot(2,2,1); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);
subplot(2,2,2); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);
subplot(2,2,3); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);
subplot(2,2,4); set(gca,'CLim',[0 cmax]); hold on; plot(180/pi*Long0_all(1),180/pi*Lat0_all(1),'Marker','x','Color',ORExLineColors(1,:),'MarkerSize',18);

%% histogram - line version - lifetime vs various parameters


% ORExcolorMap = ORExcolorMap(7:end,:);
% 
% set(groot, 'DefaultAxesColor',[38 38 38]./256, ...
% 'DefaultAxesXColor', [1 1 1], ...
% 'DefaultAxesYColor', [1 1 1], ...
% 'DefaultAxesZColor', [1 1 1], ...
% 'DefaultFigureColor',[38 38 38]./256, ...
% 'DefaultLegendTextColor',[1 1 1]);

% particle size
psize_all = repmat(radius_part,length(vmags),1);
rset = unique(radius_part(1:50));
Tend_all_reshape = reshape(Tend_all,1606*length(vmags),1);
figure; fig_rad_hist = gcf; hold on; 
ax = gca;
ax.ColorOrder = ORExcolorMap(ceil(linspace(1,32,22)),:); %jet(length(rset));
t_vs_radius = [];
for ii = 1:length(rset)
    inds = find(psize_all==rset(ii));
    figure; fig_temp = gcf; hrad = histogram(Tend_all_reshape(inds),[0:40],'Normalization','probability');
    temp = hrad.Values;
%     temp(temp==0) = .1;
    figure(fig_rad_hist); plot(.5:1:39.5,temp,'-o');
    close(fig_temp);
    t_vs_radius = [t_vs_radius; temp];
end
figure(fig_rad_hist); ylabel('$\%$ of 1606 Cases'); xlabel('Lifetime [days]');
legend(num2cell(num2str(rset),2));

t_vs_radius = t_vs_radius';
figure; bar(.5:1:39.5,t_vs_radius,1,'stacked'); ylabel('$\%$ of 1606 Cases'); xlabel('Lifetime [days]'); xlim([0 40]); title('Particle Radius');
barAxis = gca;
for ii = 1:length(barAxis.Children)
    barAxis.Children(ii).FaceColor = ax.ColorOrder(end-ii+1,:);
end

prad_escape = [];
prad_crash = [];
rset = unique(radius_part(1:50));
for ii = 1:length(vmags)
    prad_escape = [prad_escape; radius_part(Inds_direct_escape{ii}); radius_part(Inds_escape{ii})];
    prad_crash = [prad_crash; radius_part(Inds_suborbital{ii}); radius_part(Inds_orbital{ii})];
end
figure; hrad1 = histogram(prad_escape,[0, .4, .9:1:20.9]); %,'Normalization','probability');
figure; hrad2 = histogram(prad_crash,[0, .4, .9:1:20.9]); %,'Normalization','probability');
% figure; bar(.5:1:21.5,[hrad1.Values; hrad2.Values]',1,'stacked'); ylabel('$\%$ of 1606 Cases'); xlabel('Particle Radius [cm]'); xlim([0 22]); 
% barAxis = gca;
% barAxis.Children(1).FaceColor = ORExLineColors(3,:);
% barAxis.Children(2).FaceColor = ORExLineColors(5,:);
% legend('Escape','Crash')
% hold on;
% line([0 22],[401.5],'Color',ORExLineColors(1,:),'LineStyle','--')
figure; plot(rset, hrad1.Values./803.*100,'Color',ORExLineColors(3,:),'Marker','o'); xlabel('Particle Radius [cm]'); ylabel('$\%$ escaping'); 

% velocity magnitude
figure; fig_vmag_hist = gcf; hold on;
ax = gca;
ax.ColorOrder = ORExcolorMap(ceil(linspace(1,32,11)),:); %jet(length(vmags));
t_vs_velmag = [];
for ii = 1:length(vmags)
    inds = (ii-1)*1606 + 1:ii*1606;
    figure; fig_temp = gcf; hrad = histogram(Tend_all_reshape(inds),[0:40],'Normalization','probability');
    temp = hrad.Values;
%     temp(temp==0) = .1;
    figure(fig_vmag_hist); plot(.5:1:39.5,temp,'-o');
    close(fig_temp);
    t_vs_velmag = [t_vs_velmag; temp];
end
figure(fig_vmag_hist); ylabel('$\%$ of 1606 Cases'); xlabel('Lifetime [days]');
legend(num2cell(num2str(vmags'),2)');

t_vs_velmag = t_vs_velmag';
figure; bar(.5:1:39.5,t_vs_velmag,1,'stacked'); ylabel('$\%$ of 1606 Cases'); xlabel('Lifetime [days]'); xlim([0 40]); title('Velocity Magnitude')
barAxis = gca;
for ii = 1:length(barAxis.Children)
    barAxis.Children(ii).FaceColor = ax.ColorOrder(end-ii+1,:);
end

energy_escape = [];
energy_crash = [];
rho_p = 2e12; % 1, 2 g/cc
M_part = rho_p.*4/3*pi.*(radius_part.*1e-5).^3; % [kg]
for ii = 1:length(vmags)
    energy_escape = [energy_escape; M_part(Inds_direct_escape{ii}).*((vmags(ii)*1e-2)^2)./2; M_part(Inds_escape{ii}).*((vmags(ii)*1e-2)^2)./2];
    energy_crash = [energy_crash; M_part(Inds_suborbital{ii}).*((vmags(ii)*1e-2)^2)./2; M_part(Inds_orbital{ii}).*((vmags(ii)*1e-2)^2)./2];
end
energy_all = [energy_crash; energy_escape];
[Nenergy, edge_energy] = histcounts(energy_all,logspace(-8,1,100));
edge_use = [edge_energy(1) edge_energy(19) edge_energy(44) edge_energy(55) edge_energy(61) edge_energy(66) edge_energy(69) edge_energy(72) edge_energy(74) edge_energy(76) edge_energy(78:90) edge_energy(92) edge_energy(end)];
edge_use([7 9 12 14 17]) = [];
figure; hrad1 = histogram(energy_escape,edge_use); %,'Normalization','probability'); 4e-3  [10^-8 10^-6 5e-5 4e-4 2e-2 10^-1 1 3.5] [10^-8 10^-6 5e-5 tempspace(3) tempspace(6:end) 3.5]
figure; hrad2 = histogram(energy_crash,edge_use); %,'Normalization','probability');
figure; plot(hrad1.BinEdges(1:end-1), hrad1.Values./sum([hrad1.Values; hrad2.Values]).*100,'Color',ORExLineColors(3,:),'Marker','o'); xlabel('Energy [J]'); ylabel('$\%$ escaping'); set(gca,'XScale','log')

vel_escape = zeros(length(vmags),1);
for ii = 1:length(vmags)
    vel_escape(ii) = length(Inds_direct_escape{ii}) + length(Inds_escape{ii});
end
figure; plot(vmags, vel_escape./1606.*100,'Color',ORExLineColors(3,:),'Marker','o'); xlabel('Launch Velocity [cm/s]'); ylabel('$\%$ escaping');


% velocity az
az_all = repmat(vel_az,length(vmags),1);
azset = unique(vel_az);
figure; fig_az_hist = gcf; hold on; 
ax = gca;
ax.ColorOrder = ORExcolorMap(ceil(linspace(1,32,12)),:); %jet(length(azset));
for ii = 1:length(azset)
    inds = find(az_all==azset(ii));
    figure; fig_temp = gcf; hrad = histogram(Tend_all_reshape(inds),[0:40],'Normalization','probability');
    temp = hrad.Values;
%     temp(temp==0) = .1;
    figure(fig_az_hist); plot(.5:1:39.5,temp,'-o');
    close(fig_temp);
end
figure(fig_az_hist); ylabel('$\%$ of 1606 Cases'); xlabel('Lifetime [days]');
legend(num2cell(num2str(azset.*180/pi),2));

% velocity el
el_all = repmat(vel_el,length(vmags),1);
elset = unique(vel_el);
figure; fig_el_hist = gcf; hold on; 
ax = gca;
ax.ColorOrder = ORExcolorMap(ceil(linspace(1,32,7)),:); %jet(length(elset));
for ii = 1:length(elset)
    inds = find(el_all==elset(ii));
    figure; fig_temp = gcf; hrad = histogram(Tend_all_reshape(inds),[0:40],'Normalization','probability');
    temp = hrad.Values;
%     temp(temp==0) = .1;
    figure(fig_el_hist); plot(.5:1:39.5,temp,'-o');
    close(fig_temp);
end
figure(fig_el_hist); ylabel('$\%$ of 1606 Cases'); xlabel('Lifetime [days]');
legend(num2cell(num2str(elset.*180/pi),2));

%% Pull group of trajectories with same AMR
AMR_match = SRP_coeffs_part(1); % choose index 1-22, corresponding to particle radius
set_search = 'all'; %'short' '1to5' '5to40' '40plus'


Inds_matchAMR = cell(length(vmags),1);
switch set_search
    
    case 'all'
        Inds = find(SRP_coeffs_part == AMR_match);
        for ii = 1:length(vmags)
            Inds_matchAMR{ii} = Inds;
        end
        
    case 'short'
        
        Inds = find(Am_short == AMR_match./(1+.04*4/9).*1e6);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_short{ii})) - (startInd-1);
            Inds_matchAMR{ii} = Inds_short{ii}(Inds_sub);
            startInd = startInd+length(Inds_short{ii});
        end
        
        
    case '1to5'
        
        Inds = find(Am_1_5days == AMR_match./(1+.04*4/9).*1e6);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_1_5days{ii})) - (startInd-1);
            Inds_matchAMR{ii} = Inds_1_5days{ii}(Inds_sub);
            startInd = startInd+length(Inds_1_5days{ii});
        end
        
        
    case '5to40'
        
        Inds = find(Am_5_40days == AMR_match./(1+.04*4/9).*1e6);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_5_40days{ii})) - (startInd-1);
            Inds_matchAMR{ii} = Inds_5_40days{ii}(Inds_sub);
            startInd = startInd+length(Inds_5_40days{ii});
        end
        
        
    case '40plus'
        
        Inds = find(Am_max40days == AMR_match./(1+.04*4/9).*1e6);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_max40days{ii})) - (startInd-1);
            Inds_matchAMR{ii} = Inds_max40days{ii}(Inds_sub);
            startInd = startInd+length(Inds_max40days{ii});
        end
        
        
end

%% Pull group of trajectories with same initial velocity MAGNITUDE

% DON'T ACTUALLY HAVE TO DO THIS AS ITS JUST THE SAME AS INDICES ABOVE
% JUST CHANGING THE NAME HERE EFFECTIVELY
vmag_ind = 1; % choose index in vmags to match
set_search = 'all'; %'short' '1to5' '5to40' '40plus'


Inds_matchVmag = cell(length(vmags),1);
switch set_search
    
    case 'all'
        Inds_matchVmag{vmag_ind} = [1:1606]';
        
    case 'short'
        
        Inds_matchVmag{vmag_ind} = Inds_short{vmag_ind};
        
    case '1to5'
        
        Inds_matchVmag{vmag_ind} = Inds_1_5days{vmag_ind};
        
    case '5to40'
        
        Inds_matchVmag{vmag_ind} = Inds_5_40days{vmag_ind};
        
    case '40plus'
        
        Inds_matchVmag{vmag_ind} = Inds_max40days{vmag_ind};
        
end



%% Pull group of trajectories with same initial velocity DIRECTION
% ie matching vel az and el
velAz_match = vel_az(1); % choose 
velEl_match = vel_el(1); 

set_search = 'all'; %'short' '1to5' '5to40' '40plus'

Inds_matchVelDir = cell(length(vmags),1);
Inds_matchVelAz = cell(length(vmags),1);
Inds_matchVelEl = cell(length(vmags),1);
switch set_search
    
    case 'all'
        IndsAz = find(vel_az == velAz_match);
        IndsEl = find(vel_el == velEl_match);
        Inds = intersect(IndsAz, IndsEl);
        for ii = 1:length(vmags)
            Inds_matchVelDir{ii} = Inds;
            Inds_matchVelAz{ii} = IndsAz;
            Inds_matchVelEl{ii} = IndsEl;
        end
        
    case 'short'
        
        IndsAz = find(velaz_short == velAz_match);
        IndsEl = find(velel_short == velEl_match);
        Inds = intersect(IndsAz, IndsEl);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_short{ii})) - (startInd-1);
            Inds_matchVelDir{ii} = Inds_short{ii}(Inds_sub);
            Inds_sub = find(IndsAz>=startInd & IndsAz<=startInd-1+length(Inds_short{ii})) - (startInd-1);
            Inds_matchVelAz{ii} = Inds_short{ii}(Inds_sub);
            Inds_sub = find(IndsEl>=startInd & IndsEl<=startInd-1+length(Inds_short{ii})) - (startInd-1);
            Inds_matchVelEl{ii} = Inds_short{ii}(Inds_sub);
            startInd = startInd+length(Inds_short{ii});
        end
        
        
    case '1to5'
        
        IndsAz = find(velaz_1_5days == velAz_match);
        IndsEl = find(velel_1_5days == velEl_match);
        Inds = intersect(IndsAz, IndsEl);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_1_5days{ii})) - (startInd-1);
            Inds_matchVelDir{ii} = Inds_1_5days{ii}(Inds_sub);
            Inds_sub = find(IndsAz>=startInd & IndsAz<=startInd-1+length(Inds_1_5days{ii})) - (startInd-1);
            Inds_matchVelAz{ii} = Inds_1_5days{ii}(Inds_sub);
            Inds_sub = find(IndsEl>=startInd & IndsEl<=startInd-1+length(Inds_1_5days{ii})) - (startInd-1);
            Inds_matchVelEl{ii} = Inds_1_5days{ii}(Inds_sub);
            startInd = startInd+length(Inds_1_5days{ii});
        end
        
        
    case '5to40'
        
        IndsAz = find(velaz_5_40days == velAz_match);
        IndsEl = find(velel_5_40days == velEl_match);
        Inds = intersect(IndsAz, IndsEl);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_5_40days{ii})) - (startInd-1);
            Inds_matchVelDir{ii} = Inds_5_40days{ii}(Inds_sub);
            Inds_sub = find(IndsAz>=startInd & IndsAz<=startInd-1+length(Inds_5_40days{ii})) - (startInd-1);
            Inds_matchVelAz{ii} = Inds_5_40days{ii}(Inds_sub);
            Inds_sub = find(IndsEl>=startInd & IndsEl<=startInd-1+length(Inds_5_40days{ii})) - (startInd-1);
            Inds_matchVelEl{ii} = Inds_5_40days{ii}(Inds_sub);
            startInd = startInd+length(Inds_5_40days{ii});
        end
        
        
    case '40plus'
        
        IndsAz = find(velaz_max40days == velAz_match);
        IndsEl = find(velel_max40days == velEl_match);
        Inds = intersect(IndsAz, IndsEl);
        startInd = 1;
        for ii = 1:length(vmags)
            Inds_sub = find(Inds>=startInd & Inds<=startInd-1+length(Inds_max40days{ii})) - (startInd-1);
            Inds_matchVelDir{ii} = Inds_max40days{ii}(Inds_sub);
            Inds_sub = find(IndsAz>=startInd & IndsAz<=startInd-1+length(Inds_max40days{ii})) - (startInd-1);
            Inds_matchVelAz{ii} = Inds_max40days{ii}(Inds_sub);
            Inds_sub = find(IndsEl>=startInd & IndsEl<=startInd-1+length(Inds_max40days{ii})) - (startInd-1);
            Inds_matchVelEl{ii} = Inds_max40days{ii}(Inds_sub);
            startInd = startInd+length(Inds_max40days{ii});
        end
        
end



%% Pull group of trajectories with same initial velocity vector or same v and az or el
% So this is just one vmag set of the velocity direction case
vmag_ind = 1; % choose index in vmags to match

Inds_matchVvec = Inds_matchVelDir{vmag_ind};
Inds_matchVAz = Inds_matchVelAz{vmag_ind};
Inds_matchVEl = Inds_matchVelEl{vmag_ind};


%% Initial condition OEs
sma_init = zeros(1606,length(vmags));
ecc_init = zeros(1606,length(vmags));
inc_init = zeros(1606,length(vmags));
argp_init = zeros(1606,length(vmags));
RAAN_init = zeros(1606,length(vmags));
nu_init = zeros(1606,length(vmags));
edtoArgX_init = zeros(1606,length(vmags));
sma_init_sub = [];
sma_init_orb = [];
sma_init_dir = [];
sma_init_esc = [];

ecc_init_sub = [];
ecc_init_orb = [];
ecc_init_dir = [];
ecc_init_esc = [];

inc_init_sub = [];
inc_init_orb = [];
inc_init_dir = [];
inc_init_esc = [];

argp_init_sub = [];
argp_init_orb = [];
argp_init_dir = [];
argp_init_esc = [];

RAAN_init_sub = [];
RAAN_init_orb = [];
RAAN_init_dir = [];
RAAN_init_esc = [];

nu_init_sub = [];
nu_init_orb = [];
nu_init_dir = [];
nu_init_esc = [];

edotX_init_sub = [];
edotX_init_orb = [];
edotX_init_dir = [];
edotX_init_esc = [];
if caseNum == 20 || caseNum == 21
    etStr = '2019 JAN 06 20:58:28.000 UTC';
elseif caseNum == 19 || caseNum == 23
    etStr = '2019 Jan 19 00:53:40.900 UTC';
elseif caseNum == 22
    etStr = '2019 Feb 11 23:27:28.480 UTC';
else
    error('Bad case number');
end
for ii = 1:length(vmags)
    load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(ii)) 'cms_lowTol' caseStr],'T','X');
    for jj = 1:73
        currInd = (jj-1)*22 + 1;
        [oe_init, ~, ~, edotX] = ECI2OE_Bennu(X(1,(currInd-1)*6 +1:currInd*6), T(1), etStr);
        sma_init(currInd:currInd+21,ii) = oe_init(1).*ones(22,1);
        ecc_init(currInd:currInd+21,ii) = oe_init(2).*ones(22,1);
        inc_init(currInd:currInd+21,ii) = oe_init(3).*ones(22,1);
        argp_init(currInd:currInd+21,ii) = oe_init(4).*ones(22,1);
        RAAN_init(currInd:currInd+21,ii) = oe_init(5).*ones(22,1);
        nu_init(currInd:currInd+21,ii) = oe_init(6).*ones(22,1);
        edtoArgX_init(currInd:currInd+21,ii) = edotX.*ones(22,1);
    end
    sma_init_sub = [sma_init_sub; sma_init(Inds_suborbital{ii},ii)];
    sma_init_orb = [sma_init_orb; sma_init(Inds_orbital{ii},ii)];
    sma_init_dir = [sma_init_dir; sma_init(Inds_direct_escape{ii},ii)];
    sma_init_esc = [sma_init_esc; sma_init(Inds_escape{ii},ii)];
    
    ecc_init_sub = [ecc_init_sub; ecc_init(Inds_suborbital{ii},ii)];
    ecc_init_orb = [ecc_init_orb; ecc_init(Inds_orbital{ii},ii)];
    ecc_init_dir = [ecc_init_dir; ecc_init(Inds_direct_escape{ii},ii)];
    ecc_init_esc = [ecc_init_esc; ecc_init(Inds_escape{ii},ii)];
    
    inc_init_sub = [inc_init_sub; inc_init(Inds_suborbital{ii},ii)];
    inc_init_orb = [inc_init_orb; inc_init(Inds_orbital{ii},ii)];
    inc_init_dir = [inc_init_dir; inc_init(Inds_direct_escape{ii},ii)];
    inc_init_esc = [inc_init_esc; inc_init(Inds_escape{ii},ii)];
    
    argp_init_sub = [argp_init_sub; argp_init(Inds_suborbital{ii},ii)];
    argp_init_orb = [argp_init_orb; argp_init(Inds_orbital{ii},ii)];
    argp_init_dir = [argp_init_dir; argp_init(Inds_direct_escape{ii},ii)];
    argp_init_esc = [argp_init_esc; argp_init(Inds_escape{ii},ii)];
    
    RAAN_init_sub = [RAAN_init_sub; RAAN_init(Inds_suborbital{ii},ii)];
    RAAN_init_orb = [RAAN_init_orb; RAAN_init(Inds_orbital{ii},ii)];
    RAAN_init_dir = [RAAN_init_dir; RAAN_init(Inds_direct_escape{ii},ii)];
    RAAN_init_esc = [RAAN_init_esc; RAAN_init(Inds_escape{ii},ii)];
    
    nu_init_sub = [nu_init_sub; nu_init(Inds_suborbital{ii},ii)];
    nu_init_orb = [nu_init_orb; nu_init(Inds_orbital{ii},ii)];
    nu_init_dir = [nu_init_dir; nu_init(Inds_direct_escape{ii},ii)];
    nu_init_esc = [nu_init_esc; nu_init(Inds_escape{ii},ii)];
    
    edotX_init_sub = [edotX_init_sub; edtoArgX_init(Inds_suborbital{ii},ii)];
    edotX_init_orb = [edotX_init_orb; edtoArgX_init(Inds_orbital{ii},ii)];
    edotX_init_dir = [edotX_init_dir; edtoArgX_init(Inds_direct_escape{ii},ii)];
    edotX_init_esc = [edotX_init_esc; edtoArgX_init(Inds_escape{ii},ii)];
    
    
end
% figure; plot(reshape(Tend_all,1606*length(vmags),1),reshape(sma_init,1606*length(vmags),1),'o'); xlabel('tf [days]'); ylabel('$a$ [km]','Interpreter','latex'); xlim([0 40]);
% figure; plot(reshape(Tend_all,1606*length(vmags),1),reshape(ecc_init,1606*length(vmags),1),'o'); xlabel('tf [days]'); ylabel('$e$ [-]','Interpreter','latex'); xlim([0 40]);
% figure; plot(reshape(Tend_all,1606*length(vmags),1),reshape(inc_init,1606*length(vmags),1).*180/pi,'o'); xlabel('tf [days]'); ylabel('$i$ [deg]','Interpreter','latex'); xlim([0 40]);
% figure; plot(reshape(Tend_all,1606*length(vmags),1),reshape(argp_init,1606*length(vmags),1).*180/pi,'o'); xlabel('tf [days]'); ylabel('$\omega$ [deg]','Interpreter','latex'); xlim([0 40]);
% figure; plot(reshape(Tend_all,1606*length(vmags),1),reshape(RAAN_init,1606*length(vmags),1).*180/pi,'o'); xlabel('tf [days]'); ylabel('$\Omega$ [deg]','Interpreter','latex'); xlim([0 40]);
% 
% figure; scatter3(reshape(sma_init,1606*length(vmags),1), reshape(ecc_init,1606*length(vmags),1), reshape(Tend_all,1606*length(vmags),1));

figure; hold on; 
t_temp = [];
for ii = 1:11
    t_temp = [t_temp; Tend_all(Inds_suborbital{ii},ii)];
end
plot(t_temp, edotX_init_sub, 'Marker','o','MarkerSize',18,'MarkerFaceColor',ORExLineColors(7,:),'LineStyle','none','LineWidth',1,'Color',ORExLineColors(7,:)); 
t_temp = [];
for ii = 1:11
    t_temp = [t_temp; Tend_all(Inds_orbital{ii},ii)];
end
plot(t_temp, edotX_init_orb, 'Marker','^','MarkerSize',15,'MarkerFaceColor',ORExLineColors(5,:),'LineWidth',0.5,'LineStyle','none','Color',ORExLineColors(5,:)); 
t_temp = [];
for ii = 1:11
    t_temp = [t_temp; Tend_all(Inds_direct_escape{ii},ii)];
end
plot(t_temp, edotX_init_dir, 'Marker','+','MarkerSize',12,'Color',ORExLineColors(2,:),'LineStyle','none'); 
t_temp = [];
for ii = 1:11
    t_temp = [t_temp; Tend_all(Inds_escape{ii},ii)];
end
plot(t_temp, edotX_init_esc, 'Marker','x','MarkerSize',12,'Color',ORExLineColors(4,:),'LineStyle','none'); 
xlabel('tf [days]'); ylabel('edotX [deg]','Interpreter','latex'); xlim([0 40]);
legend('Suborbital','Orbital','Direct Escape','Escape')

% a vs e
figure; hold on; plot(sma_init_sub, ecc_init_sub, 'Marker','o','MarkerSize',18,'MarkerFaceColor',ORExLineColors(7,:),'LineStyle','none','LineWidth',1,'Color',ORExLineColors(7,:)); plot(sma_init_orb, ecc_init_orb, 'Marker','^','MarkerSize',15,'MarkerFaceColor',ORExLineColors(5,:),'LineWidth',0.5,'LineStyle','none','Color',ORExLineColors(5,:)); plot(sma_init_dir, ecc_init_dir, 'Marker','+','MarkerSize',12,'Color',ORExLineColors(2,:),'LineStyle','none'); plot(sma_init_esc, ecc_init_esc, 'Marker','x','MarkerSize',12,'Color',ORExLineColors(4,:),'LineStyle','none'); 
xlabel('a [km]'); ylabel('e [-]');
legend('Suborbital','Orbital','Direct Escape','Escape')

% r_p vs e
figure; hold on; plot(sma_init_sub.*(1-ecc_init_sub), ecc_init_sub, 'Marker','o','MarkerSize',18,'MarkerFaceColor',ORExLineColors(7,:),'LineStyle','none','LineWidth',1,'Color',ORExLineColors(7,:)); plot(sma_init_orb.*(1-ecc_init_orb), ecc_init_orb, 'Marker','^','MarkerSize',15,'MarkerFaceColor',ORExLineColors(5,:),'LineWidth',0.5,'LineStyle','none','Color',ORExLineColors(5,:)); plot(sma_init_dir.*(1-ecc_init_dir), ecc_init_dir, 'Marker','+','MarkerSize',12,'Color',ORExLineColors(2,:),'LineStyle','none'); plot(sma_init_esc.*(1-ecc_init_esc), ecc_init_esc, 'Marker','x','MarkerSize',12,'Color',ORExLineColors(4,:),'LineStyle','none'); 
xlabel('$r_p$ [km]'); ylabel('e [-]');
legend('Suborbital','Orbital','Direct Escape','Escape')

% Om vs inc
figure; hold on; plot(RAAN_init_sub.*180./pi, inc_init_sub.*180./pi, 'Marker','o','MarkerSize',18,'MarkerFaceColor',ORExLineColors(7,:),'LineStyle','none','LineWidth',1,'Color',ORExLineColors(7,:)); plot(RAAN_init_orb.*180./pi, inc_init_orb.*180./pi, 'Marker','^','MarkerSize',15,'MarkerFaceColor',ORExLineColors(5,:),'LineWidth',0.5,'LineStyle','none','Color',ORExLineColors(5,:)); plot(RAAN_init_dir.*180./pi, inc_init_dir.*180./pi, 'Marker','+','MarkerSize',12,'Color',ORExLineColors(2,:),'LineStyle','none'); plot(RAAN_init_esc.*180./pi, inc_init_esc.*180./pi, 'Marker','x','MarkerSize',12,'Color',ORExLineColors(4,:),'LineStyle','none'); 
xlabel('$\Omega$ [deg]'); ylabel('i [deg]');
legend('Suborbital','Orbital','Direct Escape','Escape')

% a vs argp
figure; hold on; plot(nu_init_sub.*180./pi, argp_init_sub.*180./pi, 'Marker','o','MarkerSize',18,'MarkerFaceColor',ORExLineColors(7,:),'LineStyle','none','LineWidth',1,'Color',ORExLineColors(7,:)); plot(nu_init_orb.*180./pi, argp_init_orb.*180./pi, 'Marker','^','MarkerSize',15,'MarkerFaceColor',ORExLineColors(5,:),'LineWidth',0.5,'LineStyle','none','Color',ORExLineColors(5,:)); plot(nu_init_dir.*180./pi, argp_init_dir.*180./pi, 'Marker','+','MarkerSize',12,'Color',ORExLineColors(2,:),'LineStyle','none'); plot(nu_init_esc.*180./pi, argp_init_esc.*180./pi, 'Marker','x','MarkerSize',12,'Color',ORExLineColors(4,:),'LineStyle','none'); 
xlabel('$\nu$ [deg]'); ylabel('$\omega$ [deg]');
legend('Suborbital','Orbital','Direct Escape','Escape')

% a vs edotX (edotX < 0 gives eccentricity decrease at apoapse)
figure; hold on; plot(nu_init_sub.*180./pi, edotX_init_sub, 'Marker','o','MarkerSize',18,'MarkerFaceColor',ORExLineColors(7,:),'LineStyle','none','LineWidth',1,'Color',ORExLineColors(7,:)); plot(nu_init_orb.*180./pi, edotX_init_orb, 'Marker','^','MarkerSize',15,'MarkerFaceColor',ORExLineColors(5,:),'LineWidth',0.5,'LineStyle','none','Color',ORExLineColors(5,:)); plot(nu_init_dir.*180./pi, edotX_init_dir, 'Marker','+','MarkerSize',12,'Color',ORExLineColors(2,:),'LineStyle','none'); plot(nu_init_esc.*180./pi, edotX_init_esc, 'Marker','x','MarkerSize',12,'Color',ORExLineColors(4,:),'LineStyle','none'); 
xlabel('$\nu$ [Deg]'); ylabel('edotX [-]');
legend('Suborbital','Orbital','Direct Escape','Escape')


%% OEs: Get trajectory information for sets of interest

% 40 day orbits only to start - will have to shorten time spans for other
% orbits in future
figure; fig_sma = gcf; hold on; xlabel('t [days]'); ylabel('$a$ [km]','Interpreter','latex'); xlim([0 10]);
figure; fig_ecc = gcf; hold on;  xlabel('t [days]'); ylabel('$e$ [-]','Interpreter','latex'); xlim([0 10]);
figure; fig_inc = gcf; hold on;  xlabel('t [days]'); ylabel('$i$ [rad]','Interpreter','latex'); xlim([0 10]);
figure; fig_w = gcf; hold on;  xlabel('t [days]'); ylabel('$\omega$ [rad]','Interpreter','latex'); xlim([0 10]);
figure; fig_Om = gcf; hold on;  xlabel('t [days]'); ylabel('$\Omega$ [rad]','Interpreter','latex'); xlim([0 40]);
figure; fig_rp = gcf; hold on; xlabel('t [days]'); ylabel('$r_p$ [km]','Interpreter','latex'); xlim([0 10]);
% figure; fig_ra = gcf; hold on; xlabel('t [days]'); ylabel('$r_a$ [km]','Interpreter','latex'); xlim([0 40]);
% figure; fig_r = gcf; hold on; xlabel('t [days]'); ylabel('$r$ [km]','Interpreter','latex'); xlim([0 40]);
% figure; fig_hx = gcf; hold on; xlabel('t [days]'); ylabel('$J_2$, $J_3$ balance $< 0$','Interpreter','latex'); xlim([0 40]);
% figure; fig_edx = gcf; hold on; xlabel('t [days]'); ylabel('$\dot{e} <$ 0','Interpreter','latex'); xlim([0 10]);
% figure; fig_plane = gcf; hold on; xlabel('$\Omega$ [deg]','Interpreter','latex'); ylabel('$i$ [deg]','Interpreter','latex'); 
counter = 1;
count_to_v_case = zeros(0,2);
tlastind_vec = [];
if caseNum == 20 || caseNum == 21
    etStr = '2019 JAN 06 20:58:28.000 UTC';
elseif caseNum == 19 || caseNum == 23
    etStr = '2019 Jan 19 00:53:40.900 UTC';
elseif caseNum == 22
    etStr = '2019 Feb 11 23:27:28.480 UTC';
else
    error('Bad case number');
end
for ii = 1:length(vmags)
    load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(ii)) 'cms_lowTol' caseStr],'T','X');
%     Inds = Inds_max40days{ii};
    Inds = Inds_orbital{ii};
    for jj = 1:length(Inds)
        if Tend_all(Inds(jj),ii) < 1
            continue
        end
        count_to_v_case = [count_to_v_case; ii Inds(jj)];
        
        tlastind = find(T(1:961)./(3600*24)==Tend_all(Inds(jj),ii),1);
        if isempty(tlastind)
            tlastind = 961;
        end
        tlastind_vec = [tlastind_vec; tlastind];
        
        [oe_out{counter}, ~, ~, edotArgX_all{counter}] = ECI2OE_Bennu(X(1:tlastind,(Inds(jj)-1)*6 +1:Inds(jj)*6), T(1:tlastind), etStr);
        
%         figure(fig_sma); plot(T(1:tlastind)./(3600*24),oe_out{counter}(:,1));
%         figure(fig_ecc); plot(T(1:tlastind)./(3600*24),oe_out{counter}(:,2));
%         figure(fig_inc); plot(T(1:tlastind)./(3600*24),oe_out{counter}(:,3));
%         figure(fig_w); plot(T(1:tlastind)./(3600*24),oe_out{counter}(:,4));
%         figure(fig_Om); plot(T./(3600*24),oe_out{counter}(:,5));
%         figure(fig_rp); plot(T(1:tlastind)./(3600*24), oe_out{counter}(:,1).*(1-oe_out{counter}(:,2)));
%         figure(fig_ra); plot(T./(3600*24), oe_out{counter}(:,1).*(1+oe_out{counter}(:,2)));
%         figure(fig_r); plot(T./(3600*24), rmag_all{counter}(:,1));
%         figure(fig_hx); plot(T./(3600*24), hhat_all{counter}(:,1));
%         figure(fig_edx); plot(T(1:tlastind)./(3600*24), edotArgX_all{counter}(:,1));
%         figure(fig_plane); plot(oe_out{counter}(:,5), oe_out{counter}(:,3), '-o')
        counter = counter + 1;
    end
end

figure(fig_sma);
for ii = 1:size(count_to_v_case,1)
    plot(T(1:tlastind_vec(ii))./(3600*24),oe_out{ii}(:,1));
end
figure(fig_ecc);
for ii = 1:size(count_to_v_case,1)
    plot(T(1:tlastind_vec(ii))./(3600*24),oe_out{ii}(:,2));
end
figure(fig_inc);
for ii = 1:size(count_to_v_case,1)
    plot(T(1:tlastind_vec(ii))./(3600*24),oe_out{ii}(:,3).*180/pi);
end
figure(fig_w);
for ii = 1:size(count_to_v_case,1)
    plot(T(1:tlastind_vec(ii))./(3600*24),oe_out{ii}(:,4).*180/pi);
end
figure(fig_Om);
for ii = 1:size(count_to_v_case,1)
    plot(T(1:tlastind_vec(ii))./(3600*24),oe_out{ii}(:,5).*180/pi);
end
figure(fig_rp);
for ii = 1:size(count_to_v_case,1)
    plot(T(1:tlastind_vec(ii))./(3600*24),oe_out{ii}(:,1).*(1-oe_out{ii}(:,2)));
end


%% Periapse/apoapse for final paper cases after getting rmag_all
for ii = 1:length(rmag_all)
    drmag = diff(rmag_all{ii});
    sdrmag = sign(drmag);
    dsdrmag = diff(sdrmag);
    indp_all{ii} = find(dsdrmag==2)+1;
    inda_all{ii} = find(dsdrmag==-2)+1;    
end

plotCase = 1543;
figure; hold on; plot(T./3600/24,rmag_all{plotCase}); plot(T(indp_all{plotCase})./3600/24,rmag_all{plotCase}(indp_all{plotCase}),'o'); plot(T(inda_all{plotCase})./3600/24,rmag_all{plotCase}(inda_all{plotCase}),'o');

% compute time information (all of these cases have the same T vector)
load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(1)) 'cms_lowTol' caseStr],'T');
figure; hold on; xlabel('Case Index'); ylabel('Periods [days]');
for ii = 1:length(rmag_all)

    T_periapse{ii} = T(indp_all{ii});
    T_apoapse{ii} = T(inda_all{ii});
    Periods_periapse_40days{ii} = diff(T_periapse{ii});
    Periods_apoapse_40days{ii} = diff(T_apoapse{ii});
    if ~isempty(Periods_periapse_40days{ii})
        plot(ii.*ones(length(Periods_periapse_40days{ii}),1),Periods_periapse_40days{ii}./3600/24,'o');
    end
    if ~isempty(Periods_apoapse_40days{ii})
        plot(ii.*ones(length(Periods_apoapse_40days{ii}),1),Periods_apoapse_40days{ii}./3600/24,'o');
    end

end



%% Plot a given 40 day case from original data
counterInd = 63;
load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(count_to_v_case(counterInd,1))) 'cms_lowTol' caseStr],'T','X');
figure; plot3(X(:,(count_to_v_case(counterInd,2)-1)*6 + 1),X(:,(count_to_v_case(counterInd,2)-1)*6 + 2),X(:,(count_to_v_case(counterInd,2)-1)*6 + 3))





%% Make a movie of trajectories in sun relative frame

figure; hold on; 
currentcase = 0; 
for ii = 1:length(count_to_v_case)
    if count_to_v_case(ii,1) ~= currentcase
        disp(['ii = ' num2str(ii) ' and currentcase before = ' num2str(currentcase)]);
        load(['Particles_p9AU_12kPalmer_rho1190_Vel_' num2str(vmags(count_to_v_case(ii,1))) 'cms_lowTol' caseStr],'T','X');
        currentcase = count_to_v_case(ii,1);
    end
    X_Sun_rot{ii} = ACI2SunRot_Bennu(X(:,(count_to_v_case(ii,2)-1)*6 +1:count_to_v_case(ii,2)*6), T, etStr);
    T_Sun_rot{ii} = T;
%     plot3(X_Sun_rot{ii}(:,1),X_Sun_rot{ii}(:,2),X_Sun_rot{ii}(:,3))
end
% 
% escapes = [];
% for ii = 1:length(rmag_all)
%     if rmag_all{ii}(end) > 15
%         escapes = [escapes; ii];
%     end
% end

% escapes = [8     9    33   118   119   120   121   122   123]; 
% escapes = [119   120   121   122   123]; 

ORExLineColors = [...
%     0 0 0; ...          % 1 black
    256 256 256; ...
    210 17 71; ...      % 2 red
    255 135 65; ...     % 3 orange
    255 223 90; ...     % 4 yellow
    113 210 162; ...    % 5 light green
    64 150 181; ...     % 6 light blue
    100 85 165; ...     % 7 light purple
    ]./256; 
%     12 51 114]./256;    % 8 indigo

ShapeModelFile = '../Shapes - 75cm delivery/Palmer/g_12580mm_spc_obj_0000n00000_v020_noComments.obj';
Polygon = GenerateVertices(ShapeModelFile);

period_half_count = zeros(length(count_to_v_case),1);
for jj = 1:length(count_to_v_case)
    period_half_count(jj) = floor(2*pi*sqrt(oe_out{jj}(1,1)^3/4.892e-9)/2/900); % half a period, in steps of 900s
end

for ii = 1:length(T)
    figH = figure; hold on;
    patch('Faces', Polygon.FacetIndex, 'Vertices', Polygon.Vertex(:,2:4),'FaceColor',[.8 .8 .8],'LineStyle','none')
    axis equal; grid on; set(figH, 'Color', [38 38 38]./256); xlabel('$\hat{x}$ [km]'); ylabel('$\hat{y}$ [km]'); zlabel('$\hat{z}$ [km]'); 
    set(gca, 'Color', [38 38 38]./256,'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
    line([0 0],[0 0],[0 5],'Color','w'); %,'LineWidth',6);
    line([-5 0],[0 0],[0 0],'Color','y'); %,'LineWidth',6);
    lgd = legend('Bennu','$\omega_B$','to Sun');
    lgd.TextColor = [1 1 1];
    title(sprintf('Time = %4.2f days',T(ii)/3600/24));
    xticks([-8:2:8]);
    yticks([-8:2:8]);
    zticks([-8:2:8]);
    
    for jj = 1:length(count_to_v_case)
        if ~isempty(find(escapes==jj,1))
            continue
        end
        if ii < period_half_count(jj)
            plot3(X_Sun_rot{jj}(1:ii,1),X_Sun_rot{jj}(1:ii,2),X_Sun_rot{jj}(1:ii,3),'Color',ORExLineColors(mod(jj,7)+1,:));
            plot3(X_Sun_rot{jj}(ii,1),X_Sun_rot{jj}(ii,2),X_Sun_rot{jj}(ii,3),'Marker','o','MarkerSize',10,'Color',ORExLineColors(mod(jj,7)+1,:),'MarkerFaceColor',ORExLineColors(mod(jj,7)+1,:));
        else
            plot3(X_Sun_rot{jj}(ii-period_half_count(jj)+1:ii,1),X_Sun_rot{jj}(ii-period_half_count(jj)+1:ii,2),X_Sun_rot{jj}(ii-period_half_count(jj)+1:ii,3),'Color',ORExLineColors(mod(jj,7)+1,:));
            plot3(X_Sun_rot{jj}(ii,1),X_Sun_rot{jj}(ii,2),X_Sun_rot{jj}(ii,3),'Marker','o','MarkerSize',10,'Color',ORExLineColors(mod(jj,7)+1,:),'MarkerFaceColor',ORExLineColors(mod(jj,7)+1,:));
        end
    end
    
    
    xlim([-8 8]); ylim([-8 8]); zlim([-8 8]);
    
    drawnow;
    
    view(-130,20);
    Regframes(ii) = getframe(figH);
    
    view(0,90);
    Nframes(ii) = getframe(figH);
    
    view(60,-20);
    Oppframes(ii) = getframe(figH);
    
    close(figH);
    
end

% make movie
writerObj2 = VideoWriter('test_reg_v2.mp4','MPEG-4'); %[nameStr '_poleLong0_AngleTargeting.mp4']);
writerObj2.FrameRate=30;
writerObj2.Quality = 100;
open(writerObj2);
for i = 1:length(Regframes)
    writeVideo(writerObj2,Regframes(i));
end
close(writerObj2);

writerObj2 = VideoWriter('test_north_v2.mp4','MPEG-4'); %[nameStr '_poleLong0_AngleTargeting.mp4']);
writerObj2.FrameRate=30;
writerObj2.Quality = 100;
open(writerObj2);
for i = 1:length(Nframes)
    writeVideo(writerObj2,Nframes(i));
end
close(writerObj2);

writerObj2 = VideoWriter('test_opp_v2.mp4','MPEG-4'); %[nameStr '_poleLong0_AngleTargeting.mp4']);
writerObj2.FrameRate=30;
writerObj2.Quality = 100;
open(writerObj2);
for i = 1:length(Oppframes)
    writeVideo(writerObj2,Oppframes(i));
end
close(writerObj2);