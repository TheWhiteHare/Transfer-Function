function [CMRO2_stim] = calculate_CMRO2()

mainFilenamePath = 'F:\Vessel 2D SS 0-20Perc PV_VL_Un\Master Plots Oxygen NO Coupling\Movies and Figures\To_Davis\';
filenames = {'oxygen_LTA','Drug_pre_post'};
probe_depth = {'DV100','DV300','DV500','DV800'}; 
region = {'FC','PC'};

figure(004)
for region_index = [1 2];


ii = 1;
filepath = ([mainFilenamePath filenames{ii}]);
[oxygen_LTA] = importdata([filepath '.mat'], ' ', 8);
t = oxygen_LTA.(region{region_index}).time;

pp = 2;
hold_data = oxygen_LTA.(region{region_index}).(probe_depth{pp});
[a] = isnan(hold_data);
remove_row = find(a(:,1) == 1);
oxygen_LTA.(region{region_index}).(probe_depth{pp})(remove_row,:) = [];

threshold = 3;

for pp = 1:length(probe_depth)
    subplot(2,4,pp+(region_index-1)*4),hold on
    ii = 1;
    for kk = 1:size(oxygen_LTA.(region{region_index}).(probe_depth{pp}),1)
        O2{kk} = oxygen_LTA.(region{region_index}).(probe_depth{pp})(kk,:);
        baseline = mean(O2{kk}(1:91));
        %if any(threshold<=(O2{kk}-mean(O2{kk}(1:91))))
            h1 = plot(t,O2{kk}-baseline,'Color',[0 0 0 0.25]);
            oxygen.(region{region_index}).all{pp}(ii,:) = O2{kk};
            ii = ii + 1;
        %else
        %end
        
    end
    if any(any(isnan(oxygen.(region{region_index}).all{pp})))
    oxygen.(region{region_index}).all{pp}(isnan(oxygen.(region{region_index}).all{pp})) = [];
    else
    end
    h2 = plot(t,mean(oxygen.(region{region_index}).all{pp})-mean(oxygen.(region{region_index}).all{pp}(1:91)),'k','LineWidth',3);
    axis([-3 10 -2 12])
    title([region{region_index} ' all responses: ' probe_depth{pp}])
    xlabel('time [s]')
    ylabel('\DeltaPO_{2}')
    %h3 = plot(oxygen_data{region_index}(:,1)-buffer,oxygen_data{region_index}(:,2)-oxygen_data{region_index}(1,2),'Color',[0 0 1],'LineWidth',3);
    legend([h1 h2], 'individual trial','average','COMSOL')
end
end


for region_index = [1 2];
for pp = 1:4

switch region_index
    case 1 %FC
        M = 0.95 %5% constriction
        PO2_increase = 5.7; %mmHg
    case 2 %FL/HL
        M = 1.1 %10% dilation
        PO2_increase = 5.7; %mmHg
    otherwise
end

hold_data = mean(oxygen.(region{region_index}).all{pp})
hold_time = find(hold_data(91:) == max(hold_data))
stim_O2 = mean(oxygen.(region{region_index}).all{pp}(:,hold_time),2);
baseline_O2 = mean(oxygen.(region{region_index}).all{pp}(:,1:91),2);

CMRO2_stim{region_index,pp} = extract_CMRO2(baseline_O2, stim_O2, PO2_increase, M);

end
end


delta_CMRO2 = cellfun(@(x) x(2,:)-x(1,:),CMRO2_stim,'UniformOutput',0);

depth = [100 300 500 800];
figure, hold on
for pp = 1:4;
[h_FC] = scatter(delta_CMRO2{1,pp},-depth(pp).*ones(length(delta_CMRO2{1,pp}),1),'c')
[h_FLHL] = scatter(delta_CMRO2{2,pp},-depth(pp).*ones(length(delta_CMRO2{2,pp}),1),'m')
legend([h_FC h_FLHL],'FC','FL/HL')
axis([-0.05 0.05 -1000 0])
end

depth = [100 300 500 800];
h_FC = [];
h_FC_depth = [];
h_FLHL = [];
h_FLHL_depth = [];
for pp = 1:4;
[h_FC] = [h_FC; delta_CMRO2{1,pp}'];
[h_FC_depth] = [h_FC_depth; depth(pp).*ones(length(delta_CMRO2{1,pp}'),1)];
[h_FLHL] = [h_FLHL; delta_CMRO2{2,pp}'];
[h_FLHL_depth] = [h_FLHL_depth; depth(pp).*ones(length(delta_CMRO2{2,pp}'),1)];
end
figure, hold on
[h_FC_box] = boxplot(h_FC,h_FC_depth,'Colors','c')
[h_FLHL_box] = boxplot(h_FLHL,h_FLHL_depth,'Colors','m')
ylabel('\DeltaCMRO_{2} [\mumole/cm^{3} min]')
xlabel('depth [\mum]')
title('\Delta CMRO_{2} in FC and FL/HL')


end
% CMRO2_running = 0.6/60*0.6;
% M = 1.1;
% PO2_increase = 5.7;
% 
% PO2_artery = [40] %mmHg
% CMRO2 = 3/60; %M s^-1
% D = 2.8*10^-9; %m^2 s^-1
% alpha = 0.00139; %M mmHg^-1
% R_artery = 0.9*10^-5; %m
% Rt = 0.5*10^-4; %m
% dr = 10^-6;
% 
% 
% r = [R_artery:dr:Rt]; %m
% PO2_rest = PO2_artery + CMRO2/(4*alpha*D).*(r.^2 - R_artery^2) - CMRO2/(2*D*alpha)*Rt^2.*log(r./R_artery)
% 
% figure, hold on
% h1 = plot(r.*10^6,PO2_rest,'Color',[0 0 1],'LineWidth',2)
% 
% % now with dilation and increased CMRO2
% 
% 
% PO2_artery = PO2_artery + PO2_increase; %mmHg
% CMRO2 = 3/60; %M s^-1
% CMRO2 = CMRO2+CMRO2_running;
% D = 2.8*10^-9; %m^2 s^-1
% alpha = 0.00139; %M mmHg^-1
% R_artery = R_artery*M; %m
% Rt = 0.5*10^-4; %m
% dr = 10^-6;
% 
% 
% r = [R_artery:dr:Rt]; %m
% PO2_rest = PO2_artery + CMRO2/(4*alpha*D).*(r.^2 - R_artery^2) - CMRO2/(2*D*alpha)*Rt^2.*log(r./R_artery)
% 
% h2 = plot(r.*10^6,PO2_rest,'Color',[1 0 0],'LineWidth',2)

function [CMRO2_stim] = extract_CMRO2(baseline_O2, stim_O2, PO2_increase, M)
for ii = 1:length(stim_O2)
    
    % M = 0.9;
    % PO2_increase = 5.7;
    
    PO2_artery = [40]; %mmHg
    %CMRO2 = 3/60; %M s^-1
    D = 2.8*10^-9; %m^2 s^-1
    alpha = 0.00139; %M mmHg^-1
    R_artery = 0.9*10^-5; %m
    Rt = 0.5*10^-4; %m
    r = Rt;
    PO2_rest = baseline_O2(ii);
    
    
    a = (r^2-R_artery^2)/(4*D*alpha);
    b = (Rt^2*log(r/R_artery))/(2*D*alpha);
    CMRO2_stim(1,ii) = (PO2_rest-PO2_artery)/(a - b);
    
    
    
    PO2_artery = PO2_artery + PO2_increase; %mmHg
    R_artery = 0.9*10^-5*M; %m
    PO2_rest = stim_O2(ii);
    
    a = (r^2-R_artery^2)/(4*D*alpha);
    b = (Rt^2*log(r/R_artery))/(2*D*alpha);
    CMRO2_stim(2,ii) = (PO2_rest-PO2_artery)/(a - b);
    
    
end
end
