%Summarize and process MVCO Enivronmental Variables: Temperature, Light and Nutrients

% MAKE SURE CONNECTED TO SOSIKNAS FIRST!!!

addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools

%% Temperature:
%--------------------------------------------------------------------------------------------------------
%if have gaps in the shore-mast temperature - pad with data from the beam or perfrom a correlation correction:

load /Volumes/Lab_data/MVCO/EnvironmentalData/Tday_beam.mat %from the MVCO tower beam
eval('beam_years=yearlist;')
eval('beam_date=mdate_mat;')
load /Volumes/Lab_data/MVCO/EnvironmentalData/Tall_day.mat %from undersea node
eval('node_date=mdate;')
eval('node_years=year;')
eval('Tday_node=Tday;')

% If just want to pad missing values:
% Tday_padded=Tday;
% for j=4:11
%     ii=find(isnan(Tday(:,j)));
%     Tday_padded(ii,j)=Tday_beam(ii,j-3);
% end
%
% Tday=Tday_padded;

% Node-correlated temperature:

months=[1 30; 31 60; 61 90; 91 120; 121 150; 151 180; 181 210; 211 240; 241 270; 271 300; 301 330; 331 366];
T_corr=zeros(12,4);
%find the years where they overlap:
[temp, node_ii, beam_ii]=intersect(node_years,beam_years);
plot_flag=0;

for q=1:length(months)
    
    x=Tday_node(months(q,1):months(q,2),node_ii); x=x(:);
    X=[ones(length(x),1) x]; %add coefficent to do regression
    y=Tday_beam(months(q,1):months(q,2),beam_ii); y=y(:);
    [blin,~,~,~,stats] = regress(y,X);
    
    if plot_flag==1
        clf
        plot(x,y,'b.','markersize',10)
        hold on
        plot(0:25, blin(1)+blin(2)*(0:25),'g')
        title(num2str(months(q,1)),'fontsize',14), xlabel('Node Temp'), ylabel('Beam Temp')
        disp(stats([1 3]))
        
        pause
    end
    
    T_corr(q,:)=[blin' stats([1 3])];
    
end

% now correct with the adjusted node temperature

%add columns for missing years of Tbeam:
[temp, inode]=setdiff(node_years,beam_years); %years that are not included in the beam data...
ibeam=setxor(inode,1:16);
%%
Tbeam_corr=nan(size(Tday_node)); %this is the larger of the two...
Tbeam_corr(:,ibeam)=Tday_beam; %add with columns of beam temp
corr_years=node_years;

for q=1:12
    temp_beam=Tbeam_corr(months(q,1):months(q,2),:);
    ii=find(isnan(temp_beam));
    
    Tday_slice=Tday(months(q,1):months(q,2),:); %matches months and years of current Tday_beam
    temp_beam(ii)=T_corr(q,1)+T_corr(q,2)*Tday_slice(ii);
    
    Tbeam_corr(months(q,1):months(q,2),:)=temp_beam;
end
%
Tbeam_corr=Tbeam_corr(:,3:end); %just as only have starting 2003 data and syn data only goes to 2015!
corr_years=corr_years(3:end);
time_Tcorr=node_date(:,3:end);

% Bin and average temperature:
Tcorr_avg=nanmean(Tbeam_corr,2);
Tcorr_std=nanstd(Tbeam_corr,0,2);

[Tcorr_wk, time_Tcorr_wk, yd_wk]=ydmat2weeklymat(Tbeam_corr,corr_years);
[Tcorr_month, time_Tcorr_mn] = ydmat2monthlymat(Tbeam_corr,corr_years);

[Tcorr_avg_wk,Tcorr_std_wk, ~, Tcorr_avg_mn, Tcorr_avg_std, yd_mn] = dy2wkmn_climatology(Tbeam_corr,corr_years);
% Tcorr_avg_wk=nanmean(Tcorr_wk,2);
% Tcorr_std_wk=nanstd(Tcorr_wk,0,2);



%% Light Data
%--------------------------------------------------------------------------------------------------------

%Find integrated radiation and max radiation for each day:

total_light=[];

for year=2003:2016
    
    switch year
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    %
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year) '/model/solar' num2str(year) '.mat'])
    Solar(Solar <0)=0;
    date_met_local=date_met-4/24; %to change from UTC to local time for dealing with light :)
    
    daylist=unique(floor(date_met_local));
    daylist=daylist(~isnan(daylist));
    light=nan(length(daylist),5);
    
    for j=1:length(daylist)
        
        day=daylist(j);
        jj=find(dawn(:,1)==day);
        
        if ~isnan(dawn(jj,2)) %meaning that this data has already been screen - a dawn hour = good day
        
            ii=find(date_met_local >= day & date_met_local <= day+1);
        
%         %if there is a gap in the data, and this gap is more than 3 hours
%         %and occurs during the day, skip integration:
%         [mgap, im]=max(diff(date_met_local(ii)));
%         gaphour=24*(date_met_local(ii(im+1))-day);
%         
%         if any(diff(date_met_local(ii)) < 0) %meaning there is a time error in the data...
%             
%             %remove that datapoint for now...
%             %             jj=find(diff(date_met_local(ii)) < 0);
%             %             nii=setdiff(ii,ii(jj+1)); %exclude index that is out of sync
%             %             ii=nii;
%             fprintf('found day with time out of order: %6.2f! Excluding that day for now...\n',day)
%             
%             gap_flag=2;
%             light(j,:)=[day nan nan nan gap_flag nan];
%             
%         elseif  isempty(mgap) %for days that have only one time point...
%             
%             gap_flag=3;
%             light(j,:)=[day nan nan nan gap_flag nan];
%             
%         elseif (mgap > 3/24 && (gaphour < 5 || gaphour > 20)) || (mgap < 3/24) %if small gap or if gap appears during dark hours, go ahead and integrate:
%             
            gap_flag=0;
            radi=trapz(date_met_local(ii)',Solar(ii)',2);
            %Solar units are in W/m^2; if integrated over time this becomes:
            % W/m2 = J/s/m2 = ((kg m2)/s3) / m2 integrated over day time unit =
            % = ((kg m2) / s3) / m2 * d * (24hr/d)(60min/1hr)(60s/1min)= ((kg m2) / s3)*d*(86400s/d)/ m2 = ((kg m2)/s2) / m2 = J/m2
            %multiply by 1e-6 to get megajoules!
            
            light(j,:)=[day max(Solar(ii)) 86400*1e-6*radi length(ii) gap_flag]; %the units for integrated light are MJ/m2  ( -> the alt unit is w h m2 (must mulltiply by 10^6/3600) )
            
        else
            gap_flag=1;
            light(j,:)=[day nan nan nan gap_flag];
            %light(j,:)=[day nan nan length(ii) gap_flag mgap];
        end
        clear ii radi gap_flag
    end
    
%         %spot check for intergration!
%         clf     
%         %plot the data!
%         plot(date_met,Solar,'.-','color',[0.0265 0.6137 0.8135])
%         set(gca, 'layer', 'top')
%         hold on
%         plot(light(:,1)+12/24,10*light(:,3),'rp','markerfacecolor','r','markersize',18)    
%         datetick
%         pause
        
    total_light=[total_light;light];
    clear Solar date_met date_met_local light
    
end  %year

total_light=total_light(2:end,:); %removes a 2002 Dec 31st date from top!

% bin light data:

[time_light, light_int, lightyears] = timeseries2ydmat(total_light(:,1), total_light(:,3));
light_avg=nanmean(light_int,2);
light_std=nanstd(light_int,0,2);

[weekly_light, time_light_wk] = ydmat2weeklymat(light_int, lightyears);
[light_month, time_light_mn] = ydmat2monthlymat(light_int, lightyears );

[light_avg_wk,light_std_wk,~,light_avg_mn,light_std_mn,~] = dy2wkmn_climatology(light_int, lightyears);
% light_avg_wk=nanmean(weekly_light,2);
% light_std_wk=nanstd(weekly_light,0,2);



%% Nutrients!
%--------------------------------------------------------------------------------------------------------

load /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/MVCO_environmental_data/nut_data_repsB.mat

% find tower or node stations:

tower_ind=find(station4nut==4 | station4nut==5); %4 is tower, 5 is node I believe....
tower_nut=MVCO_nut_reps(tower_ind,:); %matrix subset
tower_nutdate=datenum(cell2mat(tower_nut(:,3)));

depth=cell2mat(tower_nut(:,7));

%only look at surface measurements for now:
surfind=find(depth <= 6);
surf_tower_nutdate=tower_nutdate(surfind); %more matrix subsets based on depth
surf_tower_nut=tower_nut(surfind,:);

% there may be some days were there were mulitple casts (i.e. SC trips) -
%average the nutrients that came from that day as well as across replicates:
surf_days=unique(surf_tower_nutdate);
NO3_surf=nan(length(surf_days),1);
NH4_surf=nan(length(surf_days),1);
PO4_surf=nan(length(surf_days),1);
SiOH_surf=nan(length(surf_days),1);

%remove any negative values for the nutrient measurements:
temp_NO3=cell2mat(surf_tower_nut(:,8:10));
temp_NO3(temp_NO3 < 0)=0;
temp_NH4=cell2mat(surf_tower_nut(:,11:13));
temp_NH4(temp_NH4 < 0)=0;
ii=find(temp_NH4 > 10); %too high to be real...
temp_NH4(ii)=NaN;
temp_PO4=cell2mat(surf_tower_nut(:,17:19));
temp_PO4(temp_PO4 < 0)=0;
temp_SiOH=cell2mat(surf_tower_nut(:,14:16));
temp_SiOH(temp_SiOH < 0)=0;

for j=1:length(surf_days)
    
    day=surf_days(j);
    ii=find(surf_tower_nutdate==day); %go thru and see if more than one cast for this depth or did tower and node stations same day:
    
    NO3_surf(j)=nanmean(nanmean(temp_NO3(ii,:)));
    NH4_surf(j)=nanmean(nanmean(temp_NH4(ii,:)));
    PO4_surf(j)=nanmean(nanmean(temp_PO4(ii,:)));
    SiOH_surf(j)=nanmean(nanmean(temp_SiOH(ii,:)));
    
    if length(ii) > 1
        disp(['date: ' datestr(day) ' : ' num2str(length(ii))])
    end
    
end


% bin nutrients:

[time_NO3, dy_NO3, nutyears, nutyrdy] = timeseries2ydmat(surf_days, NO3_surf);

[NO3_wk, time_NO3_wk] = ydmat2weeklymat(dy_NO3, nutyears );
[NO3_month, time_NO3_mn] = ydmat2monthlymat(dy_NO3, nutyears );
[~, ~, ~, NO3_avg_mn, NO3_std_mn, yd_mn] = dy2wkmn_climatology(dy_NO3, nutyears);
% NO3_avg_mn=nanmean(NO3_month,2);
% NO3_std_mn=nanstd(NO3_month,0,2);

[time_NH4, dy_NH4] = timeseries2ydmat(surf_days, NH4_surf);
[NH4_wk, time_NH4_wk] = ydmat2weeklymat(dy_NH4, nutyears );
[NH4_month, time_NH4_mn] = ydmat2monthlymat(dy_NH4, nutyears );
[~, ~, ~, NH4_avg_mn, NH4_std_mn, ~] = dy2wkmn_climatology(dy_NH4, nutyears);
% NH4_avg_mn=nanmean(NH4_month,2);
% NH4_std_mn=nanstd(NH4_month,0,2);

[time_PO4, dy_PO4] = timeseries2ydmat(surf_days, PO4_surf);
[PO4_wk, time_PO4_wk] = ydmat2weeklymat(dy_PO4, nutyears );
[PO4_month, time_PO4_mn] = ydmat2monthlymat(dy_PO4, nutyears );
[~, ~, ~, PO4_avg_mn, PO4_std_mn, ~] = dy2wkmn_climatology(dy_PO4, nutyears);
% PO4_avg_mn=nanmean(PO4_month,2);
% PO4_std_mn=nanstd(PO4_month,0,2);

[time_SiOH, dy_SiOH] = timeseries2ydmat(surf_days, SiOH_surf);
[SiOH_wk, time_SiOH_wk] = ydmat2weeklymat(dy_SiOH, nutyears );
[SiOH_month, time_SiOH_mn] = ydmat2monthlymat(dy_SiOH, nutyears );
[~, ~, ~, SiOH_avg_mn, SiOH_std_mn, ~] = dy2wkmn_climatology(dy_SiOH, nutyears);
% SiOH_avg_mn=nanmean(SiOH_month,2);
% SiOH_std_mn=nanstd(SiOH_month,0,2);

%Choose just one nutrient time:
if isequal(time_NO3,time_NH4,time_PO4,time_SiOH)
    time_nut=time_NH4;
end

if isequal(time_NO3_wk,time_NH4_wk,time_PO4_wk,time_SiOH_wk)
    time_nut_wk=time_NH4_wk;
end

if isequal(time_NO3_mn,time_NH4_mn,time_PO4_mn,time_SiOH_mn)
    time_nut_mn=time_NH4_mn;
end
% Need to sync time?
% synflag=1; %sync to max syn data available
%
% if synflag==1
%     synyears=2003:2015;
%     %eventually should get fancier and load from most recent processed
%     %syndata...
% end

dd=date;
%
eval(['save /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/mvco_envdata_' dd(1:2) dd(4:6) dd(8:end) '.mat *light* *Tcorr* corr_years Tbeam_corr dy_*  *_month nutyears nutyrdy time_nut* SiOH_wk NH4_wk NO3_wk PO4_wk SiOH_*_mn NH4_*_mn NO3_*_mn PO4_*_mn'])
%-regexp (?<!time)_\w*_mn