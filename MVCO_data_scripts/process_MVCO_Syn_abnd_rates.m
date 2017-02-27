%Summarize and process MVCO Synechococcus abundance, model division rate estimates, net growth and loss rate

% MAKE SURE CONNECTED TO SOSIKNAS FIRST!!!

addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools
addpath /Volumes/Lab_data/MVCO/FCB/Syn_divrate_model/

%% Synechococcus cell abundance, fluorescence and size:

allmatdate=[];
allsynconc=[];
allbeads=[];
allsynSSC=[];
allsynSSCmode=[];
allsynPE=[];
allsynPEmode=[];
allsynCHL=[];
allsynCHLmode=[];
allsynvol=[];
allsynvolmode=[];

for yearlabel=2003:2016
    
    switch yearlabel
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
    
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(yearlabel) '/data/processed/grouped/groupsum.mat'])
    
    %eval(['matdate_' num2str(yearlabel) '=cellresultsall(:,1);'])
    %eval(['synconc_' num2str(yearlabel) '=cellNUMall(:,1)./cellresultsall(:,3);'])
    to_use=exclude_data(cellresultsall,yearlabel);
    
%     clf,
%     subplot(1,2,1,'replace')
%     plot(cellresultsall(:,1),cellNUMall(:,1)./cellresultsall(:,3),'r.-')
%     hold on
%     plot(cellresultsall(to_use,1),cellNUMall(to_use,1)./cellresultsall(to_use,3),'.-','color',[0.2081    0.1663    0.5292])
%     title(num2str(yearlabel))
%     subplot(1,2,2,'replace')
%     plot(cellresultsall(:,1),cellNUMall(:,1)./cellresultsall(:,3),'r.-')
%     hold on
%     plot(cellresultsall(to_use,1),cellNUMall(to_use,1)./cellresultsall(to_use,3),'.-','color',[0.2081    0.1663    0.5292])
%     title(num2str(yearlabel)), set(gca,'yscale','log')
%     pause
%     
    allmatdate=[allmatdate; cellresultsall(to_use,1)];
    allsynconc=[allsynconc; cellNUMall(to_use,1)./cellresultsall(to_use,3)];
    
    allsynSSC=[allsynSSC; cellSSCall(to_use,1)./beadmatchall(to_use,5)];
    allsynvol = [allsynvol; cytosub_SSC2vol(cellSSCall(to_use,1)./beadmatchall(to_use,5))]; 
    allsynPE=[allsynPE; cellPEall(to_use,1)./beadmatchall(to_use,2)];
    allsynCHL=[allsynCHL; cellCHLall(to_use,1)./beadmatchall(to_use,4)];
    
    allsynSSCmode=[allsynSSCmode; cellSSCmodeall(to_use,1)./beadmatchall(to_use,5)];
    allsynvolmode = [allsynvolmode; cytosub_SSC2vol(cellSSCmodeall(to_use,1)./beadmatchall(to_use,5))]; 
    allsynPEmode=[allsynPEmode; cellPEmodeall(to_use,1)./beadmatchall(to_use,2)];
    allsynCHLmode=[allsynCHLmode; cellCHLmodeall(to_use,1)./beadmatchall(to_use,4)];
    %allbeads=[allbeads; beadmatchall(to_use,[2 5])];
    
end

%%
% remove nan's
ii=find(~isnan(allmatdate));
allmatdate=allmatdate(ii); allsynconc=allsynconc(ii);  

allsynPE=allsynPE(ii); %allbeads=allbeads(ii,:);
allsynPEmode=allsynPEmode(ii);  
allsynCHL=allsynCHL(ii);
allsynCHLmode=allsynCHLmode(ii);
allsynSSC=allsynSSC(ii);
allsynSSCmode=allsynSSCmode(ii);
allsynvol=allsynvol(ii);
allsynvolmode=allsynvolmode(ii);
%allmatdate=allmatdate-4/24; %shift from UTC time to local time

%smooth the abundance data over 48 hours, with data separated by 1 day
%treated as separate chunks:
synrunavg=mvco_running_average(allmatdate, allsynconc,48,1);

clearvars -except synrunavg allmatdate allsynconc allsynSSC* allsynPE* allsynCHL* allsynvol* %allbeads

%% bin syn cell abundance:

[time_syn_ns_dy, daily_syn_ns] = timeseries2ydmat(allmatdate, allsynconc); %raw, unsmoothed abundance
[time_syn, daily_syn, synyears, ydmu] = timeseries2ydmat(allmatdate, synrunavg); %smoothed abundance
syn_avg=nanmean(daily_syn,2);
syn_avg(end)=NaN; %this is because only one year with a leap year could be included, so the average isn't really good....
syn_std=nanstd(daily_syn,0,2);

[time_PE, daily_PE] = timeseries2ydmat(allmatdate, allsynPE); %Syn PE fluorescence
[time_CHL, daily_CHL] = timeseries2ydmat(allmatdate, allsynCHL); %Syn CHL fluorescence
[time_SSC, daily_SSC] = timeseries2ydmat(allmatdate, allsynSSC); %SSC, bead normalized
[time_vol, daily_vol] = timeseries2ydmat(allmatdate, allsynvol); %Cell volume from SSC-bead normalized

PE_avg=nanmean(daily_PE,2);
CHL_avg=nanmean(daily_CHL,2);
SSC_avg=nanmean(daily_SSC,2);
vol_avg=nanmean(daily_vol,2);
[PE_avg_wk, PE_std_wk] = dy2wkmn_climatology(daily_PE, synyears);
[CHL_avg_wk, CHL_std_wk] = dy2wkmn_climatology(daily_CHL, synyears);
[SSC_avg_wk, SSC_std_wk] = dy2wkmn_climatology(daily_SSC, synyears);
[vol_avg_wk, vol_std_wk] = dy2wkmn_climatology(daily_vol, synyears);

%for the modes:
[time_PEmode, daily_PEmode] = timeseries2ydmat(allmatdate, allsynPEmode); %Syn PE fluorescence
[time_CHLmode, daily_CHLmode] = timeseries2ydmat(allmatdate, allsynCHLmode); %Syn CHL fluorescence
[time_SSCmode, daily_SSCmode] = timeseries2ydmat(allmatdate, allsynSSCmode); %smoothed abundance
[time_volmode, daily_volmode] = timeseries2ydmat(allmatdate, allsynvolmode); %smoothed abundance

PE_mode_avg=nanmean(daily_PEmode,2);
CHL_mode_avg=nanmean(daily_CHLmode,2);
SSC_mode_avg=nanmean(daily_SSCmode,2);
vol_mode_avg=nanmean(daily_volmode,2);
[PE_mode_avg_wk, PE_mode_std_wk] = dy2wkmn_climatology(daily_PEmode, synyears);
[CHL_mode_avg_wk, CHL_mode_std_wk] = dy2wkmn_climatology(daily_CHLmode, synyears);
[SSC_mode_avg_wk, SSC_mode_std_wk] = dy2wkmn_climatology(daily_SSCmode, synyears);
[vol_mode_avg_wk, vol_mode_std_wk] = dy2wkmn_climatology(daily_volmode, synyears);

[weekly_syn, time_syn_wk, yd_wk] = ydmat2weeklymat(daily_syn, synyears);
[syn_avg_wk, syn_std_wk] = dy2wkmn_climatology(daily_syn, synyears);
% syn_avg_wk=nanmean(weekly_syn,2);
% syn_std_wk=nanstd(weekly_syn,0,2);


%% Division rates from the model:
%-----------------------------------------------------------------------------------------------------------------------------------------------

allgrowthrates=[];
modelallmatdate=[];
allMR=[];

%pathname='/Users/kristenhunter-cevera/Documents/MATLAB/Synechococcus/Fraction_Distribution_Model/Dirichlet_MN_model/Partial_hour_Dirichlet_MN_plateaugamma/';
%filelist=dir([pathname 'mvco_13par_dmn_*_beadmean_Jan2015.mat']); %for
%13param model

% pathname='/Users/kristenhunter-cevera/Documents/MATLAB/Syn_division_rate_models/Fraction_Distribution_Model/14_param_model/';
% filelist=dir([pathname 'mvco_14par_dmn_*_beadmean_Jan2015.mat']);
%
% for j=1:length(filelist)
%
%     filename=filelist(j).name;
%     eval(['load ' pathname filename])
%     year=filename(16:19);
%     eval(['MR=modelresults' year '_np;'])
%
%     % allgrowthrates=[allgrowthrates; MR(:,16)]; %13param_model
%     allgrowthrates=[allgrowthrates; MR(:,17)];
%     modelallmatdate=[modelallmatdate; MR(:,1)];
%     allMR=[allMR; MR];
%
% end

for yearlabel=2003:2016
    
    switch yearlabel
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
    
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(yearlabel) '/model/output_July2016/mvco_14par_dmn_' num2str(yearlabel) '.mat'])
    
    [days2redo, days2exclude]=exclude_modeldata(yearlabel);
    
    if ~isempty(days2exclude)
        days2exclude=str2num(cell2mat(days2exclude(:,1)));
        to_use=find(ismember(modelresults(:,1),days2exclude)==0);
    else
        to_use=find(~isnan(modelresults(:,1)) & modelresults(:,1)~=0); %just in case ;)
    end
    
%     clf,
%     plot(modelresults(:,1),modelresults(:,17),'r.','markersize',20)
%     hold on
%     plot(modelresults(to_use,1),modelresults(to_use,17),'.','color',[0.2081    0.1663    0.5292],'markersize',20)
%     title(num2str(yearlabel))
%     pause
    
    allgrowthrates=[allgrowthrates; modelresults(to_use,17)];
    modelallmatdate=[modelallmatdate; modelresults(to_use,1)];
    allMR=[allMR; modelresults];
    
end

%
% remove days that don't have a rate (model run didn't work...)
% ind=find(modelallmatdate==0);
% ii=setdiff(1:length(modelallmatdate),ind);
% modelallmatdate=modelallmatdate(ii);
% allgrowthrates=allgrowthrates(ii);
% allMR=allMR(ii,:);

% and those above 2 per day?
ind=find(allgrowthrates > 2);
allgrowthrates(ind)=NaN;

% % exclude spurious mu's based on temperature data:
% [mdate_mu, daily_mu, yearlist, ydmu ] = timeseries2ydmat(modelallmatdate, allgrowthrates);
% daily_mu=[nan(366,1) daily_mu]; %just if don't have 2003 data yet...
% jj=find((daily_mu > 0.5 & Tday < 8) | (daily_mu > 0.4 & Tday < 4)) % | (daily_mu > 0.3 & Tday < 6));
% daily_mu(jj)=nan;

clearvars 'modelresults*' 'allmodelruns*' ii ind MR filelist filename %remove individual years from workspace
%clearvars -except synrunavg allmatdate allsynconc allgrowthrates allMR modelallmatdate

% bin division rates
[time_mu, daily_mu, muyears] = timeseries2ydmat(modelallmatdate, allgrowthrates);

% %exclude suspicious mu's baesd on temperature?
[smu, new_mu_est, no_new_est]=replace_suspicious_mus(time_mu,daily_mu,0); %1 or 0 for plotflag!

%%
load('mvco_envdata_16Aug2016.mat', 'Tbeam_corr')
figure, plot(Tbeam_corr,daily_mu,'.','color',[0.4 0.4 0.4])
hold on, plot(Tbeam_corr(new_mu_est(:,1)),daily_mu(new_mu_est(:,1)),'rp')
plot(Tbeam_corr(no_new_est(:,1)),daily_mu(no_new_est(:,1)),'bp')
% plot(Tbeam_corr(new_mu_est(:,1)),new_mu_est(:,2),'kp')
plot(Tbeam_corr(no_new_est(:,1)),no_new_est(:,3),'cp')
%%
daily_mu(new_mu_est(:,1))=nan; %new_mu_est(:,2);
daily_mu(no_new_est(:,1))=nan;
%%
% adjust size for mu matrix that may not have all the years yet:
[yearsmissing, iy] = setdiff(synyears,muyears);
if ~isempty(yearsmissing)
    daily_mu(:,iy)=nan(size(daily_mu,1),length(iy));
end

mu_avg=nanmean(daily_mu,2);
mu_std=nanstd(daily_mu,0,2);

% The script saprse_weeklybin.m  bins the data according to calendar weeks. If however, a
%week only has 1 observation, then that obs is tacked onto the
%previous or following week, depending on which has closest data.


%[weekly_mu, time_mu_wk] = ydmat2weeklymat(daily_mu, synyears);
[time_mu_wk, weekly_mu]=sparse_weeklybin(time_mu,daily_mu,synyears);
[mu_avg_wk, mu_std_wk] = dy2wkmn_climatology(daily_mu, synyears);
% mu_avg_wk=nanmean(weekly_mu,2);
% mu_std_wk=nanstd(weekly_mu,0,2);


%% Compute net growth rate and loss rates:
%-----------------------------------------------------------------------------------------------------------------------------------------------
%find the syn abundance for each day that matches around dawn and use this to calculate net growth rate

%Find dawn hour from solar.mat files for each year:

dawnhours=[];
for year=2003:2016;
    
    switch year
        case 2003
            yearlabel='May';
        case 2004
            yearlabel='Apr';
        case 2005
            yearlabel='Apr';
        case 2006
            yearlabel='May';
        case 2007
            yearlabel='Mar';
        otherwise
            yearlabel='Jan';
    end
    
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' yearlabel num2str(year) '/model/solar' num2str(year) '.mat;']) %Change path here!
    dawnhours=[dawnhours; dawn];
    
end

[~, daily_dawn, ~, ~ ] = timeseries2ydmat(dawnhours(:,1), dawnhours(:,2));
dawn_avg=nanmedian(daily_dawn,2);

%Can also find it from the dielstarthr for each day in the setupdays files... see original climatology_data.m around June 2016

%% ----------------------------------------------------------------------------------------------------------
%For each day in the syn abundance data, see if there is a corresponding model day,
%record index for matching day and find the time period that best matches the
%diel period to use to calculate net growth rate:

daylist=unique(floor(allmatdate));
dawnavg=zeros(length(daylist),5);
dawn_ind=cell(length(daylist),1);
netind=[];

for j=1:length(daylist)
    
    day=daylist(j);
    
    %average over just a few hours around dawn if possible:
    dw=find(dawnhours(:,1)==day);
    dielst=dawnhours(dw,2);
    
    if isempty(dw) | isnan(dielst) %have not found a dawn value or has a nan value
        %use closest day value:
        dii=find(dawnhours(:,1) >= day-2 & dawnhours(:,1) <= day+2);
        if ~isempty(dii)
            dielst=nanmean(dawnhours(dii,2));
        else %if still empty, use climatology average:
            yd=find_yearday(day);
            dielst=dawn_avg(yd);
        end
        
    end
    
    %Now that we have an estimate of dawn for each day - Are there measurements near dawn?
    ind=find(allmatdate(:,1) >= day+((dielst(1)-2)/24) & allmatdate(:,1) <= day+((dielst(1)+2)/24));
    
    if ~isempty(ind)
        dawnavg(j,:)=[day nanmean(synrunavg(ind)) length(ind) 24*(allmatdate(ind(end))-allmatdate(ind(1))) dielst(1)];
        dawn_ind(j)={ind};
    else %if not, then need to skip...
        dawnavg(j,:)=[day nan(1,3) dielst];
    end
    
    
end

%% ---------------------------------------------------------------------------------------------------------
% Calculate the net growth rate
%there are a few ways to do this, temporarily have settled on taking the
%average cell conc around dawn (+/- 3 hours) and using this as the syn abnd
%for each day to do the net growth rate calc:

netmu_avg=[dawnavg(1:end-1,1) log(dawnavg(2:end,2)./dawnavg(1:end-1,2)) diff(dawnavg(:,1)) diff(dawnavg(:,5))]; %[day    net mu    diff in time     diff in dawnhr]

%exclude any gap longer than a day:
jj=find(netmu_avg(:,3) > 1);
netmu_avg(jj,2)=nan;

%Find the days that having matching model runs:
[tt, id, im]=intersect(netmu_avg(:,1),allMR(:,1));

%THERE ARE SOME MODEL DAYS THAT DO NOT HAVE A NET GROWTH RATE -
%LIKEY, THESE ARE DAYS WHERE DAWN HOURS CAN BE MISSING!!! 
%BUT NEED TO CHECK THIS!

%INDEED, all days come back as either missing 1-5 of first dawn hours, or
%abundance was screened out later due to other problems...

% Pause here for sanity figures!
% lightblue=[191 	239 	255]./[255 255 255];
% 
% figure
% plot(allmatdate,synrunavg,'k.-')
% hold on
% for j=1:length(dawnavg)
%     plot(allmatdate(dawn_ind{j}),synrunavg(dawn_ind{j}),'r.--','markersize',10)
% end
% 
% temp=netmu_avg(id,:);
% qq=find(isnan(temp(:,2)));
% missing_days=temp(qq,1);
% 
% for q=1:length(missing_days);
%     
%     day=missing_days(q);
%     w1=find(dawnavg(:,1)==day);
%     dielst=dawnavg(w1,5);
%     
%     w2=find(netmu_avg(:,1)==day);
%     
%     %How many hours of data are past dawn?
%     %  ind=find(allmatdate(:,1) > day+((dielst-2)/24) & allmatdate(:,1) < day+((dielst+3)/24)+1);
%     %  fprintf('Length of ind: %f\n',length(ind))
%     
%     f1=fill([day; day+1; day+1; day]+dielst/24,[0 0 5e5 5e5],lightblue);hold on;
%     set(f1,'linestyle', 'none'), xlim([day-3 day+3]),  uistack(f1,'bottom')
%     ylim([0.5*dawnavg(w1+1,2) 1.2*dawnavg(w1+1,2)])
%     
%     datetick('x','mm/dd/yy','keeplimits')
%     disp([num2str(q) ' out of ' num2str(length(missing_days))])
%     pause
%     
% end

%% Older version (pre June2016):
%netmu_avg=[dawnavg(1:end-1,1) log(dawnavg(2:end,2)./dawnavg(1:end-1,2))./(diff(dawnavg(:,1))) diff(dawnavg(:,1))];
%for those days that start a large gap - treat at NaNs:
%[gap jj]=sort(netmu_avg(:,3),'descend');
%ii=gap > 5; %greater than 5 days don't use;
%netmu_avg(jj(ii), 2)=nan;
%daily_mu=[nan(366,1) daily_mu]


%bin net growth rate data:

[time_net, daily_net, netyears] = timeseries2ydmat(netmu_avg(:,1), netmu_avg(:,2));
net_avg=nanmean(daily_net,2);
net_std=nanstd(daily_net,0,2);


[weekly_net, time_net_wk] = ydmat2weeklymat(daily_net, netyears );
[net_avg_wk, net_std_wk] = dy2wkmn_climatology(daily_net, netyears);
% net_avg_wk=nanmean(weekly_net,2);
% net_std_wk=nanstd(weekly_net,0,2);


%% Back calculate loss rates:
%---------------------------------------------------------------------------------------------------------

daily_loss=daily_mu-daily_net;
loss_avg=nanmean(daily_loss,2);
loss_std=nanstd(daily_loss,0,2);

[time_loss_wk, weekly_loss]=sparse_weeklybin(time_net,daily_loss,synyears);
[loss_avg_wk, loss_std_wk] = dy2wkmn_climatology(daily_loss, synyears);
% loss_avg_wk=nanmean(weekly_loss,2);
% loss_std_wk=nanstd(weekly_loss,0,2);

%%
dd=date;
eval(['save /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/syndata_' dd(1:2) dd(4:6) dd(8:end) '.mat *syn* *mu* *net* *loss* *PE* *SSC* *vol* *CHL* allgrowthrates allMR modelallmatdate allmatdate'])
