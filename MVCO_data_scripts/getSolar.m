if year2do == 2016
    eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99B.txt'',''ascii'');'])
elseif year2do == 2017
    eval(['X=load(''' envpath num2str(year2do) '_MetDat_sB.C99.txt'',''ascii'');'])
else
    eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99'',''ascii'');'])
end

yd_met=X(:,2);
Hour = X(:,5);
Solar=X(:,14);  %Solar_campmt_median
% Year = X(:,1);
% Day = X(:,4);

clear X

%fix for duplicated data, sometimes interspersed with good data?
[tt, it]=sort(yd_met);
yd_met=yd_met(it);
Solar=Solar(it);
Hour=Hour(it);

[tt2, iu]=unique(yd_met);
yd_met=yd_met(iu);
Solar=Solar(iu);
Hour=Hour(iu);

t = find(diff(yd_met) < 0);
fprintf('Found "out of order" time for %1.0f instances!\n',length(t))
yd_met(t+1) = NaN;
Solar(t+1) = NaN;
Hour(t+1) = NaN;
% Day(t+1) = NaN;
% Year(t+1) = NaN;

%exclude nan's:
ind = find(~isnan(Solar) & ~isnan(yd_met));
yd_met = yd_met(ind);
Solar = Solar(ind);
Hour = Hour(ind);

date_met = yd_met + datenum(['1-0-' num2str(year2do)]);  %should be UTC
dawnlevel = 5;   %threshold for light

    
%% find dawn and dusk of each day:
unqday = unique(floor(date_met - 5/24)); %approx number of local days
unqday = unqday(~isnan(unqday));
yrdy=find_yearday(unqday);
dawnhr=nan(366,1);
duskhr=nan(366,1);

for count = 1:length(unqday),
    ind = find(floor(date_met - 5/24) == unqday(count));
    ind2 = find(Solar(ind) > dawnlevel);
    if ~isempty(ind2),
        dawnhr(yrdy(count)) = Hour(ind(ind2(1))) - 1;  %1 h before dawnlevel?
        duskhr(yrdy(count)) = Hour(ind(ind2(end))) + 1;  %next hour after end of light
        if dawnhr(yrdy(count)) == -1, keyboard, end;
    else
        dawnhr(yrdy(count)) = NaN;
        duskhr(yrdy(count)) = NaN;
    end;
end;

% Need a good estimate of where dawn/dusk is to detect gaps:

%constructed from raw_data with get_average_dawn.m, 
%load expected UTC values of dawn and dusk over 13 years :)
if ~isempty(strfind(computer,'WIN'))
    load \\sosiknas1\lab_data\mvco\FCB\Syn_divrate_model\median_dawn.mat
else
    load /Users/kristenhunter-cevera/phyto-division-rate-model/setup_and_postprocessing_src/median_dawn.mat  
end

%replace dawn and dusk hours if deviate from median:
ind = find(abs(dawnhr-dawn_median) > 1);
dawnhr(ind) = dawn_median(ind);
dawnhr=dawnhr(yrdy);
ind = find(abs(duskhr-dusk_median) > 1);
duskhr(ind) = dusk_median(ind);
duskhr=duskhr(yrdy);

temp = datevec(date_met(~isnan(date_met)));
lhour = NaN*ones(length(date_met),1);
lhour(~isnan(date_met)) = temp(:,4); %UTC hour
clear temp

%% Find gaps or not good days and flag accordingly:

for count = 1:length(unqday),
    %
    disp(datestr(unqday(count)))
    ind = find(floor(date_met - 5/24) == unqday(count));  %rough estimate for local day
    
    if ~isnan(dawnhr(count)),
        %
        ind2 = find(lhour(ind) >= dawnhr(count) & lhour(ind) <= duskhr(count)); %grab indexes that correspond to light hours!
        if ~isempty(ind2)
            %
            if count > 1 & count < length(unqday), % include pts just before and just after light to check for start and end gaps
                ind3 = [ind(ind2(1))-1; ind(ind2); ind(ind2(end))+1];  % include pts just before and just after light
            elseif count == 1 %first day
                ind3 = [ind(ind2); ind(ind2(end))+1];
            else %last day
                ind3 = [ind(ind2(1))-1; ind(ind2)];
            end;
            
            %ind3 is current day just before and after dawn light hours :)
            gaps = find(diff(yd_met(ind3)) >= 2/24);  %find gaps > 2 h during light
            
            if ~isempty(gaps),
                %
                if gaps(end) == length(ind3) - 1,
                    if lhour(ind(end)) < duskhr(count) - 3, %big gap around dusk
                        dawnhr(count) = NaN;
                    else
                        date_met = [date_met(1:ind3(end-1)); unqday(count) + (duskhr(count)+.001)/24; date_met(ind3(end-1)+1:end)];
                        Solar = [Solar(1:ind3(end-1)); 0; Solar(ind3(end-1)+1:end)];
                        lhour  = [lhour(1:ind3(end-1)); duskhr(count)+.001; lhour(ind3(end-1)+1:end)];
                        yd_met  = [yd_met(1:ind3(end-1)); floor(yd_met(ind3(2))) + (duskhr(count)+4.001)/24; yd_met(ind3(end-1)+1:end)];
                    end;
                end;
                %
                if gaps(1) == 1 && ~isnan(dawnhr(count)) %second logic is for if previous dusk gap was detected and now this is flagged
                    if lhour(ind(1)) >= dawnhr(count) + 3,  %big gap before dawn
                        dawnhr(count) = NaN;
                    else
                        date_met = [date_met(1:ind3(1)); unqday(count) + (dawnhr(count)+.001)/24; date_met(ind3(2):end)];
                        Solar = [Solar(1:ind3(1)); 0; Solar(ind3(2):end)];
                        lhour  = [lhour(1:ind3(1)); dawnhr(count)+.001; lhour(ind3(2):end)];
                        yd_met  = [yd_met(1:ind3(1)); floor(yd_met(ind3(2))) + (dawnhr(count)+4.001)/24; yd_met(ind3(2):end)];
                    end;
                end;
                
                %With edits..now refind any gaps
                gaps = find(diff(yd_met(ind(ind2))) >= 2/24);
                if ~isempty(gaps)
                    dawnhr(count) = NaN;
                end;
            end;
        end;  %indexes for light hours in day
    end;   %indexes for local day
end;  %count


%% %%%%%%% identify missing days, and for years 2005-2007, 2010-2013, check if buoy data can be substituted! %%%%%%%%%%%%%%%%%%%%%%%%

% For most of 2010, need to use nantucket buoy data!
%Buoy data is formatted as:
%#YY  MM DD hh mm  SRAD1  SWRAD  LWRAD
%#yr  mo dy hr mn   w/m2   w/m2   w/m2
%but appears to use different instruments for different years...

if ismember(year2do,[2005:2007 2010:2013]) && buoy_flag==1; %meaning, yes - you'd like to splice in buoy data
    
    disp('Checking Nantucket Buoy data for any missing gaps in MVCO record...')
    filename = fullfile(envpath,'Nantucket_buoy44008_lightdata',['44008r' num2str(year2do)]);
    startRow = 3;
    formatSpec = '%4f%3f%3f%3f%3f%7f%7f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
    buoy_data=cell2mat(dataArray(:,1:8));
    buoy_date=datenum([buoy_data(:,1:5) zeros(length(buoy_data(:,1)),1)]);
    buoy_date=buoy_date-((2/3)/24); %seems to better align with MVCO data if shifted
    %data is reported for ending hour, but not sure when in the hour data was taken!
    clearvars filename startRow formatSpec fileID dataArray ans;  
    
    if ismember(year2do,2010:2013)
        sw_rad=buoy_data(:,7);       
    elseif ismember(year2do,2005:2007)
        sw_rad=buoy_data(:,6);
    end
          
    %% splice in nantucket buoy data:
    % go through each day and see when data is available - then splice together:
    missing_days=setxor(find_yearday(unqday(~isnan(dawnhr))),1:366);
    buoy_added=[];
    mvco_to_remove=[];
    duskhr2=[];
    dawnhr2=[];
    unqday2=[];
    
    for q=1:length(missing_days) %check if buoy data has them!
        tempday=missing_days(q)+datenum(['1-0-' num2str(year2do)]);     
        qq=find(buoy_date-5/24 >=tempday & buoy_date-5/24 <tempday+1); %local frame of reference!
        
        if ~isempty(qq) && length(qq) >=22 %buoy data is available and seems complete!            
            %add data from buoy            
            buoy_added=[buoy_added; buoy_date(qq) sw_rad(qq)];   
            
            %add to dawn:
            ind2 = find(sw_rad(qq) > dawnlevel);
            
            tempdawnhr=buoy_data(qq(ind2(1)),4)-1-1; %1 h before dawnlevel? (-1 is for offset of data recording)
            if abs(tempdawnhr-dawn_median) > 1, tempdawnhr = dawn_median(missing_days(q)); end
            tempduskhr=buoy_data(qq(ind2(end)),4)+1-1; %next hour after end of light plus offset
            if abs(tempduskhr-dusk_median) > 1, tempduskhr = dusk_median(missing_days(q)); end
            
            dawnhr2 = [dawnhr2; tempdawnhr];  
            duskhr2 = [duskhr2; tempduskhr];  
            unqday2 = [unqday2; tempday];
            
            %data to remove from MVCO track:
            mm=find(date_met >= buoy_date(qq(1)) & date_met <= buoy_date(qq(end)));
            mvco_to_remove=[mvco_to_remove; mm];
        end   
    end
    
    % now add in and remove data - yikes!  
    [ind]=setdiff(1:length(date_met),mvco_to_remove); %mvco data to keep, that won't be replaced
    date_met_temp=[date_met(ind); buoy_added(:,1)];
    solar_temp=[Solar(ind); buoy_added(:,2)];  
    [~, is]= sort(date_met_temp);
    date_met=date_met_temp(is);
    Solar=solar_temp(is);
    
    unqtemp=[unqday; unqday2];
    dawnhr_temp=[dawnhr; dawnhr2];
    duskhr_temp=[duskhr; duskhr2];
    [~,is2]=sort(unqtemp);
    unqday=unqtemp(is2);
    dawnhr=dawnhr_temp(is2);
    duskhr=duskhr_temp(is2);
    
    yrdy=find_yearday(unqday);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  check for nighttime noise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

date_met_local=date_met-5/24; %for easier handling

for j=1:length(unqday)
    
    day=unqday(j);
    ii=find(date_met_local >= day & date_met_local <= day+1);
    
    %find any light levels during nightime hours:
    nn=find(date_met_local(ii) <= day+dawn_median(yrdy(j))/24-5/24 | date_met_local(ii) >= day+dusk_median(yrdy(j))/24-5/24); %min normal dawn and max normal dusk
    
    %expected night:
    if ~isempty(nn)
        if any(Solar(ii(nn)) > 20)
            disp(['Found "light" readings during dark period for day: ' num2str(day) ': ' datestr(day) ' correcting...'])
            Solar(ii(nn)) = 0;
            
            %             plot(date_met_local(ii(nn)),Solar(ii(nn)),'.-')
            %             xlim([day-3 day+3])
            %             datetick('x','mm dd','keeplimits')
            %             pause
            
        end
    end
end

%% %%%%%% a bit more clean up.... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%skip last local time hours of 12-31 for previous year:
dawn = [unqday dawnhr];
[yy,~,~,~,~,~]=datevec(dawn(:,1));
jj=find(yy==year2do);
dawn=dawn(jj,:); %occassionally Dec 31st of previous year is included!


%% %%%%%%%%%%%%%%%%%%%%%%   AND PLOT!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if solarplotflag
    clf
    for j=1:366
        day=datenum(['1-0-' num2str(year2do)])+j;
        %color in night:
        f1=fill([day-1+dusk_median(j)/24; day+dawn_median(j)/24; day+dawn_median(j)/24; day-1+dusk_median(j)/24],[0 0 1000 1000],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    
    %plot the data!
    h1=plot(date_met,Solar,'.-','color',[0.0265 0.6137 0.8135]);
    set(gca, 'layer', 'top')
    if exist('buoy_added')
        hold on
        h2=plot(buoy_added(:,1),buoy_added(:,2),'o','linewidth',2,'markersize',4,'color',[0 0 0.7]);
        legend([h1(1); h2(1)],'MVCO data','Buoy data','location','NorthOutside')
        title([num2str(year2do) ' Solar data: ' num2str(size(unqday2,1)) ' days added from buoy data'])
    end
    
    %add dawn lines:
    for i=1:length(dawn)
        line([dawn(i,1)+dawn(i,2)/24 dawn(i,1)+dawn(i,2)/24],[0 1000],'linewidth',2,'color',[0.8 0.5 0]) %first matlab default color: [0 0.5 0.8]
    end
end

xlim([datenum(['1-1-' num2str(year2do)]) datenum(['12-31-' num2str(year2do)])])
ylim([-10 max(Solar)])
datetick('x','keeplimits')

%% %SAVE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('buoy_added')
    eval(['save ' solarsavepath 'solar' num2str(year2do) '.mat Solar date_met dawn buoy_added'])
else
    eval(['save ' solarsavepath 'solar' num2str(year2do) '.mat Solar date_met dawn'])
end
%% screen days, one by one:
%  for count = 1:length(unqday),
%     ind = find(floor(date_met) == unqday(count));
%     figure(1), clf
%     plot(yd_met(ind), Solar(ind), '.-')
%     hold on
%     line([dawnhr(count) dawnhr(count)]/24+unqday(count)-datenum(['1-0-' num2str(year2do)]), [0 1000], 'color', 'r')
%     line([duskhr(count) duskhr(count)]/24+unqday(count)-datenum(['1-0-' num2str(year2do)]), [0 1000], 'color', 'g')
%     title([num2str(unqday(count)) ' ; ' datestr(unqday(count))])
%     pause
%  end;
%