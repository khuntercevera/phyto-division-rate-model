% if year2do <= 2015
%     eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99'',''ascii'');'])
% else
%     eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99B.txt'',''ascii'');'])
% end
eval(['X=load(''' envpath num2str(year2do) '_MetDat_s.C99'',''ascii'');'])
 
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
fprintf('Found out of order time for %1.0f instances!\n',length(t))
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

%% For 2010, need to use nantucket buoy data!
if year2do == 2010
    
    filename = '/Volumes/Lab_data/MVCO/FCB/MVCO_Jan2010/model_code/buoy_solar.txt';
    startRow = 3;
    formatSpec = '%4f%3f%3f%3f%3f%7f%7f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
    buoy_data=cell2mat(dataArray(:,1:8));  
    buoy_date=datenum([buoy_data(:,1:5) zeros(length(buoy_data(:,1)),1)]);
    clearvars filename startRow formatSpec fileID dataArray ans;
        
    %% merge mvco and nantucket buoy data:
    % go through each day and see where data is available - then splice together:
    merged_solar=[];
    merged_date=[];
    merged_hour=[];
    for yearday=1:365
        
        ind_mvco=find(floor(date_met)==(yearday+datenum('1-0-2010')));
        ind_buoy=find(floor(buoy_date)==(yearday+datenum('1-0-2010')));
        if (~isempty(ind_mvco) & length(ind_mvco) > 60) | (~isempty(ind_mvco) & yearday < 121) %close enough to a full day but still use mvco as doesn't have buoy data until 5/1
            merged_date=[merged_date; date_met(ind_mvco)];
            merged_solar=[merged_solar; Solar(ind_mvco)];
            merged_hour=[merged_hour; Hour(ind_mvco)];
        elseif ~isempty(ind_buoy)
            merged_date=[merged_date; buoy_date(ind_buoy)];
            merged_solar=[merged_solar; buoy_data(ind_buoy,7)];
            merged_hour=[merged_hour; buoy_data(ind_buoy,4)];
        end
        
    end
    
    date_met=merged_date;
    Solar=merged_solar;
    Hour=merged_hour;
    yd_met=date_met-datenum('1-0-10');
    
    
end


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

%load ~/Documents/phyto-division-rate-model/median_dawn.mat %from raw data, exected UTC values of dawn and dusk over 13 years :)
if ~isempty(strfind(computer,'WIN'))
    load \\sosiknas1\lab_data\mvco\FCB\Syn_divrate_model\median_dawn.mat 
else
    load /Volumes/Lab_Data/MVCO/FCB/Syn_divrate_model/median_dawn.mat  %from raw data, exected UTC values of dawn and dusk over 13 years :)
end

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
    disp(unqday(count))
    ind = find(floor(date_met - 5/24) == unqday(count));  %rough estimate for local day
    if ~isnan(dawnhr(count)),
        %
        ind2 = find(lhour(ind) >= dawnhr(count) & lhour(ind) <= duskhr(count)); %grab indexes that correspond to light hours!
        if ~isempty(ind2)
            %
            if count > 1 & count < length(unqday), % include pts just before and just after light to check for start and end gaps
                ind3 = [ind(ind2(1))-1; ind(ind2); ind(ind2(end))+1];  % include pts just before and just after light
            elseif count == 1,
                ind3 = [ind(ind2); ind(ind2(end))+1];
            else %last day
                ind3 = [ind(ind2(1))-1; ind(ind2)];
            end;
            
            %ind3 is current days just before and after dawn light hours :)
            %
            gaps = find(diff(yd_met(ind3)) > 3/24);  %find gaps > 2 h during light
            %
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
                gaps = find(diff(yd_met(ind(ind2))) > 3/24);
                if ~isempty(gaps)
                    dawnhr(count) = NaN;
                end;
            end;
        end;  %indexes for light hours in day
    end;   %indexes for local day
end;  %count


%% now check for nighttime noise"
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

%skip last local time hours of 12-31 for previous year
if year2do == 2003
    dawn = [unqday(2:end) dawnhr(2:end)];
else
    dawn = [unqday dawnhr];
end


if solarplotflag
    clf
    for j=1:366
        day=datenum(['1-0-' num2str(year2do)])+j;
        %color in night:
        f1=fill([day-1+dusk_median(j)/24; day+dawn_median(j)/24; day+dawn_median(j)/24; day-1+dusk_median(j)/24],[0 0 1000 1000],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    
    %plot the data!
    plot(date_met,Solar,'.-','color',[0.0265 0.6137 0.8135])
    set(gca, 'layer', 'top')
    
    %add dawn lines:
    for i=1:length(dawn)
        line([dawn(i,1)+dawn(i,2)/24 dawn(i,1)+dawn(i,2)/24],[0 1000],'linewidth',2,'color',[0.8 0.5 0]) %first matlab default color: [0 0.5 0.8]
    end
end
datetick

%%

eval(['save ' solarsavepath 'solar' num2str(year2do) '.mat Solar date_met dawn'])
