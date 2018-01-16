% a more manual version of setup days: :)

%setup days script for model runs for FCB1tank data - time is local and
%relfects T0 of the diluiton series experiments

mergedpath = '/Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/data_picoeuks_copy/processed2b/grouped/merged/';  %paths changed for Swallowtail
groupedpath = '/Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/data_picoeuks_copy/processed2b/grouped/';
modelpath = '/Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/model_picoeuks/input2b/';

culture = '3'; % '6 %'3'
port=culture;

%filelist = dir([groupedpath 'bFCB_200*' culture '.mat']);
filelist = dir([groupedpath 'FCB1_2012*' culture '.mat']);
%
load /Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/all_light
%Light data from MVCO:
Eall_MV = [date_met-4/24 Solar];  %date - now in local time (MVCO data comes in UTC) (local, matlab), relative PAR
Eall_MV(Eall_MV(:,2)<0,2) = 0;  %no negatives allowed
Eall_MV(:,1)=Eall_MV(:,1);

%Light data from the radiometer:
Eall_R = [sg_time sg_light];  %date -(local, matlab), relative PAR
Eall_R(Eall_R(:,2)<0,2) = 0;  %no negatives allowed
Eall_R(:,1)=Eall_R(:,1);

load /Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/DS_datekey
%volbins = 2.^[-3:0.125:4];
volbins=2.^[-5:1/6:4.5]; %58 bins
%keyboard
filenum = 1;
filename = filelist(filenum).name;
disp(['loading...' filename])
eval(['load ' groupedpath filename])
eval(['load ' mergedpath filename(1:end-8) 'merged' filename(end-7:end)])

cellclass = 4;  %cell group to model (1 = syn, 4 & 5 = euks)
numfiles = length(filelist);
%load /Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/DS_dawn_hours %contains matlab startday and then the hour when the DS experiments started
load /Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/model_code/DS_dawn_start6.mat
%%

% For the first file and index:

daylist = floor(min(cellresults(:,1))):floor(max(cellresults(:,1)));
DS_ind=7; %start here for good DS experiments
dayind = find(daylist==DS_datekey(DS_ind,1));
day = daylist(dayind);

%%
dielstarthr = DS_dawn(find(DS_dawn(:,1) == day),2);
if isnan(dielstarthr), dielstarthr = 0; end;
if isempty(dielstarthr), dielstarthr = 0; end;
offsetdate = cellresults(:,1) - dielstarthr/24;
dielind = find(floor(offsetdate) == day);


%% after the first day:
DS_ind=DS_ind+1;

dayind = find(daylist==DS_datekey(DS_ind,1));
day = daylist(dayind);
dielstarthr = DS_dawn(DS_dawn(:,1) == day,2);

offsetdate = cellresults(:,1) - dielstarthr/24;  %reset so diel start = 0 hour of day
clear cellvol
dielind = find(floor(offsetdate) == day);
disp(['day: ' num2str(day) ' :' datestr(day) ' Dawn Hour: ' num2str(dielstarthr)])

%% in case need to change dawn start hour:
dielstarthr = 6;
DS_dawn(DS_ind,2)=dielstarthr;
save DS_dawn_start3 DS_dawn DS_datekey
offsetdate = cellresults(:,1) - dielstarthr/24;  %reset so diel start = 0 hour of day
clear cellvol
dielind = find(floor(offsetdate) == day);
disp(['day: ' num2str(day) ' :' datestr(day) ' Dawn Hour: ' num2str(dielstarthr)])

%% Once you've gotten to DS_ind 7:

filenum = 7;
filename = filelist(filenum).name;
disp(['loading...' filename])
eval(['load ' groupedpath filename])
eval(['load ' mergedpath filename(1:end-8) 'merged' filename(end-7:end)])

%% In case are in last day of file:

%NEXT WHILE LOOP ONLY FOR LAST DAY IN CURRENT FILE:
while day >= floor(max(cellresults(:,1)))-1 & filenum ~= numfiles & ~isempty(dielind),  %last day in file, keep partial day and open next file
    oldcellresults = cellresults(dielind(1):end,:); %keep last day (fix bug from dielind to dielind(1):end 12/11/04)
    oldmergedwithclass = mergedwithclass(dielind(1):end);  %keep last day
    oldbeadmatch = beadmatch(dielind(1):end,:);  %keep last day
    filenum = filenum + 1;
    filename = filelist(filenum).name;  %next file
    %         mergedpath = '/Users/kristenhunter-cevera/Documents/MATLAB/FCB1tank/data/processed2b/grouped/merged/';
    disp(['loading...' filename])
    eval(['load ' groupedpath filename])
    eval(['load ' mergedpath filename(1:end-8) 'merged' filename(end-7:end)])
    if (strcmp(filename(1:4), 'oc01') | strcmp(filename(1:4), 'oc17') | strcmp(filename(1:4), 'oc19')),
        cellresults(:,1) = cellresults(:,1) + 4/24;  %convert local to UTC
    end;
    cellresults = [oldcellresults; cellresults];  %add last day to start of next file
    mergedwithclass = [oldmergedwithclass mergedwithclass];
    beadmatch = [oldbeadmatch; beadmatch];
    if isnan(dielstarthr) , dielstarthr = 0; end;
    if isempty(dielstarthr), dielstarthr = 0; end;
    offsetdate = cellresults(:,1) - dielstarthr/24;  %recalculate all dates and indices
    dielind = find(floor(offsetdate) == day);
    dayind = 1; %restart at beginning of daylist
    clear old*
end;  %if day = floor(cellresults(end,1)),
%%

%        daylist = floor(min(offsetdate)):floor(max(offsetdate));
daylist = floor(min(cellresults(:,1))):floor(max(cellresults(:,1)));
disp(['Length: ' num2str(length(dielind))])

if ~isempty(dielind) & dielind(end) < size(offsetdate,1),
    dielind = [dielind; dielind(end)+1];
    dielhr = floor((offsetdate(dielind) - day)*24);
    dielind = dielind(find(dielhr <= 24)); %in case one is over 25?
    dielhr = floor((offsetdate(dielind) - day)*24);
else   %case where no data points for day or last date
    dielhr = NaN;
end;
%keyboard
if dielhr(1) == 1,  %repeat first hr 1 for hr 0
    dielind = [dielind(1); dielind]
    dielhr = [0; dielhr];
    disp('missing first hour')
end;
if dielhr(end) == 23,  %repeat hr 23 as hr 24
    dielind = [dielind; dielind(end)]
    dielhr = [dielhr; 24];
    disp('missing last hour')
end;

disp(['Current Length: ' num2str(length(dielind))])
%
if dielhr(1) == 0 & dielhr(end) == 24 & isempty(find(diff(dielhr) > 2)),  %full day ==> proceed, otherwise just go to next day
    
    cellresultsfordiel = cellresults(dielind,:);
    for count = 1:length(dielind),  %get all cell data
        temp = mergedwithclass{dielind(count)};
        ind = find(temp(:,end) == cellclass);
        cellvol{count} = cytosub_SSC2vol(temp(ind,5)./beadmatch(dielind(count),5));  %SSC bead units values converted to volume
    end;
    %
    temp = find(diff(dielhr) == 0);  %case where more than one sample in same hr (only at quick acq restart)
    while ~isempty(temp), %every other one matches
        temp = temp(1); %deal with first one first...loop to get subsequent
        cellvol{temp} = [cellvol{temp}; cellvol{temp+1}]; %add together cells measured in same hour
        cellvol = cellvol([1:temp, temp+2:end]); %omit second in same hour
        %                cellresultsfordiel(temp,:) = mean(cellresultsfordiel(temp:temp+1,:));
        %                cellresultsfordiel(temp,:) = [mean(cellresultsfordiel(temp:temp+1,1)) cellresultsfordiel(temp,2)+cellresultsfordiel(temp+1,2) mean(cellresultsfordiel(temp:temp+1,3))];  %add together acq time...
        cellresultsfordiel(temp,:) = [mean(cellresultsfordiel(temp:temp+1,1)) cellresultsfordiel(temp,2)+cellresultsfordiel(temp+1,2) cellresultsfordiel(temp,3)+cellresultsfordiel(temp+1,3)];  %add together acq time...new mode heidi 12/22/06
        cellresultsfordiel = cellresultsfordiel([1:temp, temp+2:end],:);
        dielhr = dielhr([1:temp, temp+2:end]);
        dielind = dielind([1:temp, temp+2:end]);
        temp = find(diff(dielhr) == 0);  %check again
    end;
    temp = find(diff(dielhr) == 2);  %case where one hour missing
    while ~isempty(temp),  %average across 1 skipped hour
        temp = temp(1);
        tempvol{1} = [cellvol{temp}; cellvol{temp+1}]; %add together cells measured in surrounding hours
        cellvol = [cellvol([1:temp]) tempvol cellvol([temp+1:end])];
        tempcellres = mean(cellresultsfordiel([temp, temp+1],:));
        tempcellres(2) = cellresultsfordiel(temp,2) + cellresultsfordiel(temp+1,2);  %add together acq times
        cellresultsfordiel = [cellresultsfordiel(1:temp,:); tempcellres; cellresultsfordiel(temp+1:end,:)];
        dielhr = [dielhr(1:temp); dielhr(temp)+1; dielhr(temp+1:end)];
        dielind = [dielind(1:temp); dielind(temp)+1; dielind(temp+1:end)];
        temp = find(diff(dielhr) == 2); %look for another...
        disp('  averaged over one missing hour')
    end;
    clear temp count
end
% Light:

if DS_ind <= 6
    Eall=Eall_R;
else Eall=Eall_MV;
end

if day==735167 %special case for day where light data is missing, but looks almost the
%same as the following day...
    
    Eall_MV2=Eall_MV;
    Eall_MV2(:,1)=Eall_MV2(:,1)-1;
    Eall2=Eall_MV2;
    
    Eind = find(floor(Eall(:,1) - dielstarthr/24) == day);
    Eind2 = find(Eall2(:,1) > Eall(Eind(end),1) & floor(Eall2(:,1) - dielstarthr/24) == day);  %
       
    if ~isempty(Eind),
        Eind2 = [Eind2; Eind2(end)+1];
        Edata = [Eall(Eind,:); Eall2(Eind2,:)];
        
        Edata(:,1) = (Edata(:,1) - day)*24 - dielstarthr;  %Radiation data for hours from 0 to 24
    else
        Edata = [NaN NaN];
    end;
    
else %no light data missing and proceed as normal:
    
    Eind = find(floor(Eall(:,1) - dielstarthr/24) == day);  %
    if ~isempty(Eind),
        Eind = [Eind; Eind(end)+1];
        Edata = Eall(Eind,:);
        if (Eall(Eind(end),1) - Eall(Eind(end-1))) > 2/24, %make sure end of day does not interpolate to middle of next light
            Edata(end,2) = 0;
        end;
        Edata(:,1) = (Edata(:,1) - day)*24 - dielstarthr;  %Radiation data for hours from 0 to 24
    else
        Edata = [NaN NaN];
    end;
end

clear Eind
% See how the distributions are looking:
%next from xiopt4
%
Vhists = zeros(length(volbins),length(cellvol));
N_dist =  Vhists;

for i = 1:length(cellvol) %calculate observed size distributions
    N_dist(:,i) = histc(cellvol{i}, volbins);  %size distribution
    Vhists(:,i) = N_dist(:,i)/sum(N_dist(:,i));    %normalized size distribution
end
%            cellsperml = [sum(N_dist)./cellresultsfordiel(:,2)'./cellresultsfordiel(:,3)'];
cellsperml = [sum(N_dist)./cellresultsfordiel(:,3)']; %heidi 12/22/06

figure(6), set(gcf,'Position',[ 124         555        1323         417])
subplot(131), pcolor(Vhists), title('Vhists','Fontsize',14)
subplot(132), pcolor(N_dist), title('N_dist','Fontsize',14)
subplot(133), plot(Edata(:,1),Edata(:,2),'.-','MarkerSize',8)
%%

eval(['save ' modelpath 'day' num2str(day) '_' num2str(port) '_picoeuk_data volbins Vhists N_dist cellsperml dielstarthr cellvol Edata'])  %save data for rerunning batch later




