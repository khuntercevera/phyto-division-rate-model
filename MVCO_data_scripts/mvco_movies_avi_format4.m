function [  ] = mvco_movies_avi_format4( year , groupedpath0, mergedpath0, savepath )

%To have MATALB make movies without exporting figures to screen, must set
%figure 'Visible' property to 'off'. In the past, it seems that getframe
%didn't play nice with this, but seems to be workign now?

%Getframe or a combination of im2frame with printing the image take a
%loooong time. For the current figure generated below it's about ~10s per
%frame :( Basicially, turning a figure into an image takes awhile and there
%doesn't seem to be a good work around for speed for this. 


%must use avifile - the trick is to turn off visible in the figure and then
%do not call getframe, which seems it must export to screen in order to
%function, but rather use this undocumented hardcopy function zbuffer_cdata
%all within the MATLAB VideoWriter functions

% IF RUNNING FROM fcbmain_field1.m , these parts aren't necessary:
% switch year
%     case 2003
%         yearlabel='May';
%     case 2004
%         yearlabel='Apr';
%     case 2005
%         yearlabel='Apr';
%     case 2006
%         yearlabel='May';
%     case 2007
%         yearlabel='Mar';
%     otherwise
%         yearlabel='Jan';
% end
%
% %mergedpath0 = ['/Volumes/Lab_data/MVCO/FCB/MVCO_' yearlabel num2str(year) '/data/processed/grouped/merged/']; %paths for Swallowtail (Kristen's machine)
% %groupedpath0 = ['/Volumes/Lab_data/MVCO/FCB/MVCO_' yearlabel num2str(year) '/data/processed/grouped/'];
%    mergedpath0 = ['\\sosiknas1\Lab_data\MVCO\FCB\MVCO_' yearlabel num2str(year) '/data/processed/grouped/merged/']; %paths for Swallowtail (Kristen's machine)
%    groupedpath0 = ['\\sosiknas1\Lab_data\MVCO\FCB\MVCO_' yearlabel num2str(year) '/data/processed/grouped/'];

%NOTE: USE OF '0' AT THE END OF MERGED AND GROUPED PATHNAMES IS IMPORTANT -
%a 'mergedpath' name is stored in the data files that are being loaded!!!

%% Create and open a writerObj:
eval(['writerobj1 = VideoWriter(''' savepath 'together_' num2str(year) '.avi'');'])
open(writerobj1);

disp(['Making combined movie for: ' num2str(year)])

%find the files:
filelist_merged=dir([mergedpath0 '*merged_*.mat']);
if year >= 2006
    filelist_grouped=dir([groupedpath0 'FCB*_*.mat']);
else
    filelist_grouped=dir([groupedpath0 '*_*.mat']);
end

% Get the dates that go along with each file:
%if load straight from filelist, the actual dates will be out of order as
%the way the files line up are in order of FCB#...
grouplist=[];
filedate=[];
%%
for q=1:length(filelist_grouped)
    
    groupedfile=filelist_grouped(q).name;
    
    eval(['load(''' groupedpath0 groupedfile ''',''cellresults'')']) %for the date
    
    ind=find(~isnan(cellresults(:,1)));
    if ~isempty(ind)
        grouplist=[grouplist; repmat(groupedfile,length(ind),1)];
        filedate=[filedate; cellresults(ind,1)];
    end
    
    clear cellresults
end


[~, chorder]=sort(filedate); %sort files into chronological order (chorder is the indexing for this)


%%
count=0;
groupedfile = grouplist(chorder(count+1),:); %see what the next file needs to be

eval(['load ' groupedpath0 groupedfile])
disp(['loaded 1st grouped and merged file...' groupedfile])

if  year > 2005
    eval(['load ' mergedpath0 groupedfile(1:end-6) 'merged' groupedfile(end-5:end-4)]) %2006 and up
else %
    eval(['load ' mergedpath0 groupedfile(1:7) 'merged' groupedfile(8:9)])  %   %2005 and below
end

%set up the figures:
fig1=figure(1);
%keyboard
set(fig1,'Position',[24         124        1849         861],'Visible','off')

maxvalue = 1e6;  %is this too high?
bins = 10.^(0:log10(maxvalue)/127:log10(maxvalue));  %make 256 log spaced bins
maxvalueSSC = 1e7;  %is this too high?
binsSSC = 10.^(0:log10(maxvalueSSC)/127:log10(maxvalueSSC));  %make 256 log spaced bins

while count <= length(filedate)-1 %have gone through all the hours...
    %
    %profile on
    count=count+1;
    %for count=1:20
    j=find(cellresults(:,1)==filedate(chorder(count))); %This should be there
    
    if isempty(j)
        disp('Whoops - cannot find the corresponding date in cellresults?')
        keyboard
    end
    
    mwc=mergedwithclass{j}; 


    if any(any(~isnan(mwc)))
        
        %FIGURE 1 - heatmap historgram of where cells are:
        set(0,'CurrentFigure',fig1)
        
        subplot(2,4,1,'replace')
        [n,x] = histmulti5(mwc(:,3:4), [bins' bins']); %forward ligth scatter and chl fluorecence
        indn = find(n == 0); n(indn) = NaN;
        surf(x(:,1), x(:,2), log10(n)')
        ylabel('CHL','Fontsize',14), xlabel('FLS','Fontsize',14)
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        shading flat
        axis([1 1e6 1 1e6])
        view(2)
        title([datestr(cellresults(j,1))],'Fontsize',14)
        %
        
        subplot(2,4,2,'replace')
        [n,x] = histmulti5([mwc(:,5), mwc(:,2)], [binsSSC' bins']); %ssc and PE fluorescence
        indn = find(n == 0); n(indn) = NaN;
        surf(x(:,1), x(:,2), log10(n)')
        ylabel('PE','Fontsize',14), xlabel('SSC','Fontsize',14)
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        shading flat
        view(2)
        axis([1 1e7 1 1e6])
        
        subplot(2,4,3,'replace')
        [n,x] = histmulti5([mwc(:,5), mwc(:,4)], [binsSSC' bins']); %SSC and Chl
        indn = find(n == 0); n(indn) = NaN;
        surf(x(:,1), x(:,2), log10(n)')
        ylabel('CHL','Fontsize',14), xlabel('SSC','Fontsize',14)
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        shading flat
        view(2)
        axis([1 1e7 1 1e6]) % axis([1e2 1e7 1e1 1e6]) for other years
        
        if year >= 2007
            title(['FCB# ' groupedfile(4)])
        else
            title('FCB1')
        end
        
        subplot(2,4,4,'replace')
        indpe=find(mwc(:,end) == 1 | mwc(:,end) == 2 | mwc(:,end) == 3 | mwc(:,end) == 6); %all PE containing classes
        
        [n,x] = histmulti5([mwc(indpe,4), mwc(indpe,2)], [bins' bins']); %PE and Chl
        indn = find(n == 0); n(indn) = NaN;
        surf(x(:,1), x(:,2), log10(n)')
        %hold on
        %loglog([1:10:1e5], 10.^(log10([1:10:1e5])*fit(1) + fit(2)), 'r') %not sure what this is...
        ylabel('PE','Fontsize',14), xlabel('CHL','Fontsize',14)
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        shading flat
        view(2)
        axis([1 1e6 1 1e6])
        %                   %        title('PE containing cells only')
        %clear maxvalue bins
        
      
        %BOTTOM PANELS: 2D plots by call of what particle is:
        numcluster = 6;
        colorstr = ['r', 'b', 'y', 'g', 'k', 'c', 'm'];
              
        subplot(2,4,5,'replace')
        hold on
        ylabel('CHL','Fontsize',14), xlabel('FLS','Fontsize',14)
        for c = 1:numcluster,
            indc = find(mwc(:,end) == c);
            eval(['loglog(mwc(indc,3),mwc(indc,4), ''' colorstr(c) 'o'', ''markersize'', 1)'])
        end;
        
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        axis([1 1e6 1 1e6])
        X = 1:1000:1e6; plot(X,X*5+50000, 'k-')
        %title([datestr(cellresults(j,1))])
        
        subplot(2,4,6,'replace')
        hold on
        ylabel('PE','Fontsize',14), xlabel('SSC','Fontsize',14)
        for c = 1:numcluster,
            indc = find(mwc(:,end) == c);
            eval(['loglog(mwc(indc,5),mwc(indc,2), ''' colorstr(c) 'o'', ''markersize'', 1)'])
        end;
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        axis([1 1e7 1 1e6])

        
        subplot(2,4,7,'replace')
        hold on
        ylabel('CHL','Fontsize',14), xlabel('SSC','Fontsize',14)
        for c = 1:numcluster,
            indc = find(mwc(:,end) == c);
            eval(['loglog(mwc(indc,5),mwc(indc,4), ''' colorstr(c) 'o'', ''markersize'', 1)'])
        end;
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        axis([1 1e7 1 1e6])
        
        subplot(2,4,8,'replace')
        hold on
        ylabel('PE','Fontsize',14), xlabel('CHL','Fontsize',14)
        for c = 1:numcluster,
            indc = find(mwc(:,end) == c);
            eval(['h' num2str(c) '=loglog(mwc(indc,4),mwc(indc,2), ''' colorstr(c) 'o'', ''markersize'', 1);'])
        end;
        %loglog([1:10:1e5], 10.^(log10([1:10:1e5])*fit(1) + fit(2)), 'r')
        set(gca, 'xscale', 'log', 'yscale', 'log','Fontsize',14)
        axis([1 1e6 1 1e6])
        if ~isempty(h1)
            hleg=legend([h1 h2 h3 h4 h5 h6],'syn', 'bright cryptos', 'junk w/pe', 'euks (i.e., no pe)', 'euk junk', 'dim cryptos','Orientation','Horizontal');
            set(hleg,'Position',[0.25    0.48    0.5526    0.0268])
        end;
        
    
         
        F1=getframe(fig1);
%        F1 = im2frame(zbuffer_cdata(fig1)); %an older fix for earlier
%        matlab versions, getframe previously didn't play nice with
%        unvisible objects, I think....
        writeVideo(writerobj1,F1);
        
        
    end
    
    %Histograms
    %Only Syn:
    % figure(3)
    
    % inds = find(mwc(:,end) == 1);
    % hist(log(mwc(inds,5)),256);
    % %hist(log(partialdatmerged(datbins,5)),256);
    % axis([6 10 0 inf])
    % xlabel('SSC')
    % title([datestr(cellresultsall(ind(j),1)) ' culture# ' num2str(port)])
    
    if count ~= length(chorder) %if not last hour, see if another file needs to be loaded:
        
        groupedfile2 = grouplist(chorder(count+1),:); %see what the next file needs to be
        
        %If the next file is not the current file loaded:
        if ~strcmp(groupedfile2,groupedfile) %go ahead and load the necessary grouped and merged files to get the right date
            
            groupedfile=groupedfile2;
            disp(['loading grouped and merged files...' groupedfile '...current frame: ' num2str(count) ' out of ' num2str(length(filedate))])
            eval(['load ' groupedpath0 groupedfile])
            
            if  year > 2005
                eval(['load ' mergedpath0 groupedfile(1:end-6) 'merged' groupedfile(end-5:end-4)]) %2006 and up
            else %
                eval(['load ' mergedpath0 groupedfile(1:7) 'merged' groupedfile(8:9)])  %   %2005 and below
            end
        end
    end
   % end %for for loop for trouble shooting purposes!
end %end of while loop


close(writerobj1);

end %end of function