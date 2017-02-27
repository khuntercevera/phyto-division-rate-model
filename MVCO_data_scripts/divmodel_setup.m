
clear
close all

for year2do = 2006:2012;

    % addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/ %has cytosub_SSC2vol.m
    % addpath /Users/kristenhunter-cevera/Documents/phyto-division-rate-model/MVCO_FCB_PROCESSING/

    if ~isempty(strfind(computer,'WIN'))
        envpath='\\sosiknas1\lab_data\mvco\EnvironmentalData\';
        rootpath='\\sosiknas1\lab_data\MVCO\FCB\';
    else
        envpath='/Volumes/Lab_data/MVCO/EnvironmentalData/';
        rootpath='/Volumes/Lab_data/MVCO/FCB/';
    end

    do_Solar = 0;
    %    do_moreSolar = 0; %not ready yet, but one day soon - try to recoop days lost by using Nantucket buoy data :)
    do_setupdays = 0;
    do_setupdays_movie = 0;
    do_model = 1;
    do_modelfit_movie = 0;
    redo_model=0;

    switch year2do
        case 2003
            datapath = fullfile(rootpath,'MVCO_May2003/');
        case 2004
            datapath = fullfile(rootpath,'MVCO_Apr2004/');
        case 2005
            datapath = fullfile(rootpath,'MVCO_Apr2005/');
        case 2006
            datapath = fullfile(rootpath,'MVCO_May2006/');
        case 2007
            datapath = fullfile(rootpath,'MVCO_Mar2007/');
        otherwise
            datapath =fullfile(rootpath,['MVCO_Jan' num2str(year2do) '/']);
    end

    if do_Solar
        solarsavepath=fullfilte(datapath, '/model/');
        solarplotflag=1;
        if ~exist(solarsavepath, 'dir'), mkdir(solarsavepath), end;
        getSolar
        pause
    end

    %     if do_moreSolar %future goal - pad Solar data with buoy solar...
    %         solarsavepath=[datapath '/model/'];
    %         getmoreSolar
    %     end

    if do_setupdays
        mergedpath0 = fullfile([datapath,'data/processed/grouped/merged/']); %the '0' designation is important - mergedpath is a save variable name in the files that will be downloaded....
        groupedpath = fullfile(datapath,'data/processed/grouped/');
        beadpath=fullfile(datapath, 'data/processed/beads/');
        modelpath = fullfile(datapath, 'model/input_picoeuk_Sept2016/');
        if ~exist(modelpath, 'dir'), mkdir(modelpath), end; %where daily input will go....
        plotflag=1;
        setup_days_all
    end

    % make movies of input days:
    if do_setupdays_movie
        MVCO_setup_days_QC
    end

    %run the model:
    if do_model

        disp(num2str(year2do))
        pathname=fullfile(datapath,'\model\input_beadmean_July2016\');
        savepath=fullfile(datapath,'\model\output_Aug2016_onecomp\'); %Change path to sosiknas here!
        if ~exist(savepath,'dir'), mkdir(savepath), end % make directory if doesn't exist yet
        addpath \\sosiknas1\lab_data\mvco\FCB\Syn_divrate_model\onecomp_model_code
        call_to_opt
        %cd onecomp_model_code
        %call_to_opt_field_onecomp
    end

    %make model result figures:
    if do_modelfit_movie
        modelres_path=fullfile(datapath,'/model/output_July2016/');  %path to model results
        setupdays_path=fullfile(datapath,'/model/input_beadmean_July2016/'); %path to setup_days input
        addpath(fullfile(rootpath,'Syn_divrate_model/model_code/')); %path to access code to project model day forward
        addpath(fullfile(rootpath,'Syn_divrate_model/subplot_tight/')); %path for making 'tighter' subplots :)

        modelfit_movies
    end

    %rerun specified additional days
    if redo_model
        disp(['Redoing days in ' num2str(year2do)])
        input_path=fullfile(datapath,'/model/input_beadmean_July2016/'); %path to setup_days input
        savepath=fullfile(datapath,'/model/output_July2016/');  %path to model results
        if ~exist(savepath,'dir'), mkdir(savepath), end %make directory if doesn't exist yet
        [days2redo]=exclude_modeldata(year2do); % get the list of days to redo

        if ~isempty(days2redo)
            days2redo=str2num(cell2mat(days2redo(:,1))); %turn cell array into day list
            call_to_opt_redo120
        end
    end

    %clearvars -except year2do
    % close all
end