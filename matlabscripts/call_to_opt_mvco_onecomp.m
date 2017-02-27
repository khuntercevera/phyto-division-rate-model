%main script to process both benchtop and field data using the dirichlet
%multinomial distribution:

%script is set up to run batches of random startpoints using multistart as
%this seems to be faster than running fmincon one solver run at time to
%find the answer:

%some front matter:
hr1=7; hr2=25;
ms=MultiStart('Display','off','TolX',1e-6,'UseParallel','always','StartPointsToRun','bounds');
opts=optimset('Display','off','TolX',1e-8,'Algorithm','interior-point','UseParallel','always','MaxIter', 3000,'MaxFunEvals',10000);
icsTol=0.2;
tolvec=[0.01 0.01 100 0.005 0.5 0.5 10];

restitles={'day';'gmax';'b';'E*';'dmax';'m1';'sigma';'s';'-logL';'mu';'exitflag';'number solver runs'};
notes='E* bounds are from 0 to max(Einterp), 7 param-one component model, piece-wise linear gamma function';

% For field data:

filelist = dir([pathname 'day*data.mat']);

modelresults=zeros(length(filelist),12);
allmodelruns=cell(length(filelist),2);


for filenum=1:length(filelist)

    filename=filelist(filenum).name;
    day=str2num(filename(4:9));
    disp(['optimizing day: ' num2str(day) ' file#: ' num2str(filenum)])

    eval(['load ' pathname filename])


    if size(N_dist,2) < 25
        m=size(N_dist,2);
        N_dist=[nan(57,25-m) N_dist];
        Vhists=[nan(57,25-m) Vhists];
    end

    %Interpolate Light Data:
    time=0:(1/6):25;
    nnind = find(~isnan(Edata(:,2)));
    Einterp = interp1(Edata(nnind,1),Edata(nnind,2),time);
    Einterp(find(isnan(Einterp))) = 0;

    %Set bounds:
    lb=-[1e-4 1e-4 1e-4 1e-4 10 1 1e-4];
    ub=[1 15 max(Einterp) 1 50 10 1e4];

    a1=-1*eye(7);
    a2=eye(7);
    A=zeros(14,7);
    A(1:2:13,:)=a1;
    A(2:2:14,:)=a2;

    B=zeros(14,1);
    B(1:2:13)=lb;
    B(2:2:14)=ub;


    x0=[0.2*rand 6*rand max(Einterp)*rand 0.1*rand 30*rand+20 10*rand+2 1e4*rand];

    tpoints = CustomStartPointSet([0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 30*rand(40,1)+20 10*rand(40,1)+2 1e4*rand(40,1)]);

    problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) loglike_DMN_7params(Einterp,N_dist,theta,volbins,hr1,hr2),'Aineq',A,'bineq',B,'options',opts);
    [xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);

    %open up the soln sturcutre:
    temp=zeros(40,10);
    start_points=zeros(40,7);
    c=1;

    for j=1:length(soln)
        %check to see if all start points led to an individual solution or
        %not (MultiSTart will only return unique solutions)
        g=cell2mat(soln(j).X0);
        if length(g)==7 %only one start_point led to that solution
            start_points(c,:)=g;
            temp(c,1:7)=soln(j).X;
            temp(c,8)=soln(j).Fval;
            temp(c,9)=growth_rate_phours_6params(Einterp,volbins,N_dist,temp(c,1:6),hr1,hr2);
            temp(c,10)=soln(j).Exitflag;
            c=c+1;
        else
            num=length(g)/7;
            start_points(c:c+num-1,:)=squeeze(reshape(g',1,7,num))';
            temp(c:c+num-1,1:7)=repmat(soln(j).X,num,1);
            temp(c:c+num-1,8)=repmat(soln(j).Fval,num,1);
            temp(c:c+num-1,9)=repmat(growth_rate_phours_6params(Einterp,volbins,N_dist,temp(c,1:6),hr1,hr2),num,1);
            temp(c:c+num-1,10)=repmat(soln(j).Exitflag,num,1);
            c=c+num;
        end
    end
    %just in case have rows left as zeros
    qq=find(temp(:,1)~=0);
    temp=temp(qq,:);

    modelfits=temp;
    start_points=start_points(qq,:);
    allstarts=start_points;

    %let's now ask, in the first batch run, did the solver "converge"?
    [sortlogL, ii]=sort(temp(:,8));

    if abs(sortlogL(5)-sortlogL(1)) < icsTol
        flag1 = 0;
    else
        disp(num2str(sortlogL(1:5)))
        flag1 = 1;
    end;

    partol=max(modelfits(ii(1:5),1:7))-min(modelfits(ii(1:5),1:7));
    if sum(abs(partol) < tolvec)==7 || sum((abs(partol./modelfits(ii(1),1:7)) < 0.05))==7 %either the modelfits are within an absolute tolerance or within a relative tolerance
        flag2 = 0;
    else
        flag2 = 1;
    end


    disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

    k=1; %batch number
    while (flag1 || flag2) && k <= 5

        disp(['k: ' num2str(k)])
        k=k+1;
        x0=[0.2*rand 6*rand max(Einterp)*rand 0.1*rand 30*rand+20 10*rand+2 1e4*rand];

        tpoints = CustomStartPointSet([0.2*rand(40,1) 6*rand(40,1) max(Einterp)*rand(40,1) 0.1*rand(40,1) 30*rand(40,1)+20 10*rand(40,1)+2 1e4*rand(40,1)]);

        problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) loglike_DMN_7params(Einterp,N_dist,theta,volbins,hr1,hr2),'Aineq',A,'bineq',B,'options',opts);
        [xmin,fmin,exitflag,~,soln] = run(ms,problem,tpoints);

        %open up the soln sturcutre:
        temp=zeros(40,10);
        start_points=zeros(40,7);
        c=1;
        temp=zeros(40,10);
        zeros(40,10);
        c=1;
        for j=1:length(soln)
            %check to see if all start points led to an individual solution or
            %not (MultiSTart will only return unique solutions)
            g=cell2mat(soln(j).X0);
            if length(g)==7 %only one start_point led to that solution
                start_points(c,:)=g;
                temp(c,1:7)=soln(j).X;
                temp(c,8)=soln(j).Fval;
                temp(c,9)=growth_rate_phours_6params(Einterp,volbins,N_dist,temp(c,1:6),hr1,hr2);
                temp(c,10)=soln(j).Exitflag;
                c=c+1;
            else
                num=length(g)/7;
                start_points(c:c+num-1,:)=squeeze(reshape(g',1,7,num))';
                temp(c:c+num-1,1:7)=repmat(soln(j).X,num,1);
                temp(c:c+num-1,8)=repmat(soln(j).Fval,num,1);
                temp(c:c+num-1,9)=repmat(growth_rate_phours_6params(Einterp,volbins,N_dist,temp(c,1:6),hr1,hr2),num,1);
                temp(c:c+num-1,10)=repmat(soln(j).Exitflag,num,1);
                c=c+num;
            end
        end
        %just in case have rows left as zeros
        qq=find(temp(:,1)~=0);
        temp=temp(qq,:);
        start_points=start_points(qq,:);

        modelfits=[modelfits; temp];
        allstarts=[allstarts; start_points];

        %okay, now see after this batch run, did the solver "converge"?
        [sortlogL, ii]=sort(modelfits(:,8));

        if abs(sortlogL(5)-sortlogL(1)) < icsTol
            flag1 = 0;
        else
            disp(num2str(sortlogL(1:5)))
            flag1 = 1;
        end;

        partol=max(modelfits(ii(1:5),1:7))-min(modelfits(ii(1:5),1:7));
        if sum(abs(partol) < tolvec)==7 || sum((abs(partol./modelfits(ii(1),1:7)) < 0.05))==7 %either the modelfits are within an absolute tolerance or within a relative tolerance
            flag2 = 0;
        else
            flag2 = 1;
        end

    end
    %
    [s, jj]=sort(modelfits(:,8));
    xmin=modelfits(jj(1),1:7);
    fmin=modelfits(jj(1),8);
    exitflag=modelfits(jj(1),10);

    [mu]=growth_rate_phours_6params_plt(Einterp,volbins,N_dist,xmin(1:6),hr1,hr2);

    modelresults(filenum,:)=[day xmin fmin mu exitflag length(modelfits)];
    allmodelruns{filenum,1}=modelfits;
    allmodelruns{filenum,2}=allstarts;

    eval(['save mvco_7par_dmn_' num2str(year2do) ' modelresults allmodelruns'])
end

eval(['modelresults_one' num2str(year2do) '=modelresults;'])
eval(['allmodelruns_one' num2str(year2do) '=allmodelruns;'])
eval(['save ' savepath 'mvco_7par_dmn_' num2str(year2do) ' modelresults_one' num2str(year2do) ' allmodelruns_one' num2str(year2do) ' restitles notes'])
