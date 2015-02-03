function [prob]=negloglike_calc(Einterp,N_dist,theta,volbins,hr1,hr2)

%Calculates negative log likelihood of a set of parameters given a day of
%observations (counts of cells in each size class for each hour). Model version is a two subpopulation model structure 
% as described in Hunter-Cevera et al. 2014. Subpopulations are
% distinguished based on starting mean volume. The subpopulation with the
% smaller volume is referred to subpopn 1

%The inputs are as follows:

%Einterp - interpolated light data for every 10 min of the day (W/m2)
%N_dist - number of counts of cells in each size class as specified by volbins
%volbins - cell size classes (micrometers cubed)
%hr1 and hr2 refer to the starting and ending hour of the portion of day
%you want to fit. In the paper, we've used hr1=7 hours after dawn and
%hr2=25 (run till end of day). 
%theta - set of parameters, described below:

gmax1=theta(1); %max fraction of cells growing into next size class, subpopn 1
b1=theta(2);  %shape parameter division function, subpopn 1
E_star1=theta(3); %shape parameter of growth function (point where function switches from linear to constant), subpopn 1
dmax1=theta(4); %max fraction of cells able to divide in a given size class, subpopn 1
gmax2=theta(5); %max fraction of cells growing into next size class, subpopn 2
b2=theta(6); %shape parameter division function, subpopn 2
E_star2=theta(7); %shape parameter of growth function (point where function switches from linear to constant), subpopn 2
dmax2=theta(8); %max fraction of cells able to divide in a given size class, subpopn 2
f=theta(9); %proportion parameter, specifies starting fraction of subpopn 1 
m1=theta(10); %mean volume for starting cell size distribution, subpopn 1
m2=theta(11); %mean volume for starting cell size distribution, subpopn 2
sigma=theta(12); %variance parameter for starting cell size distributions for both popn 1 and 2
s=theta(13); %overdispersion parameter for the Dirichlet-multinomial distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create transition matrices for all hours (hr2-hr1 or 24hrs):

q=hr2-hr1; 
B_day1=zeros(57,57,q);   
B_day2=zeros(57,57,q);

for t=(hr1-1):(hr2-2)
     B1=matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1); %matrix construction function
     B_day1(:,:,t-hr1+2)=B1;
     B2=matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2);
     B_day2(:,:,t-hr1+2)=B2;
end


%% Project forward each subcomponent: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup starting size distributions: first hour size distribution is a mixture of two log normal distributions
%Note construction using normal distribution as size classes are logarithmically spaced
y1=normpdf(1:57,m1,sigma); 
y2=normpdf(1:57,m2,sigma);

%proportion size distributions accordingly with starting cell number from observations and proportion parameter
Nt1=f*sum(N_dist(:,hr1))*(y1./sum(y1)); 
Nt2=(1-f)*sum(N_dist(:,hr1))*(y2./sum(y2));

%Nt1 and Nt2 contain count data, whereas the variable simdist (below) are
%the cell size distributions used in the likelihood calculation. We project
%using numbers and then calculate the distribution.

%rearrange for matrix multiplication
Nt1=Nt1'; %counts in each size bin for starting hour
Nt2=Nt2';
simdist(:,1)=(Nt1+Nt2)./sum(Nt1+Nt2); %first hour cell size distribution, 

%create model day by projecting forward the counts for each hour
for t=1:q
    
    Nt1(:,t+1)=B_day1(:,:,t)*Nt1(:,t);           
    Nt2(:,t+1)=B_day2(:,:,t)*Nt2(:,t);
    
    simdist(:,t+1)=(Nt1(:,t+1)+Nt2(:,t+1))./sum(Nt1(:,t+1)+Nt2(:,t+1)); %normalize to get distribution for likelihood
     
    %just in case something is awry:
%     if any(isnan(simdist(:,t+1))) 
%         disp(['simdist contains a nan... theta:' num2str(theta)])
% %         keyboard 
%     end

end

%% Calculate log likelihood using the Dirichlet Multinomial distribution: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We directly calculate the log likelihood value, which allows us to avoid mulitpying small likelihood values:

%specify the expected distribution from simdist:
alpha=s*simdist; %s is overdispersion parameter
TotN=sum(N_dist); %Total cell counts for each hour

logL=zeros(q+1,1);
for t=1:q+1
    C = gammaln(s) - gammaln(TotN(:,t+hr1-1)+s); %constant out in front
    logL(t)=C+sum(gammaln(N_dist(:,t+hr1-1)+alpha(:,t)) - gammaln(alpha(:,t)));
end

%Note that there is another constant that is technically part of the likelihood function 
%that we don't incorporate into this calcution (N(t)!/n1(t)!n2(t)!...n57(t)!, where N(t) is 
%total number of cells obs at time t and n(t) is number of cells in each size class for time t) 
%as this value remains the same over the course of the day. Since it does does not depend on theta, 
%we can exclude it and it will not impact on the search for the best fit parameters

prob=-sum(logL); %return negative logL for the optimization routine fmincon (finds minima, not maxima)

