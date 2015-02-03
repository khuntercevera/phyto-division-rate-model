function [mu mu1 mu2 p1 p2]=growth_rate(Einterp,volbins,N_dist,theta,hr1,hr2)

%Calculate the composite growth rate (division rate) of two subpopulations from specified parameter values in theta
%Also calculates growth rates of each individual subpopulations and their ending proportions

%Einterp - interpolated light data for every 10 min of the day (W/m2)
%N_dist - observed number of counts of cells in each size class as sepcified by volbins
%volbins - cell size classes (micrometers cubed)
%hr1 and hr2 refer to the starting and ending hour of the portion of day
%In the paper, we've been using hr1=7 hours after dawn and hr2=25 (run till end of day). 
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


%% Create all transition matrices for the day: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=hr2-hr1;
B_day1=zeros(57,57,q);   
B_day2=zeros(57,57,q);

for t=(hr1-1):(hr2-2)
     B1=matrix_const(t,Einterp,volbins,b1,dmax1,E_star1,gmax1);
     B_day1(:,:,t-hr1+2)=B1;
     B2=matrix_const(t,Einterp,volbins,b2,dmax2,E_star2,gmax2);
     B_day2(:,:,t-hr1+2)=B2;
end

%% Project forward starting size distributions: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup starting size distributions: first hour size distribuiton is a mixture of two log normal distributions
%Note construction using normal distribtuion as size classes are logarithically spaced
y1=normpdf(1:57,m1,sigma); 
y2=normpdf(1:57,m2,sigma);

%proportion size distributions accordingly with starting cell number from observations and proportion parameter
Nt1=f*sum(N_dist(:,hr1))*(y1./sum(y1)); 
Nt2=(1-f)*sum(N_dist(:,hr1))*(y2./sum(y2));

%Nt1 and Nt2 contain count data
%rearrange for matrix multiplication
Nt1=Nt1'; %counts in each size bin for starting hour
Nt2=Nt2';

%create model day by projecting forward the counts for each hour
for t=1:q
    
    Nt1(:,t+1)=B_day1(:,:,t)*Nt1(:,t);           
    Nt2(:,t+1)=B_day2(:,:,t)*Nt2(:,t);
    
end

%note that because we do not allow division to occur for the first 6 hours
%of the day and that this matches the time span with which we have modeled the cell size
%distributions, we have implicitly divided by t=1 day in the below growth rate calculation.
%If you were to shorten the part of the day that is modeled past the time when cell division is
%allowed, you would need to divide by the appropriate fraction of the day.

mu=(log(sum((Nt1(:,end)+Nt2(:,end)))./sum(Nt1(:,1)+Nt2(:,1))));
mu1=log(sum((Nt1(:,end)))./sum(Nt1(:,1)));
mu2=log(sum((Nt2(:,end)))./sum(Nt2(:,1)));

p1=sum(Nt1(:,end))./sum(Nt1(:,end)+Nt2(:,end));
p2=sum(Nt2(:,end))./sum(Nt1(:,end)+Nt2(:,end));
