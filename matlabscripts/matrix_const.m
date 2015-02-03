function B=matrix_const(hr,Einterp,volbins,b,dmax,E_star,gmax)

%Returns hourly transition matrix B for hour hr, which is in turn
%constructed from 10-min interval transition matrices: 
%A(t) is transition matrix for each 10 minute time interval within the hour 

%hr is hour 
%Einterp is interpolated light data over the entire day in 10 min intervals
%b, dmax and E_star and gmax are the parameters that specify the division and growth functions for a population

%Other constants needed for transition matrix construction:
m=length(volbins);
j = find(2*volbins(1) == volbins); %size class for which cells in which division would result in cells smaller than vmin; j=1+(1/dv);  dv=0.125
dt=(1/6); %corresponds to 10 mins (units are hours)
ts=6; %hours after dawn in which division is allowed to occur

%-----------Division function---------------------
%del=frac of cells that divide in each size class per time step dt
del=(dmax.*volbins.^b)./(1+(volbins.^b)); 

%-----------Growth function---------------------
%y is fraction of cells that grow into next class, depends on incident radiation (E), which is dependent on time
y=gmax*ones(size(Einterp));
ind=find(Einterp < E_star);
y(ind)=(gmax/E_star) * Einterp(ind);

%-----------%construct trastion matrix---------------%
A=spalloc(m,m,(m+2*(m-1))); %allocate sparse matrix for A, (size of matrix (m x m) and how many nonzeros will be in A
stasis_ind=1:(m+1):m^2;  %indexing for matrix - goes down a column and then right to next column
growth_ind=2:(m+1):m^2;  %corresponds to growth assignments
div_ind=((j-1)*m)+1:(m+1):m^2;  %division assignments 

for t=1:(1/dt)
    
     if hr <= ts %don't allow cells to divide if before ts
     delta=zeros(1,m);
     else
     delta=del; %allow cells to divide
     end 
     
     A(stasis_ind)=(1-delta)*(1-y(t+6*hr));         % stasis, the 6*hr part in the indexing is because each hour is to match Einterp, which is broken up into 10min segments
     A(m,m)=(1-delta(m));
                 
     A(growth_ind)=y(t+6*hr)*(1-delta(1,1:m-1));    %growth on sub-diagonal
        
     A(1,1:j-1)=A(1,1:j-1)+2*delta(1:j-1);           %division on superdiagonal and on first row partial
     A(div_ind)=2*delta(j:m);
     
     if t ==1
         B=A;
     else
        B=A*B;
     end
     A = A.*0; %reinitialize A for next loop
%      if isnan(B),
%          disp('B has a NaN!')
% %                 keyboard, 
%      end
end

