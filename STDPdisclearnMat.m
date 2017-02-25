function [W,Tsp,indT] = STDPdisclearnMat(W,t,ct,fr,Tsp,indT,tminus,tplus,etam,etap,wdown,wup)
%INPUTS:
%W: weight matrix
%t: time vector
%ct: current time step
%Tsp: spike matrix
%indT: vector with index of next spike in matrix
%etam: multiplicative constant for depression
%etap: multiplicative constant for potentiation
%tminus: window time constant for depression (post->pre)
%tplus: window time constant for potentiation (pre->post)
%wdown: lower weight bound
%wup: upper weight bound
%OUTPUTS
%W: weight matrix
%Tsp: spike matrix
%indT: vector with index of next spike in matrix

%% Initialise:
%STDP - windowed time matrix
N = length(W);
n = length(t);
ExpT = zeros(N,n);
inc = false(N,1);
exce = inc;
Aus = W~=0; %Unsigned Adjacency Matrix
dw = zeros(N,1); %weight increase

%% Do:
%Store times of spiking:
      Tsp(indT(fr)) = t(ct);
      %Increase index vector
      indT(fr) = indT(fr) + N; 
      %Goes to next entry in the same row - linear indexing with N rows.
      
      %Compute weight:
      for j = 1:N
        for fj = 1:find(Tsp(j,:)<Inf,1,'last')
          T = Tsp - Tsp(j,fj); 
          %Compute difference between presynaptic and postsynaptic spike
          %times
          ExpT(T >= 0 & T < Inf & Aus(:,j)&(ones(1,n))) = etap*exp(-T(T >= 0 & T < Inf & Aus(:,j)&(ones(1,n)))/tplus);
          ExpT(T < 0 & Aus(:,j)&(ones(1,n))) = - etam*exp(T(T < 0 & Aus(:,j)&(ones(1,n)))/tminus);
          dw = dw + sum(ExpT,2);
        end
        
        %Initialise increment indexer with neurons postsynaptic to j:
        inc(Aus(:,j)) = true;
        
        %Update only rows that wouldn't exceed limit, push to limit
        %otherwise:
        exce = W(:,j) + dw > wup; %Weights exceeding upper bound after increment
        inc(exce) = false;
        exce(j) = false; %Don't update self weight
        W(exce,j) = wup; %Push to limit
        
        exce = W(:,j) + dw < wdown; %Weights exceeding lower bound after increment
        inc(exce) = false;
        exce(j) = false;
        W(exce,j) = wdown; %Push to limit
        
        %Update only the weights that were nonzero in original connectivity
        %matrix
        W(inc,j) = W(inc,j) + dw(inc);
        
        dw(:) = 0; %Reset increment
        ExpT(:) = 0;%Reset window matrix
      end
      
   %MAKE SURE ONLY GOOD WEIGHTS ARE UPDATED

end

