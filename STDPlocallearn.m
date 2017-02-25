function [W,x,y] = STDPlocallearn(W,x,y,tc,fr,dt,tminus,tplus,aplus,aminus,wdown,wup)
%INPUTS:
%W: weight matrix
%x: vector of presynaptic traces (Nx1)
%y: vector of postsynaptic traces (Nx1)
%tc: current time step
%fr: neurons which fired at this step
%dt: integration time step length
%aminus: multiplicative constant for postsynaptic spikes
%aplus: multiplicative constant for presynaptic spikes
%tminus: window time constant for depression (post->pre)
%tplus: window time constant for potentiation (pre->post)
%wdown: lower weight bound
%wup: upper weight bound
%OUTPUTS:
%W: weight matrix
%x: vector of presynaptic traces (Nx1)
%y: vector of postsynaptic traces (Nx1)

%% Initialise
N = length(W); %# of neurons
Aus = W~=0; %Unsigned adjacency matrix

%% Do
  %Update traces (Euler):
  x(:,tc) = x(:,tc-1) + dt/tplus*(-x(:,tc-1) + aplus.*fr); %Numerical Integ
  y(:,tc) = y(:,tc-1) + dt/tminus*(-y(:,tc-1) + aminus.*fr);
  %Compute weight update (Euler):
  dW = dt*(-y(:,tc).*aplus.*fr*ones(1,N).*aplus.*fr + ones(N,1)*(x(:,tc).*aplus.*fr)' )...
  .*heaviside(wup - W).*heaviside(W - wdown).*Aus;
%Code works because fr*fr' gives symmetric matrix
  %Check for bound-surpassing increments:
  W(W + dW > wup) = wup; %Push to bound;
  W(W + dW < wdown) = wdown; %Push to bound;
  good = W + dW < wup & W + dW > wdown;
  W(good) = W(good) + dW(good);
%introduce Aplus and Aminus? Even more parameters?
%Introduce Euler/RK4 option?
end

