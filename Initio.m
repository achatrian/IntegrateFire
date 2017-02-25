function [Ne,Ni,t,vth,ures,u,a1,Isp,ft,fr,frate,Wsto,x,y,apl,ami,Tsp,indT,wchange] = Initio(W,Pe,n,dt,Vth,Ures,asp,ta1,Apl,Ami)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%Number of neurons
N = length(W);

Ne = floor(N*Pe); %Number of excitatory neurons in network
Ni = N - Ne; %Number of inhibitory neurons in network

t = linspace(0,n*dt,n); %time - for defining nonconstant Iext's

%Voltage reset and thresholds
ures = Ures*ones(N,1);
vth = Vth*ones(N,1);

%Voltage matrix:
u = zeros(N,n);
u(:,1) = ures;

%Spike function:
%Amplitude init
a1 = @(t,tf) asp*exp(-(t-tf)/ta1); %spike function

%Spike-induced current, initially white noise:
Isp = zeros(N,n);

%Matrix ticking the spiking times for each neuron
ft = false(N,n);

%Logical array, true for entries of neurons which did not fire at previous step
fr = false(1,N); %a true corresponds to not having fired at the previosu step

%Firing rate calculation:
frate = zeros(N,1);
    
%W storage:
Wsto = zeros(N,N,n);
Wsto(:,:,1) = W;

%STDPon - two local variables:
x = zeros(N,n); %presynaptic traces
y = zeros(N,n); %postsynaptic traces
apl = Apl*ones(N,1); %increment of x trace after spike
ami = Ami*ones(N,1); %increment of y trace after spike

%STDP - Spike matrix:
Tsp = Inf(N,n);
indT = 1:N; %indT/N gives the entries of each row

end

