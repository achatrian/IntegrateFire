clear
clc
close all
%%
%Simulation:
tic
%N: Total number of network neurons;
%Pe: Probability that a neuron is excitatory;
%n: number of timesteps;
%dt: timestep length;
%Random weight matrix:
%pcon: average connectivity;
%wmin, wmax: weight bounds of initial matrix;
%Iamp: Amplitude (average) of external input current (pA);
%tau: time constant = RC
%intmet: integration method (4thRK, Euler,...)
%Ures: reset voltage (mV)
%Vth: threshold voltage (mV)
%R: resistance (MOmhs)
%f: evolution function;
%ta1: spike decay time constant;
%asp = initial amplitude of spike;
%a1: spike function;
%Nin: number of input and output neurons
%Nm: number of neurons in pattern;
%npl: length of pattern;
%nps: pause length;
%unc: uncertainty in pattern
%pattype: pattern type
%Nco: number of coactivated neurons
%Ntco: number of stimualted neurons in test
%learn: evolution function handle - @name of evolution function;
%b: Hebbian learning rate
%wdown: lower weight limit in learning;
%wup: upper weight limit in learning;
%tminus: STDP window time constant -
%tplus: STDP window time constant +
%etam: STDP coeff +
%etap: STDP coeff -

%Choose simulation parameters:
N = 500; Pe = 1; n = 610; dt = 0.01; pcon = 1; wmin = 0.0001;  wmax = 0.06;
Iamp = 100; tau = 28; intmet = 'RK4th'; Ures = -72; Vth = -38; R = 188; 
f = @(u,I) -(u-Ures) + R*I; ta1 = 0.001; asp = 1; Nin = 100; Nm = Nin; 
npl = 300; nps = n/100; unc = 0.001; pattype = 'altern'; Nco = 25; Ntco = 23; 
learn = 'Hblearn'; b = 1000; wdown = 0; wup = 1; tminus = 1; tplus = 1;
etam = 1; etap = 1; Apl = 3; Ami = 3; spec = 'hier';
%Load network parameters
%Programme learning rule sweep:

%Option to interrupt simulation if neuron has not learned pattern to save
%time

%Online computation of statistics
    %average weight change
    %similarity to wanted activation pattern
    
%Initialise variables separately
[W,A,Aus] = RandomW(N,Pe,pcon,wmin,wmax,spec);
Wini = W; %store initial matrix

%Initialise:
[Ne,Ni,t,vth,ures,u,a1,Isp,ft,fr,frate,Wsto,x,y,apl,ami,Tsp,indT] = ...
    Initio(W,Pe,n,dt,Vth,Ures,asp,ta1,Apl,Ami);

%Form pattern
[Iext,Itest,PattNeur] = PattGen(Iamp,W,Pe,n,Nm,npl,nps,unc,pattype,Nco,Ntco);

%Set input and output neurons and change pattern format:
%[W,inpout,Iext,Itest] = MakeInpOut(W,Nin,Iext,PattNeur);

%% Simulation:
for nst=2:n 
      %-------Dynamics--------------------------------------------------
      [u,Isp,fr,ft,frate] = NetDyn(u,vth,ures,W,A,ft,t,frate,tau,nst,Isp,Iext,f,a1,intmet);
      %-------Learning:-------------------------------------------------
      %Hebbian:
      if strcmp(learn,'Hblearn')
        [W,frate] = Hblearn(W,b,frate,fr,nst,wdown,wup);
      elseif strcmp(learn,'STDPsum')
      %STDP-sum:
        [W,Tsp,indT] = STDPdisclearn4for(W,t,nst,fr,Tsp,indT,tminus,tplus,etam,etap,wdown,wup);
      elseif strcmp(learn,'STDPloc')
      %STDP-local
        [W,x,y] = STDPlocallearn(W,x,y,nst,fr,dt,tminus,tplus,apl,ami,wdown,wup);
      end
      %-----------------------------------------------------------------
      %Place W in film:
      Wsto(:,:,nst) = W;
      %Track weight update:
      wchange(nst) = 1/(sum(sum(Aus)))*sqrt(sum(sum((W - Wsto(:,:,nst-1))^2)));
end

ftPatt = ft; %Activation in training
%% Test:
%Initialise:
[Ne,Ni,t,vth,ures,u,a1,Isp,ft,fr,frate,Wsto,x,y,apl,ami,Tsp,indT] = ...
    Initio(W,Pe,n,dt,Vth,Ures,asp,ta1,Apl,Ami); %reset matrices

for nst=2:n 
    [u,Isp,fr,ft,frate] = NetDyn(u,vth,ures,W,A,ft,t,frate,tau,nst,Isp,Itest,f,a1,intmet);
end

timespent = toc; %running time

%% Statistics

[E] = IFstats(W,Wini,ft,ftPatt);

%% Graphics
run('ProcIF.m')
