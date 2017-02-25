function [Iext,Itest,PattNeur] = PattGen(Iamp,W,Pe,n,Nm,npl,nps,unc,PattType,Nco,Ntco)
%INPUTS:
%1-Iamp: external current amplitude
%2-N: number of neurons;
%3-Pe: probability that a neuron is excitatory;
%4-n: number of timesteps
%5-Nm: number of neurons in pattern
%6-npl: mean pattern length
%7-nps: mean pause between pattern repetitions (in timesteps)
%8-unc: variance factor of time parameter
%9-PattType: patterntype - string
%10-Nco: - Number of coactivated neurons neurons

%OUTPUTS:
%1-Iext: external input current to train neurons
%2-Itest: external input current to test neurons
%3-PattNeur: neurons in pattern in order of activation

%Keeps percentages of excitatory and inhibitory neurons in pattern as from
%the neural network


%ideas:
%could do n_or_t and have time vector input instead of total number of
%timesteps
%Add noise?

%% Checks
N = length(W); %For checks

if(Pe > 1)
    error('Pe is in [0,1]')
end

if(Nm > N)
    error('Nm is < N')
end

if(nps > n)
    error('tps is < n')
end

if(n <= 0)
    error('n is > 0')
end

if(npl > n)
    error('tpl is < n')
end
%% Initialise
Iext = zeros(N,n);
Itest = Iext;
loopn = 1; %loop neuron
stl = 0; %loop step length
l = 1; %loop counter
cpl = 1; %loop pattern length
pll = 1; %limit to pattern length in current loop
psl = 1; %pause length in loop

%Compute average stimulus length:
nst = max([floor(npl/Nm) 1]);
%Number of neurons in patter shapes length of step activation

%Compute std's for pseudorandom extractions:
sst = nst*unc; %Std for step length
sps = nps*unc; %Std for pause length
spl = npl*unc; %Std for pattern length

%Extract pattern length:
%pll = floor(normrnd(npl,spl)); %Define pattern length
pll = npl;

%Pattern info:
PattNeur = zeros(Nco,floor(pll/nst));

%% Produce Pattern
%PATTERN 1st REPETITION:
if PattType == 'simple'
    loopn = randperm(N,Nco);
    i = 1; %loop variable
    stl = nst;
    while cpl < pll
        PattNeur(:,i) = loopn';
        i = i + 1;
        %stl = ceil(normrnd(nst,sst)); %extract length of injection
        Iext(loopn,cpl:cpl+stl-1) = Iamp; %Write Iext
        cpl = cpl + stl; %update pattern length
        loopn = randperm(N,Nco);
    end
elseif PattType == 'seqcon'
    loopn = randperm(N,Nco);
    i = 1; %loop variable
    while cpl < pll
        PattNeur(:,i) = loopn';
        i = i + 1;
        stl = ceil(normrnd(nst,sst)); %extract length of injection
        Iext(loopn,cpl:cpl+stl-1) = Iamp; %Write Iext
        cpl = cpl + stl; %update pattern length
        postn = mod(find(W(:,loopn)),N) + 1; %Postsynaptic neurons to loopn
        loopn = postn(randperm(length(postn),Nco));
    end
elseif PattType == 'overlp'
    loopn = randperm(N,Nco);
    i = 1; %loop variable
    stl = 2*ceil(normrnd(nst,sst)); %extract length of injection
    cpl = stl/2 + 1;
    while cpl < pll
        PattNeur(:,i) = loopn';
        i = i + 1;
        Iext(loopn,cpl-stl/2:cpl+3/2*stl-1) = Iamp; %Write Iext
        stl = 2*ceil(normrnd(nst,sst)); %extract length of injection
        cpl = cpl + stl; %update pattern length
        %double the time compared to simple pattern
        postn = mod(find(W(:,loopn)),N) + 1; %Postsynaptic neurons to loopn
        loopn = postn(randperm(length(postn),Nco));
        
    end
elseif PattType == 'altern'
    Ngroup = floor(N/Nco); %number of full groups of coactivated neurons
    loopn = 1:Nco; 
    i = 1;
    stl = nst;
    while cpl < pll
        %stl = ceil(normrnd(nst,sst)); %extract length of injection
        Iext(loopn,cpl:cpl+stl-1) = Iamp; %Write Iext
        Itest(loopn(1:Ntco),cpl:cpl+stl-1) = Iamp;
        cpl = cpl + stl; %update pattern length
        PattNeur(:,(cpl-1)/stl) = loopn';
        i = 1 + mod(i,Ngroup); %mod produces number 0 -> Nco-1, so add 1.
        loopn = ((i-1)*Nco+1):i*Nco;
    end
end



%% Propagate Pattern (and test current)
%Fills the rest of Iext with the pattern:
l = cpl; %pattern actually has length cpl - 1
while l+cpl < n
    Iext(:,l:l+cpl-1) = Iext(:,1:cpl);
    Itest(:,l:l+cpl-1) = Itest(:,1:cpl);
    l = l + cpl;
    %psl = floor(normrnd(nps,sps));
    %psl = nps;
    %l = l + psl;
    if l + cpl > n
       Iext(:,l:n) = Iext(:,1:n-l+1);
       Itest(:,l:n) = Itest(:,1:n-l+1);
       l = n; %to exit loop
    end
end

%% Make test current

if PattType ~= 'altern'
    Itest(:,1:cpl) = Iext(:,1:cpl); %Defines the current to test whether neuron has fired.
    %Fills the end with the beginning of the pattern:
    %Iext(:,l:n) = Iext(:,1:n-l+1);
end


end

