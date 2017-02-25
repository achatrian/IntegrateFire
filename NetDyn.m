function [u,Isp,fr,ft,frate] = NetDyn(u,vth,ures,W,A,ft,t,frate,tau,nst,Isp,Iext,f,a1,method)
%INPUTS:
%u: Membrane voltage matrix
%vth: Voltage threshold value
%ures: Voltage reset value
%W: Weight Matrix
%A: Signed Adjacency matrix with exc (+ve weight) and inh (-ve) neurons.
%fr: Fired neurons
%t: Time vector
%tau: Membrane time constant
%nst: current time step
%Isp:  Synaptic Input Current matrix
%Iext: External Input Current matrix
%a1: Spike decay function
%method: string describing numerical integration method, e.g. RK4th or Euler
%OUTPUTS:
%u: Membrane voltage matrix
%Isp: Update Ext Current
%fr: fire neurons
%ft: matrix containing all firing times;

%% Initialise:

dt = t(nst) - t(nst-1);
n = length(t);

if method == 'RK4th'
    N = length(W);
    k = zeros(N,4);
    %approximator matrix
end

%% Run dynamics:
      I = Iext + Isp; %total current
        
      %Voltage is updated first since it uses points at previous time point
      if method == 'Euler'
        %Update firing
        %u(:,i) = u(:,i-1) + (h/6)/tau*(k1+2*k2+2*k3+k4);
        u(:,nst) = u(:,nst-1) + dt/tau*f(u(:,nst-1),I(:,nst-1)); 
        %Compute spike-induced current:
      elseif method == 'RK4th'
        k(:,1) = f(u(:,nst-1),I(:,nst-1)); %~fire are the neurons which have not fired
        k(:,2) = f(u(:,nst-1) + dt/2*k(:,1),I(:,nst-1));
        k(:,3) = f(u(:,nst-1) + dt/2*k(:,2),I(:,nst-1));
        k(:,4) = f(u(:,nst-1) + dt*k(:,3),I(:,nst-1));
        u(:,nst) = u(:,nst-1) + dt/tau*(k(:,1) + 2*(k(:,2) + k(:,3)) + k(:,4));
      end
      
      % Update firing time record:
      fr = u(:,nst) > vth; %true values corresponding to neurons that have fired at current time step
      % NB: if no neuron exceeded the threshold, fire is all false.

      Isp(:,nst:n) = Isp(:,nst:n) + sum(W(:,fr).*A(:,fr),2)*a1(t(nst:n),t(nst)); % #nf x (n-i) matrix
      %sums over columns of W corresponding to neurons which fired

      %Store firing times: if neuron fired at nth time step, then nth entry
      %is a true:
      ft(:,nst) = fr;

      %Resets values of neurons which fired:
      u(fr,nst) = ures(fr);
      
      %Updates firing rate
      frate = frate + 1/nst*(fr - frate); 
end

