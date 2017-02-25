function [W,frate] = Hblearn(W,b,frate,fr,nst,wdown,wup)
%INPUTS:
%W: weight matrix
%b: Hebbian learning rate
%frate: vector containing firing rates of all neurons
%fr: vector that's true in entry i if neuron i fired in this timestep
%wdown: lower weight bound
%wup: upper weight bound
%OUTPUTS:
%W: weight matrix
%% Initialise
Aus = W ~= 0; %Unsigned Adjacency Matrix

%% Do
%Hebbian:
      %Update firing frequency ~ activity:
      frate = frate + 1/nst*(fr - frate);

      %%{
      dW = b*(frate*frate'); %Compute weight change
      dW(~Aus) = 0; %Set change for unconnected synapses to 0
      
      %Strengthening connections:
      Wbad = (W + dW > wup); %Increments resulting in bound surpassing
      upd = (fr&fr'); %Synapses to check increments for (NxN)
      W(Wbad&upd) = wup; %Push boundary surpassing weights to upper bound
      upd = (~Wbad)&upd; %Synapses to increment
      W(upd) = W(upd) + dW(upd); %Effect change
      
      %Weakening connections: 
      Wbad = (W - dW < wdown); %Decrements resulting in bound surpassing
      
      %Post fired and pre did not:
      upd = (fr&~fr'); %Synapses to check increments for
      W(Wbad&upd) = wdown; %Push bound-surpassing weights to lower bound
      upd = (~Wbad)&upd; %Synapses to decrement
      W(upd) = W(upd) - dW(upd); %Effect change
      
      %Pre fired and post did not:
      upd = (~fr&fr'); %Synapses to check increments for
      W(Wbad&upd) = wdown; %Push bound-surpassing weights to lower bound
      upd = (~Wbad)&upd; %Synapses to decrement
      W(upd) = W(upd) - dW(upd); %Effect change
end

