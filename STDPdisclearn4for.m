function [W,Tsp,indT] = STDPdisclearn4for(W,t,ct,fr,Tsp,indT,tminus,tplus,etam,etap,wdown,wup)
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
Aus = W~=0; %Unsigned Adjacency Matrix

%% Do:
%Store times of spiking:
      Tsp(indT(fr)) = t(ct);
      %Increase index vector
      indT(fr) = indT(fr) + N; 
      %Goes to next entry in the same row - linear indexing with N rows.
      
      %Compute weight:
          for j = 1:N
              for i = find(Aus(:,j))' %loops over neurons post-synaptic to j
                for fj = 1:find(Tsp(j,:)<Inf,1,'last')
                    for ni = 1:find(Tsp(i,:)<Inf,1,'last')
                        T = Tsp(i,ni) - Tsp(j,fj); 
                        %Compute difference between presynaptic and postsynaptic spike
                        %times
                        if T > 0
                            %dw = etap*exp(-T/tplus); %multiply by weight?
                            %dw = W(i,j)*etap*exp(-T/tplus);
                            dw = 0.001;
                        elseif T < 0
                            %dw = - etam*exp(T/tminus);
                            dw = -0.001;
                        else
                            dw = 0;
                        end

                        if W(i,j) + dw > wup
                           W(i,j) = wup;
                        elseif W(i,j) + dw < wdown
                           W(i,j) = wdown;
                        else
                            W(i,j) = W(i,j) + dw;
                        end
                   end
                end
             end
          end
      
   %MAKE SURE ONLY GOOD WEIGHTS ARE UPDATED

end

