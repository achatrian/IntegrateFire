%Plot IF

%% Processing
%Change signs in W:
%{

%Film
%Wmov = immovie(Wsto); NEED RGB

%Plot firing times:
figure
imagesc(ft);
ylabel('Neurons')
xlabel('time-step')
title('Firing times')


%Visualise end matrix
figure
subplot(1,2,2)
imagesc(A.*W)
title('Weight matrix produced by training')
ylabel('Postsynaptic')
xlabel('Presynaptic')
colorbar


%Visualis start matrix
subplot(1,2,1)
imagesc(Wini.*A)
title('Randomly generated weight matrix')
ylabel('Postsynaptic')
xlabel('Presynaptic')
colorbar

figure
subplot(2,1,1)
imagesc(ft)
title('Test Output')
ylabel('Neuron #')
xlabel('Time-step')

subplot(2,1,2)
imagesc(Iext)
title('Training input')
ylabel('Neuron #')
xlabel('Time-step')

%{
figure
title('Test input')
ylabel('Neuron #')
xlabel('Time-step')
imagesc(Itest)
%}

%{
[~,imax] = max(wchange-1);
figure
subplot(1,2,1)
imagesc(Wsto(:,:,imax))
subplot(1,2,2)
imagesc(Wsto(:,:,imax))
%}

%figure
%plot(1:length(x),x(1,:))



%Plot voltage evolution:
figure
subplot(2,1,1)
plot(t, u)
title('Voltage vs time')
ylabel('Voltage (mV)')
xlabel('time (ms)')
legend('Neuron1','Neuron2')


I = Iext + Isp;
subplot(2,1,2)
plot(t,I)
title('Input current vs time')
ylabel('Current (mA)')
xlabel('time (ms)')
legend('Neuron1','Neuron2')

figure
plot(-70:1:-20,f(-70:1:-20,0))
title('Phase plane')
ylabel('du/dt')
xlabel('u')

semilogy(1:length(wchange),abs(wchange))
    %plot(1:length(wchange),abs(wchange))
    title('Average weight change')
    ylabel('wmean')
    xlabel('time-step')
    xlim([0 n]);
%}  
    
figure
subplot(4,1,1)
imagesc(ft)
title('firing times in testing')
subplot(4,1,2)
imagesc(ftPatt)
title('firing times in training')
subplot(4,1,3)
imagesc(Itest)
title('Test current')
subplot(4,1,4)
imagesc(Iext)
title('External current')



%Visualise end matrix
figure
subplot(1,2,2)
imagesc(A.*W)
title('Weight matrix produced by training')
ylabel('Postsynaptic')
xlabel('Presynaptic')
colorbar

%Visualis start matrix
subplot(1,2,1)
imagesc(Wini.*A)
title('Randomly generated weight matrix')
ylabel('Postsynaptic')
xlabel('Presynaptic')
colorbar




