function fig = sensfig(frf,i)

    figure(i); clf;

    
    subplot(2,1,1);
    h = loglog(frf.f,abs(squeeze(frf.dydXn(1,:,:)))); hold on;
    title('Normalized parameter sensitivities');
%     subplot(4,1,2);
%     semilogx(frf.f,phase_mod(squeeze(frf.dydXn(1,:,:))))    
     ylabel('d\phi/dX'); 
     xlim([frf.f(1) frf.f(end)]);
     ylim([1e-8 1]);
    legend(h,'x_1','x_2','x_3','x_4','x_5','x_6',1);

    subplot(2,1,2);
    loglog(frf.f,abs(squeeze(frf.dydXn(2,:,:))))
%     subplot(4,1,4);
%     semilogx(frf.f,phase_mod(squeeze(frf.dydXn(2,:,:))))
    ylabel('d\delta/dX');
    xlim([frf.f(1) frf.f(end)]);      ylim([1e-8 1]);
    xlabel('Frequency (Hz)');
    sdf('Latex');
    
    fig.hf = gcf;
end