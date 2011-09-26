function fig = frffig(mod,npm,i,j,F)

figure(i); clf; 

ymod = squeeze(freqresp(mod.y,2.*pi.*npm.f)).';

subplot(2,1,1);
    loglog(npm.f,abs(ymod(:,j)),'k');
    xlim([F(1) F(end)]);
    ylabel('abs \phi');
subplot(2,1,2)
    semilogx(npm.f,unwrap(angle(ymod(:,j))),'k');
    xlim([F(1) F(end)]);
    ylabel('abs \phi');
%     ylim(2*pi*[-1 1]);
    ylabel('angle \phi');
xlabel('Frequency (Hz)');
sdf('Latex');
% legend(h,'H\phi');



fig.hf = gcf;
