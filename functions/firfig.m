function fig = firfig(fir,sys,i,j)

    Fs = 1/diff([fir.tau(1) fir.tau(2)]);
    g_sys = [zeros(fir.m,2);impulse(sys.y,fir.tau(fir.m+1:end))]/Fs;
    
    fir.ylim = 5/3*max(abs(fir.g))'*([-1 1]+1/5);
    
    figure(i); clf;
    subplot(3,1,1:2);
        h(1) = plot(fir.tau,fir.g_raw(:,j),'Color',0.8*ones(1,3)); hold on;
        h(2) = plot(fir.tau,fir.g(:,j),'Color',0.4*ones(1,3)); box off
        h(3) = plot(fir.tau,g_sys(:,j),'k:');
        title('Finite Impulse Response (FIR)');
        legend(h,[fir.legend{j} ' raw'],[fir.legend{j} ' filtered'],...
            [fir.legend{j} ' model'],2); legend('boxoff');
        xlim([fir.tau(1),fir.tau(end)]); ylim(fir.ylim(j,:)); 
        ylabel('angle (rad)');
    subplot(3,1,3); 
        h = area(fir.tau,fir.win); set(h,'FaceColor',0.9*ones(3,1));
        ylim([0 1]); xlim([fir.tau(1),fir.tau(end)]);box off
        legend(h,'window',2);legend('boxoff');
        ylabel('(-)');
    xlabel('\tau (s)');
    sdf('Latex');
    
    fig.hf = gcf;