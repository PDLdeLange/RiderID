function fig = npmfig(npm,dat,i,j,T)

    
    leg = {'\phi','\delta'};
    lim = {0.1,0.5};


% Davis npmults check
    figure(i); clf; 
    subplot(3,1,1);
        h = plot(dat.t,dat.y(:,j),'k'); hold on % Total
        legend(h,'y(t)',2); legend('boxoff');
        xlim(T); ylim(lim{j}*([-1 1]+1/3));
        title('Output decomposition: y(t) = G(q)w(t) + v(t)'); box off; ylabel('angle (rad)');
    subplot(3,1,2);
        h = plot(npm.t,npm.y(:,j),'k'); % Deterministic
        legend(h,'G(q)w(t)',2); legend('boxoff');
        xlim(T); ylim(lim{j}*([-1 1]+1/3));
        title('='); box off; ylabel('angle (rad)');
    subplot(3,1,3); 
        h = plot(npm.t,npm.v(:,j),'k');  hold on;% Noise
        legend(h,'v(t)',2); legend('boxoff');
        xlim(T); ylim(lim{j}*([-1 1]+1/3));
        title('+'); box off; ylabel('angle (rad)');
%     subplot(4,1,4);
%         h = plot(dat.t,dat.w,'k');
%         xlim(T); ylim(150*[-1 1]-100);
%         box off; ylabel('force (N)'); 
%         legend(h,'f',3);legend('boxoff');
    xlabel('time (s)');
    sdf('Latex');
    
    fig.hf = gcf;