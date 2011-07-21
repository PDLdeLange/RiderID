function fig = resfig(fil,res,i,j,T)

% Davis results check
    figure(i); clf; 
    subplot(5,1,1:2);
        h(1) = plot(res.t,fil.y(:,j),'Color',0.5*ones(3,1)); hold on % Total
        h(2) = plot(res.t,res.y(:,j),':k'); % Deterministic
        legend(h,[res.legend{j} '_{measured}'],[res.legend{j} '_{mod}'],2); legend('boxoff');
        xlim(T); ylim(res.ylim(j,:));
        title('Measurement Results'); box off; ylabel('angle (rad)');
    subplot(5,1,3:4); 
        h = plot(res.t,res.e(:,j),'k');  hold on;% Noise
        legend(h,'e',2); legend('boxoff');
        xlim(T); ylim(res.ylim(j,:));
        box off; ylabel('angle (rad)');
    subplot(5,1,5);
        h = plot(res.t,res.f,'k');
        xlim(T); ylim(150*[-1 1]-100);
        box off; ylabel('force (N)'); 
        legend(h,'f',3);legend('boxoff');
    xlabel('time (s)');
    sdf('Latex');
    
    fig.hf = gcf;