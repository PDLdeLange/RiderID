function fig = firfig(npm,mod,i,j)

    u = zeros(size(npm.t(1:npm.m))); u(1) = 1/(npm.t(2)-npm.t(1));
    ghat = lsim(mod.y,u,npm.t(1:npm.m));
    g = npm.g(1:npm.m,:);
    g_raw = npm.g_raw(1:npm.m,:);

    
    tau = npm.t(1:npm.m);
    npm.ylim = 5/3*max(abs(npm.g))'*([-1 1]+1/5);
    
    leg = {'\phi','\delta'};
    ya = {0.007,0.03};

    figure(i); clf;
        h(1) = plot(tau,g_raw(:,j),'Color',0.8*ones(1,3)); hold on;
        h(2) = plot(tau,g(:,j),'Color',0.4*ones(1,3)); box off
%         h(3) = plot(tau,ghat(:,j),'k-.');
        title('Finite Impulse Response');
        legend(h,[leg{j} ' raw'],[leg{j} ' filtered'],...
            2); legend('boxoff'); % leg{j} ' model']
        xlim([tau(1),tau(end)]); ylim(ya{j}*([-1 1]+1/3)); 
        ylabel('angle (rad)');
    xlabel('\tau (s)');
    sdf('Latex');
    
    fig.hf = gcf;