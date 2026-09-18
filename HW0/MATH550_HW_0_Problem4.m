clc, clear, close all

font = 24; %fontsize
lw = 4.5; %linewidth
lwlrg = 6; %linewidth large
lwsml = 3; %linewidth small
window_size_for_export = [0, 50, 1920, 975];

[xdomain,ydomain] = deal([0,5]);
N_range = [10, 100];
N_plot = [10, 25, 100];
N = linspace(N_range(1),N_range(2),(N_range(2)-N_range(1))/5+1); % grid sizes
h = (xdomain(2)-xdomain(1))./(N-1);                              % step sizes
[x,y,X,Y] = deal(cell(length(N),1));                             % setup x scale for each N

[DDM, f, u, u_actual] = deal(cell(length(N),1)); % discrete derivative matrices and f for each N and solution matrix

for n = 1:length(N)
   DDM{n} = sparse(N(n)^2);
   f{n} = zeros(N(n)^2,1);
   x{n} = linspace(xdomain(1),xdomain(2),N(n));
   y{n} = linspace(ydomain(1),ydomain(2),N(n));
   [X{n},Y{n}] = meshgrid(x{n},y{n});
   u_actual{n} = zeros(length(x{n}));
   for j = 1:length(y{n})
       u_actual{n}(j,:) = sin(x{n})*cos(y{n}(j));
   end
end

for n = 1:length(N)
    DDM_row = 1;
    for j = 1:length(y{n})
        for i = 1:length(x{n})
            if i == 1 || j == 1 || i == length(x{n}) || j == length(y{n}) % on boundary
                DDM{n}(DDM_row, global_index(j,i,N(n))) = -2*h(n)^2;
            else
                DDM{n}(DDM_row, global_index(j  ,i  ,N(n))) = -4;
                DDM{n}(DDM_row, global_index(j  ,i+1,N(n))) =  1;
                DDM{n}(DDM_row, global_index(j+1,i  ,N(n))) =  1;
                DDM{n}(DDM_row, global_index(j  ,i-1,N(n))) =  1;
                DDM{n}(DDM_row, global_index(j-1,i  ,N(n))) =  1;
            end
            f{n}(DDM_row) = -2*h(n)^2*sin(x{n}(i))*cos(y{n}(j));
            DDM_row = DDM_row + 1;
        end
    end
    u{n} = DDM{n}\f{n};
    u{n} = reshape(u{n},N(n),N(n))';
end


function g = global_index(j, i, N)
    g = N*(j-1) + i;
end

figure(1)
i=1;
plots_out_of_order = [1,3,2];
z_label = [strcat("$u_{N=", num2str(N_plot(1)), "}$"), strcat("$u_{N=", num2str(N_plot(2)), "}$"), strcat("$u_{N=", num2str(N_plot(3)), "} \approx u(x,y)$")];
for n = N_plot
    subplot(2,2,plots_out_of_order(i))
    plot_i = find(N == n);
    surf(X{plot_i},Y{plot_i},u{plot_i},FaceAlpha=0.75)
    colormap(slanCM('bamo')); 
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    zlabel(z_label(i), 'Interpreter', 'latex');
    xticks(0:1:5)
    yticks(0:1:5)
    zticks(-1:0.5:1)
    fontsize(font,"points")
    xlim([0,5])
    ylim([0,5])
    zlim([-1,1])
    set(gca, 'TickLabelInterpreter', 'latex')
    ax1 = gca;
    ax1.LineWidth = lwsml*.5;
    i = i + 1;
end



%% Error calcs

[Err_p_2,Err_p_inf] = deal(zeros(length(N),1)); % setup arrays to hold errors

for j = 1:length(N)
    top_diff = u_actual{j} - u{j}; % top of error equation
    top_diff = top_diff(:); %flattening
    Err_p_2(j) = norm(top_diff)/norm(u_actual{j}(:)); %P-2 Norm Error
    Err_p_inf(j) = max(abs(top_diff))/max(abs(u_actual{j}(:))); %P-inf Norm Error
end

sgtitle('Analytic and Finite Difference Solutions to Poisson Equation $\nabla^2 u = f(x,y)$ with the Test Solution $u(x,y)=\sin(x)\cos(y)$', 'Interpreter', 'latex'); 

subplot(2,2,4)
loglog(N,Err_p_inf,LineWidth=lw,Color='#7E2F8E')
hold on
loglog(N,Err_p_2,LineWidth=lw,Color='#EDB120')
% axis square
legend('$P=\infty$','$P=2$','Interpreter', 'latex',Location='southwest')
xlabel('$N$', 'Interpreter', 'latex');
ylabel('Error', 'Interpreter', 'latex');
title('Relative Error Convergence', 'Interpreter', 'latex');
fontsize(font,"points")
grid on
box on
set(gca, 'TickLabelInterpreter', 'latex')
ax2 = gca;
ax2.LineWidth = lwsml;

set(gcf, 'Position', window_size_for_export);
exportgraphics(gcf,"HW0_4.jpg")
