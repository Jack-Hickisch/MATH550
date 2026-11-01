clc, clear, close all

font = 24; %fontsize
lw = 4.5; %linewidth
lwlrg = 6; %linewidth large
lwsml = 3; %linewidth small
window_size_for_export = [0, 50, 1920, 975];

a_s = 1.2; % "a scalar" which is variance scalar for 2D Gaussian 
[xdomain,ydomain] = deal([0,2*pi]);
N_range = [10, 100];
N_plot = [10, 30, 100];
N = linspace(N_range(1),N_range(2),(N_range(2)-N_range(1))/10+1); % grid sizes
h = (xdomain(2)-xdomain(1))./(N-1);                              % step sizes
[x,y,X,Y,theta,xib,yib,Sigma] = deal(cell(length(N),1));         % setup x scale for each N

[Lap, Sib, Jib, LHSMat, RHSVec, ResVec, q, f, f2, u, u_actual] = deal(cell(length(N),1)); % discrete derivative matrices and f for each N and solution matrix

for n = 1:length(N)
   Lap{n} = sparse(N(n)^2);
   f{n} = zeros(N(n)^2,1);
   x{n} = linspace(xdomain(1),xdomain(2),N(n));
   y{n} = linspace(ydomain(1),ydomain(2),N(n));
   [X{n},Y{n}] = meshgrid(x{n},y{n});
   u_actual{n} = zeros(length(x{n}));
   for j = 1:length(y{n})
       u_actual{n}(j,:) = cos(2*x{n})*sin(2*y{n}(j));
   end
   theta{n} = linspace(0,2*pi,floor(2*pi/h(n)));
   xib{n} = pi + cos(theta{n});
   yib{n} = pi + sin(theta{n});
   xib{n} = xib{n}(1:end-1); % 0 = 2*pi and only need one coordinate
   yib{n} = yib{n}(1:end-1); % 0 = 2*pi and only need one coordinate
   Sigma{n} = [(a_s*h(n))^2,0;0,(a_s*h(n))^2];
   f2{n} = zeros(length(xib{n}),1);
end

for n = 1:length(N)
    Lap_row = 1;
    for j = 1:length(y{n})
        for i = 1:length(x{n})
            % Laplacian Computation
            if i == 1 || j == 1 || i == length(x{n}) || j == length(y{n}) % on boundary
                Lap{n}(Lap_row, global_index(j,i,N(n))) = -8*h(n)^2;
            else
                Lap{n}(Lap_row, global_index(j  ,i  ,N(n))) = -4;
                Lap{n}(Lap_row, global_index(j  ,i+1,N(n))) =  1;
                Lap{n}(Lap_row, global_index(j+1,i  ,N(n))) =  1;
                Lap{n}(Lap_row, global_index(j  ,i-1,N(n))) =  1;
                Lap{n}(Lap_row, global_index(j-1,i  ,N(n))) =  1;
            end
            % f{n}(Lap_row) = -8*h(n)^2*cos(2*x{n}(i))*sin(2*y{n}(j));
            f{n}(Lap_row) = 0;

            % Sib Computataion
            Sib{n}(Lap_row,:) = s_ib(x{n}(i),y{n}(j),xib{n},yib{n},Sigma{n},h(n));

            Lap_row = Lap_row + 1;
        end
    end

    for i = 1:length(xib{n})
        Jib{n}(i,:) = j_ib(x{n},y{n},xib{n}(i),yib{n}(i),Sigma{n},h(n));
        f2{n}(i) = cos(2*xib{n}(i))*sin(2*yib{n}(i));
    end

    LHSMat{n} = [Lap{n},-Sib{n};Jib{n},zeros(length(xib{n}))];
    RHSVec{n} = [f{n};f2{n}];
    
    % ResVec{n} = LHSMat{n}\RHSVec{n};
    ResVec{n} = gmres(LHSMat{n},RHSVec{n},1000,1e-6,1);
    u{n} = reshape(ResVec{n}(1:N(n)^2),N(n),N(n))';
    q{n} = ResVec{n}(N(n)^2+1:end);
end

function delta_a = s_ib(x,y,xib,yib,Sigma,h)
    delta_a = zeros(1,length(xib));
    for i = 1:length(delta_a)
        delta_a(i) = (2*pi*sqrt(det(Sigma)))^-1*exp(-0.5*([xib(i);yib(i)]-[x;y])'*Sigma^-1*([xib(i);yib(i)]-[x;y]))*h^2; % used AI to find 2D Gaussian equation
    end
end


function j_res = j_ib(x,y,xib,yib,Sigma,h)
    j_res = zeros(1,length(x)*length(y));

    u_row = 1;
    for j = 1:length(y)
        for i = 1:length(x)
            j_res(u_row) = (2*pi*sqrt(det(Sigma)))^-1*exp(-0.5*([xib;yib]-[x(i);y(j)])'*Sigma^-1*([xib;yib]-[x(i);y(j)]))*h^2;
            u_row = u_row + 1;
        end
    end
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
    % xticks(0:1:5)
    % yticks(0:1:5)
    % zticks(-1:0.5:1)
    fontsize(font,"points")
    xlim([0,2*pi])
    ylim([0,2*pi])
    % zlim([-1,1])
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

sgtitle('Problem 3 Solution Plot', 'Interpreter', 'latex'); 

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
exportgraphics(gcf,"HW2_3.jpg")