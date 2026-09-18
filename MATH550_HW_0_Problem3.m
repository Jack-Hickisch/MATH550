clc, clear, close all

font = 24; %fontsize
lw = 4.5; %linewidth
lwlrg = 6; %linewidth large
lwsml = 3; %linewidth small
window_size_for_export = [0, 50, 1920, 975];
analytic_color = 'k'; %black
colors = {'#d300a5', '#00a5d3', '#a5d300'}; %triadic colors

domain = [0,5];
N_range = [5, 100];
N = linspace(N_range(1),N_range(2),(N_range(2)-N_range(1))/5+1); % grid sizes
N_plot = [5, 10, 100]; % choose 3 values of N to show on the plot (must be in N_range and divisible by 5)
h = (domain(2)-domain(1))./(N-1);                                % step sizes OFF BY 1 ORGINALLY!!!! NEEDED TO BE N-1
x = cell(length(N),1);                                           % setup x scale for each N

[DDM, f, u, u_actual] = deal(cell(length(N),1)); % discrete derivative matrices and f for each N and solution matrix

for i = 1:length(N)
   DDM{i} = zeros(N(i));
   f{i} = zeros(N(i),1);
   x{i} = linspace(domain(1),domain(2),N(i));
   DDM{i}(1,1) = 1;
   DDM{i}(end,end) = cos(domain(end));
   u_actual{i} = sin(x{i});
end

for i = 1:length(N)
    for j = 2:N(i)-1
        DDM{i}(j,j-1) = 1/h(i)^2 - sin(x{i}(j))/(2*h(i));
        DDM{i}(j,j) = 1 - 2/h(i)^2;
        DDM{i}(j,j+1) = 1/h(i)^2 + sin(x{i}(j))/(2*h(i));
    end
    for j = 1:N(i)
        f{i}(j) = sin(x{i}(j))*cos(x{i}(j));
    end
end

for i = 1:length(N)
    u{i} = DDM{i}\f{i};
end

%% Plotting
fig = figure(1);
set(gcf, 'Position', window_size_for_export);

t = tiledlayout(1,20);
nexttile([1, 13]);
hold on
fplot(@(x) sin(x),domain,"Color",analytic_color,LineWidth=lwlrg)

for i = 1:length(N_plot)
    index = find(N == N_plot(i));
    plot(x{index},u{index,1},"Color",colors{i},LineWidth=lw,LineStyle=":")
end

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u(x)$', 'Interpreter', 'latex');
title('Analytic \& Numerical Solutions for $\frac{d^2 u}{dx^2} + \sin(x)\frac{du}{dx} + u(x) = f(x)$ using $u(x)=\sin(x)$', 'Interpreter', 'latex');
legend('$u(x)$' ...
              ,strcat("$u_n(x)$ for N = ",num2str(N_plot(1))) ...
              ,strcat("$u_n(x)$ for N = ",num2str(N_plot(2))) ...
              ,strcat("$u_n(x)$ for N = ",num2str(N_plot(3))) ...
              , 'Interpreter', 'latex',Location='northeast');
ylim([-1.5 1.5])
yticks(-1.5:0.25:1.5)
xticks(0:0.5:5)
ax1 = gca;
ax1.LineWidth = lwsml;
ax1.TickLabelInterpreter = 'latex';
fontsize(font,"points")
grid on
box on

%% Error calcs

[Err_p_2,Err_p_inf] = deal(zeros(length(N),1)); % setup arrays to hold errors

for j = 1:length(N)
    top_diff = u_actual{j} - u{j}'; % top of error equation
    Err_p_2(j) = norm(top_diff)/norm(u_actual{j}); %P-2 Norm Error
    Err_p_inf(j) = max(abs(top_diff))/max(abs(u_actual{j})); %P-inf Norm Error
end

nexttile([1, 7]);
loglog(N,Err_p_inf,LineWidth=lw,Color='#7E2F8E')
hold on
loglog(N,Err_p_2,LineWidth=lw,Color='#EDB120')
legend('$P=\infty$','$P=2$','Interpreter', 'latex',Location='southwest')
xlabel('$N$', 'Interpreter', 'latex');
xlim([10,100])
ylabel('Error', 'Interpreter', 'latex');
title('Relative Error Convergence', 'Interpreter', 'latex');
fontsize(font,"points")
grid on
box on
set(gca, 'TickLabelInterpreter', 'latex')
ax2 = gca;
ax2.LineWidth = lwsml;

set(gcf, 'Position', window_size_for_export);
exportgraphics(gcf,"HW0_3.jpg")
% exportgraphics(gca, 'HW0_3.pdf', 'ContentType', 'vector');