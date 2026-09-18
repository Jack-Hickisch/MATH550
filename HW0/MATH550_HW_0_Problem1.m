clc, clear, close all
window_size_for_export = [0, 50, 1920, 975];

colors = {'#d300a5', '#00a5d3', '#a5d300'}; %triadic colors
background_color = [1,1,1]; %white
analytic_color = 'k'; %black
font = 24; %fontsize
lw = 4.5; %linewidth
lwlrg = 6; %linewidth large
lwsml = 3; %linewidth small

%% Problem Parameters

domain = [0, 2*pi];                                           % interval
f = @(x) exp(sin(x));                                         % funtion
df_dx = @(x) cos(x).*exp(sin(x));                             % funtion 1st derivative (analytic)
d2f_dx2 = @(x) -sin(x).*exp(sin(x)) + cos(x).^2.*exp(sin(x)); % funtion 2nd derivative (analytic)

% Jack's Parameters
N_range = [10, 100];
N_plot = [10, 25, 100]; % choose 3 values of N to show on the plot (must be in N_range and divisible by 5)

%% Discretization and Numberical Solution

N = linspace(N_range(1),N_range(2),(N_range(2)-N_range(1))/5+1); % grid sizes
h = (domain(2)-domain(1))./N;                                    % step sizes
x = cell(length(N),1);                                           % setup x scale for each N
f_primes_n = cell(length(N),2);                                  % setup results (first column is first derivative, 2nd is 2nd)

for j = 1:length(N)
    x{j} = linspace(domain(1), domain(2), N(j)+1); % discretize the domain for each N
    [f_primes_n{j,1}, f_primes_n{j,2}] = deal(zeros(N(j)+1,1)); % setup cell arrays to store results for f'(x) and f''(x)
end

for j = 1:length(N)
    f_primes_n{j,1}(1) =   (f(x{j}(1)+h(j)) - f(x{j}(1)    )) / h(j); % forward difference
    f_primes_n{j,1}(end) = (f(x{j}(end)) - f(x{j}(end)-h(j))) / h(j); % backward difference
end

for j = 1:length(N)
    for i = 2:length(f_primes_n{j,1})-1
        f_primes_n{j,1}(i) = (f(x{j}(i)+h(j)) - f(x{j}(i)-h(j))) / (2*h(j)); % central difference (1st derivative)
    end

    for i = 1:length(f_primes_n{j,2}) %looping through entire array since using wrap around indexing
        %techincally I'm not wrapping around becuase I'm lazy but this will
        %still work and yield the same answer
        f_primes_n{j,2}(i) = (f(x{j}(i)+h(j)) - 2*f(x{j}(i)) + f(x{j}(i)-h(j))) / h(j)^2; % central difference (2nd derivative)
    end
end

%% Part 1 Plot

fig = figure(1);
set(gcf, 'Position', window_size_for_export);

t = tiledlayout(1,20);
nexttile([1, 13]); 
fplot(df_dx,domain,"Color",analytic_color,LineWidth=lwlrg)
hold on
fplot(d2f_dx2,domain,"Color",[.5,.5,.5],LineWidth=lwlrg)
for i = 1:length(N_plot)
    index = find(N == N_plot(i));
    plot(x{index},f_primes_n{index,1},"Color",colors{i},LineWidth=lw,LineStyle=":")
    plot(x{index},f_primes_n{index,2},"Color",colors{i},LineWidth=lw,LineStyle="--")
end
legend('NumColumns', 2, 'Orientation', 'horizontal');
lgd = legend('$f^{\prime}(x)$','$f^{\prime\prime}(x)$',"$f_n^{\prime}(x)$" ...
                                                      ,strcat("$f_n^{\prime\prime}(x)$ for N = ",num2str(N_plot(1))) ...
                                                      ,"$f_n^{\prime}(x)$" ...
                                                      ,strcat("$f_n^{\prime\prime}(x)$ for N = ",num2str(N_plot(2))) ...
                                                      ,"$f_n^{\prime}(x)$" ...
                                                      ,strcat("$f_n^{\prime\prime}(x)$ for N = ",num2str(N_plot(3))) ...
                                                      , 'Interpreter', 'latex',Location='southeast');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f^{(n)}(x)$', 'Interpreter', 'latex');
title('Analytic \& Numerical Solutions for Derivatives of $f(x) = e^{\sin(x)}$', 'Interpreter', 'latex');
ylim([-2.75 1.5])
yticks(-2.75:0.25:1.5)
ax1 = gca;
ax1.LineWidth = lwsml;
ax1.XTick = 0 : pi/4 : 2*pi;
ax1.TickLabelInterpreter = 'latex';
ax1.XTickLabel = {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$', '$\frac{5\pi}{4}$', '$\frac{3\pi}{2}$', '$\frac{7\pi}{4}$', '$2\pi$'};
fontsize(font,"points")
grid on
box on
set(gca, 'Color', background_color)
lgd.Color = background_color;
hText = findall(gcf, 'Type', 'text');
set(hText, 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold');

%% Part 2

[Err_p_2{1},Err_p_inf{1},Err_p_2{2},Err_p_inf{2}] = deal(zeros(length(N),1)); % setup arrays to hold errors

for j = 1:length(N)
    top_diff{1} = df_dx(x{j})' - f_primes_n{j,1}; % top of error equation
    Err_p_2{1}(j) = norm(top_diff{1})/norm(df_dx(x{j})'); %P-2 Norm Error
    Err_p_inf{1}(j) = max(abs(top_diff{1}))/max(abs(df_dx(x{j}))); %P-inf Norm Error

    top_diff{2} = d2f_dx2(x{j})' - f_primes_n{j,2}; % top of error equation
    Err_p_2{2}(j) = norm(top_diff{2})/norm(d2f_dx2(x{j})'); %P-2 Norm Error
    Err_p_inf{2}(j) = max(abs(top_diff{2}))/max(abs(d2f_dx2(x{j}))); %P-inf Norm Error
end

nexttile([1, 7]);
loglog(N,Err_p_inf{1},LineWidth=lw,Color='#7E2F8E',LineStyle=':')
hold on
loglog(N,Err_p_2{1},LineWidth=lw,Color='#EDB120',LineStyle=':')
loglog(N,Err_p_inf{2},LineWidth=lw,Color='#7E2F8E',LineStyle='--')
loglog(N,Err_p_2{2},LineWidth=lw,Color='#EDB120',LineStyle='--')

xlabel('$N$', 'Interpreter', 'latex');
ylabel('Error', 'Interpreter', 'latex');
title('Relative Error Convergence', 'Interpreter', 'latex');
fontsize(font,"points")
grid on
box on
set(gca, 'TickLabelInterpreter', 'latex')
legend('AutoUpdate', 'off');
legend('NumColumns', 2, 'Orientation', 'horizontal');
lgd2 = legend('$P=\infty$','$P=2$ for $f_n^{\prime}(x)$','$P=\infty$','$P=2$ for $f_n^{\prime\prime}(x)$','Interpreter', 'latex',Location='southwest');
ax2 = gca;
ax2.LineWidth = lwsml;
ylim([0.0005, 0.25])

exportgraphics(gcf,"HW0_1and2.jpg")
%was trying to get it to render in vecor form but no luck:
% drawnow;
% exportgraphics(fig, 'HW0_1and2.pdf', 'ContentType', 'vector','PreserveAspectRatio', 'on');
% saveas(fig, 'my_figure.pdf');