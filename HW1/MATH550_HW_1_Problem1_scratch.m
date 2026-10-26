clc, clear, close all
window_size_for_export = [0, 50, 1920, 975];

font = 24; %fontsize
lw = 4.5; %linewidth
lwlrg = 6; %linewidth large
lwsml = 3; %linewidth small

%% Problem Parameters

mu = 1; % viscosity
N_range = [15, 30]; % grid size
N = linspace(N_range(1),N_range(2),(N_range(2)-N_range(1))/5+1); % grid sizes

x_domain = [0,1];
y_domain = [0,1];

h_x = (x_domain(2)-x_domain(1))./N;
h_y = (y_domain(2)-y_domain(1))./N;

[x_grid_t1,y_grid_t1,x_grid_t2,y_grid_t2,f_b,g_b,h_b,A,first_A_row_indices,second_A_row_indices,third_A_row_indices,s_V,s_P,b,sol,U,V,P,X,Y,X2,Y2,mag_u] = deal(cell(length(N),1));

for n = 1:length(N)
    x_grid_t1{n} = linspace(x_domain(1),x_domain(2),N(n)+1);                 % x type 1 grid
    y_grid_t1{n} = linspace(y_domain(1),y_domain(2),N(n)+1);                 % y type 1 grid
    x_grid_t2{n} = linspace(x_domain(1)+h_x(n)/2,x_domain(2)-h_x(n)/2,N(n)); % x type 2 grid
    y_grid_t2{n} = linspace(y_domain(1)+h_y(n)/2,y_domain(2)-h_y(n)/2,N(n)); % y type 2 grid
    
    [f_b{n}, g_b{n}] = deal(zeros((N(n)+1)^2,1)); % forcing function vector for first and second equations, will get assembled into b vector
    h_b{n} = zeros((N(n)-1)^2,1);                 % forcing function vector for third            equation , will get assembled into b vector   
    
    A{n} = zeros(2*(N(n)+1)^2+(N(n)-1)^2); % massive matrix we are solving
    % A = spalloc(2*(N+1)^2+(N-1)^2,2*(N+1)^2+(N-1)^2,100000);
    % b = zeros(length(A),1);   % RHS matrix
    
    first_A_row_indices{n}  = num2cell(1:(N(n)+1)^2);                           % mu*L_U 0 -G_x row
    second_A_row_indices{n} = num2cell((N(n)+1)^2+1:2*(N(n)+1)^2);              % 0 mu*L_U -G_y row
    third_A_row_indices{n}  = num2cell(2*(N(n)+1)^2+1:2*(N(n)+1)^2+(N(n)-1)^2); % D_x D_y 0     row
    
    %s_U = 0            so not needed
     s_V{n} = (N(n)+1)^2;   % V shift in A array
     s_P{n} = 2*(N(n)+1)^2; % P shift in A array
end

for n = 1:length(N)
    for eqnum = 1:length(A{n})
        switch eqnum
            case first_A_row_indices{n}
                case_i = eqnum;
                [i,j] = l_i(case_i,N(n)+1);
    
                if isBoundary(i,j,N(n)+1)
                    % force setting boundary
                    A{n}(eqnum,case_i) = 1;
                    f_b{n}(case_i) = sin(2*pi*x_grid_t1{n}(i))*sin(2*pi*y_grid_t1{n}(j));
                else
                    % Laplacian on U Grid
                    A{n}(eqnum,case_i) = -2/h_x(n)^2 - 2/h_y(n)^2;
                    A{n}(eqnum,g_i(i,j+1,N(n)+1)) = 1/h_y(n)^2;
                    A{n}(eqnum,g_i(i,j-1,N(n)+1)) = 1/h_y(n)^2;
                    A{n}(eqnum,g_i(i+1,j,N(n)+1)) = 1/h_x(n)^2;
                    A{n}(eqnum,g_i(i-1,j,N(n)+1)) = 1/h_x(n)^2;
    
                    % x_comp of Gradient on P grid
                    A{n}(eqnum,s_P{n}+g_i(i  ,j,N(n))) = -1/h_x(n);
                    A{n}(eqnum,s_P{n}+g_i(i-1,j,N(n))) = 1/h_x(n);
    
                    % Forcing function
                    f_b{n}(case_i) = mu*(2*pi - 8*pi^2)*sin(2*pi*x_grid_t1{n}(i))*sin(2*pi*y_grid_t1{n}(j));
                end
    
            case second_A_row_indices{n}
                case_i = eqnum-s_V{n};
                [i,j] = l_i(case_i,N(n)+1);
                c2pix = cos(2*pi*x_grid_t1{n}(i));
                c2piy = cos(2*pi*y_grid_t1{n}(j));
    
                if isBoundary(i,j,N(n)+1)
                    % force setting boundary
                    A{n}(eqnum,s_V{n}+case_i) = 1;
                    g_b{n}(case_i) = -3.5 + c2pix*(c2piy-1);
                else
                    % Laplacian on V Grid
                    A{n}(eqnum,s_V{n}+case_i) = -2/h_x(n)^2 - 2/h_y(n)^2;
                    A{n}(eqnum,s_V{n}+g_i(i,j+1,N(n)+1)) = 1/h_y(n)^2;
                    A{n}(eqnum,s_V{n}+g_i(i,j-1,N(n)+1)) = 1/h_y(n)^2;
                    A{n}(eqnum,s_V{n}+g_i(i+1,j,N(n)+1)) = 1/h_x(n)^2;
                    A{n}(eqnum,s_V{n}+g_i(i-1,j,N(n)+1)) = 1/h_x(n)^2;
    
                    % y_comp of Gradient on P grid
                    A{n}(eqnum,s_P{n}+g_i(i,j,N(n))) = -1/h_x(n);
                    A{n}(eqnum,s_P{n}+g_i(i,j-1,N(n))) = 1/h_x(n);
    
                    % Forcing function
                    g_b{n}(case_i) = 4*mu*pi^2*c2pix-2*pi*c2pix*c2piy;
                end
    
            case third_A_row_indices{n}
                case_i = eqnum-s_P{n};
                [i,j] = l_i_nb(case_i,N(n)+1); % on type 1 grid but exluding boundaries since not important here
    
                if isBoundary(i,j,N(n)+1)
                    % no more boundaries to set :)
                else
                    % x_comp of Divergence on U grid
                    A{n}(eqnum,g_i(i+1,j,N(n)+1)) =  1/(2*h_x(n));
                    A{n}(eqnum,g_i(i,j,N(n)+1)) = -1/(2*h_x(n));
    
                    % y_comp of Divergence on V grid
                    A{n}(eqnum,s_V{n}+g_i(i,j+1,N(n)+1)) =  1/(2*h_y(n));
                    A{n}(eqnum,s_V{n}+g_i(i,j,N(n)+1)) = -1/(2*h_y(n));
    
                    % Forcing function
                    %h_b = zero for all since divergence of velocity vector field is zero so no need to change it 
                end
        end
    end

    b{n} = [f_b{n};g_b{n};h_b{n}]; % RHS forcing vector

    sol{n} = A{n}\b{n};
    
    U{n} = reshape(sol{n}(1:(N(n)+1)^2),N(n)+1,N(n)+1);
    V{n} = reshape(sol{n}((N(n)+1)^2+1:2*(N(n)+1)^2),N(n)+1,N(n)+1);
    P{n} = reshape(sol{n}(2*(N(n)+1)^2+1:end),N(n),N(n));
    
    [X{n},Y{n}] = meshgrid(x_grid_t1{n},y_grid_t1{n});
    [X2{n},Y2{n}] = meshgrid(x_grid_t2{n},y_grid_t2{n});
    
    mag_u{n} = sqrt(U{n}.^2 + V{n}.^2);
end

figure(1)
surf(X{end},Y{end},U{end});

figure(2)
surf(X{end},Y{end},V{end});

figure(3)
surf(X2{end},Y2{end},P{end});

figure(4)
surf(Y{end},X{end},mag_u{end});
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');







function [i,j] = l_i(g_i, N) %local index funciton, accepts a global index and NxN matrix size; returns i and j
    j = floor((g_i-1)/N)+1;
    i = mod(g_i-1,N)+1;
end

function g_i = g_i(i,j,N) %global index function, accepts i and j and NxN matrix size; returns global index position
    g_i = i + (j-1)*N;
end

function bound = isBoundary(i,j,N)
    bound = i == 1 || j == 1 || i == N || j == N;
end

function [i,j] = l_i_nb(g_i,N) % EXCLUDES BOUNDARIES! local index funciton, accepts a global index and NxN matrix size; returns i and j
    j = floor((g_i-1)/(N-2))+2;
    i = mod(g_i-1,N-2)+2;
end







% look at the 