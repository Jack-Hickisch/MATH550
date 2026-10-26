clc, clear, close all
window_size_for_export = [0, 50, 1920, 975];
font = 24; 
lw = 4;
lwlrg = 2; 
lwsml = 1; 
%% Problem Parameters
mu = 1; % viscosity
N_range = [10, 100]; % grid size
N = linspace(N_range(1),N_range(2),(N_range(2)-N_range(1))/15+1); % grid sizes
x_domain = [0,1];
y_domain = [0,1];
h_x = (x_domain(2)-x_domain(1))./N;
h_y = (y_domain(2)-y_domain(1))./N;
global Ugd Vgd Pgd
[x_grid_t1,y_grid_t1,x_grid_t2,y_grid_t2,y_grid_t2_SHORT,x_grid_t2_LONG,Ugd,Vgd,Pgd,Cgd,Upnts,Vpnts,Ppnts,Tpnts,...
SU,SV,CM,f_b,g_b,h_b,A,b,sol,U,V,P,PC,X,Y,XU,YU,XV,YV,XP,YP,Ua,Va,Pa,Ua_shift,Va_shift] = deal(cell(length(N),1));
interpolation_type = {'linear','cubic','spline'};
[UC,VC,mag_u] = deal(cell(length(N),length(interpolation_type)));
%% Matrix Setup and Solving
for n = 1:length(N)
    Ugd{n} = [N(n)+1,N(n)+2]; 
    Vgd{n} = [N(n)  ,N(n)+1]; 
    Pgd{n} = [N(n)  ,N(n)  ]; 
    Cgd{n} = [N(n)+1,N(n)+1]; 
    x_grid_t1{n} =       linspace(x_domain(1)         ,x_domain(2)         ,Cgd{n}(1)  ); 
    y_grid_t1{n} =       linspace(y_domain(1)         ,y_domain(2)         ,Cgd{n}(2)  ); 
    x_grid_t2{n} =       linspace(x_domain(1)+h_x(n)/2,x_domain(2)-h_x(n)/2,Cgd{n}(1)-1); 
    x_grid_t2_LONG{n} =  linspace(x_domain(1)-h_x(n)/2,x_domain(2)+h_x(n)/2,Cgd{n}(1)+1); 
    y_grid_t2{n} =       linspace(y_domain(1)-h_y(n)/2,y_domain(2)+h_y(n)/2,Ugd{n}(2)  ); 
    y_grid_t2_SHORT{n} = linspace(y_domain(1)+h_y(n)/2,y_domain(2)-h_y(n)/2,Ugd{n}(2)-2); 
    Upnts{n} = prod(Ugd{n});                   
    Vpnts{n} = prod(Vgd{n});                   
    Ppnts{n} = prod(Pgd{n});                   
    Tpnts{n} = Upnts{n} + Vpnts{n} + Ppnts{n}; 
    SU{n} = sparse(Upnts{n},Tpnts{n}); 
    SV{n} = sparse(Vpnts{n},Tpnts{n}); 
    CM{n} = sparse(Ppnts{n},Tpnts{n}); 
    
    f_b{n} = zeros(Upnts{n},1); 
    g_b{n} = zeros(Vpnts{n},1); 
    h_b{n} = zeros(Ppnts{n},1); 
end
for n = 1:length(N)
    for j = 1:Ugd{n}(2)
        for i = 1:Ugd{n}(1)
            
            k = lcl_i(i,j,Ugd{n});
            if j == 1             % BTM Boundary - ghost nodes
                SU{n}(k,g_i(i  ,j  ,'U',n)) =  .5;
                SU{n}(k,g_i(i  ,j+1,'U',n)) =  .5;
                f_b{n}(k) = 0;
            elseif j == Ugd{n}(2) % TOP Boundary - ghost nodes
                SU{n}(k,g_i(i  ,j  ,'U',n)) =  .5;
                SU{n}(k,g_i(i  ,j-1,'U',n)) =  .5;
                f_b{n}(k) = 0;
            elseif i == 1         % LHS Boundary - Periodic
                SU{n}(k,g_i(i  ,j  ,'U',n)) = mu*(-2/h_x(n)^2 - 2/h_y(n)^2);
                SU{n}(k,g_i(i+1,j  ,'U',n)) = mu*( 1/h_x(n)^2);
                SU{n}(k,g_i(Ugd{n}(1)-1,j  ,'U',n)) = mu*( 1/h_x(n)^2); 
                SU{n}(k,g_i(i  ,j+1,'U',n)) = mu*( 1/h_y(n)^2);
                SU{n}(k,g_i(i  ,j-1,'U',n)) = mu*( 1/h_y(n)^2);
                % x_comp of Gradient on P grid wrapping to P_N
                SU{n}(k,g_i(i  ,j-1,'P',n)) = -1/h_x(n);
                SU{n}(k,g_i(Pgd{n}(1),j-1,'P',n)) =  1/h_x(n);
                % Forcing function
                f_b{n}(k) = (2*pi - 8*mu*pi^2) * sin(2*pi*x_grid_t1{n}(i))*sin(2*pi*y_grid_t2{n}(j));
            elseif i == Ugd{n}(1) % RHS Boundary - Copy left bound (periodic)
                SU{n}(k,g_i(i  ,j  ,'U',n)) =  1;
                SU{n}(k,g_i(1  ,j  ,'U',n)) = -1;
                f_b{n}(k) = 0;
            else                  % INT (interior point)
                SU{n}(k,g_i(i  ,j  ,'U',n)) = mu*(-2/h_x(n)^2 - 2/h_y(n)^2);
                SU{n}(k,g_i(i+1,j  ,'U',n)) = mu*( 1/h_x(n)^2);
                SU{n}(k,g_i(i-1,j  ,'U',n)) = mu*( 1/h_x(n)^2);
                SU{n}(k,g_i(i  ,j+1,'U',n)) = mu*( 1/h_y(n)^2);
                SU{n}(k,g_i(i  ,j-1,'U',n)) = mu*( 1/h_y(n)^2);
                % x_comp of Gradient on P grid
                SU{n}(k,g_i(i  ,j-1,'P',n)) = -1/h_x(n);
                SU{n}(k,g_i(i-1,j-1,'P',n)) =  1/h_x(n);
                % Forcing function
                f_b{n}(k) = (2*pi - 8*mu*pi^2) * sin(2*pi*x_grid_t1{n}(i))*sin(2*pi*y_grid_t2{n}(j));
            end
        end
    end
    for j = 1:Vgd{n}(2)
        for i = 1:Vgd{n}(1)
            c2pix = cos(2*pi*x_grid_t2{n}(i));
            c2piy = cos(2*pi*y_grid_t1{n}(j));
            k = lcl_i(i,j,Vgd{n}); 
            if j == 1             % BTM Boundary - Direcelet
                SV{n}(k,g_i(i  ,j  ,'V',n)) =  1;
                g_b{n}(k) = -3.5;
            elseif j == Vgd{n}(2) % TOP Boundary - Direcelet
                SV{n}(k,g_i(i  ,j  ,'V',n)) =  1;
                g_b{n}(k) = -3.5;
            elseif     i == 1         % LHS Boundary - periodic
                imx = Vgd{n}(1);                      
                SV{n}(k,g_i(i  ,j  ,'V',n)) = mu*(-2/h_x(n)^2 - 2/h_y(n)^2);
                SV{n}(k,g_i(i+1,j  ,'V',n)) = mu*( 1/h_x(n)^2);
                SV{n}(k,g_i(imx,j  ,'V',n)) = mu*( 1/h_x(n)^2);
                SV{n}(k,g_i(i  ,j+1,'V',n)) = mu*( 1/h_y(n)^2);
                SV{n}(k,g_i(i  ,j-1,'V',n)) = mu*( 1/h_y(n)^2);
                % y_comp of Gradient on P grid
                SV{n}(k,g_i(i  ,j  ,'P',n)) = -1/h_y(n);
                SV{n}(k,g_i(i  ,j-1,'P',n)) =  1/h_y(n);
                % Forcing function
                g_b{n}(k) = mu*(4*pi^2*c2pix-8*pi^2*c2pix*c2piy)-2*pi*c2pix*c2piy;
            elseif i == Vgd{n}(1) % RHS Boundary - periodic
                SV{n}(k,g_i(i  ,j  ,'V',n)) = mu*(-2/h_x(n)^2 - 2/h_y(n)^2);
                SV{n}(k,g_i(1  ,j  ,'V',n)) = mu*( 1/h_x(n)^2);
                SV{n}(k,g_i(i-1,j  ,'V',n)) = mu*( 1/h_x(n)^2);
                SV{n}(k,g_i(i  ,j+1,'V',n)) = mu*( 1/h_y(n)^2);
                SV{n}(k,g_i(i  ,j-1,'V',n)) = mu*( 1/h_y(n)^2);
                % y_comp of Gradient on P grid
                SV{n}(k,g_i(i  ,j  ,'P',n)) = -1/h_y(n);
                SV{n}(k,g_i(i  ,j-1,'P',n)) =  1/h_y(n);
                % Forcing function
                g_b{n}(k) = mu*(4*pi^2*c2pix - 8*pi^2*c2pix*c2piy)-2*pi*c2pix*c2piy;
            else                  % INT (interior point)
                SV{n}(k,g_i(i  ,j  ,'V',n)) = mu*(-2/h_x(n)^2 - 2/h_y(n)^2);
                SV{n}(k,g_i(i+1,j  ,'V',n)) = mu*( 1/h_x(n)^2);
                SV{n}(k,g_i(i-1,j  ,'V',n)) = mu*( 1/h_x(n)^2);
                SV{n}(k,g_i(i  ,j+1,'V',n)) = mu*( 1/h_y(n)^2);
                SV{n}(k,g_i(i  ,j-1,'V',n)) = mu*( 1/h_y(n)^2);
                % y_comp of Gradient on P grid
                SV{n}(k,g_i(i  ,j  ,'P',n)) = -1/h_y(n);
                SV{n}(k,g_i(i  ,j-1,'P',n)) =  1/h_y(n);
                % Forcing function
                g_b{n}(k) = mu*(4*pi^2*c2pix-8*pi^2*c2pix*c2piy)-2*pi*c2pix*c2piy;
            end
        end
    end
    for j = 1:Pgd{n}(2)
        for i = 1:Pgd{n}(1)
                k = lcl_i(i,j,Pgd{n}); 
                CM{n}(k,g_i(i  ,j+1,'U',n)) = -1/h_x(n);
                CM{n}(k,g_i(i+1,j+1,'U',n)) =  1/h_x(n);
                CM{n}(k,g_i(i  ,j+1,'V',n)) =  1/h_y(n);
                CM{n}(k,g_i(i  ,j  ,'V',n)) = -1/h_y(n);
        end
    end
    A{n} = [SU{n};SV{n};CM{n}];
    b{n} = [f_b{n};g_b{n};h_b{n}];
    
    col_lambda = [zeros(Upnts{n}+Vpnts{n}, 1); ones(Ppnts{n}, 1)];
    row_lambda = [zeros(1, Upnts{n}+Vpnts{n}), ones(1, Ppnts{n}), 0];
    
    A_ext = [A{n}, col_lambda; row_lambda];
    b_ext = [b{n}; 0];
    sol_ext = A_ext \ b_ext;
    sol{n} = sol_ext(1:end-1); 
    U{n} = reshape(sol{n}(1                  :Upnts{n}         ),Ugd{n}(1),Ugd{n}(2))';
    V{n} = reshape(sol{n}(Upnts{n}+1         :Upnts{n}+Vpnts{n}),Vgd{n}(1),Vgd{n}(2))';
    P{n} = reshape(sol{n}(Upnts{n}+Vpnts{n}+1:end              ),Pgd{n}(1),Pgd{n}(2))';
    
    [X{n} ,Y{n} ] = meshgrid(x_grid_t1{n},y_grid_t1{n});
    [XU{n},YU{n}] = meshgrid(x_grid_t1{n},y_grid_t2{n});
    [XV{n},YV{n}] = meshgrid(x_grid_t2{n},y_grid_t1{n});
    [XV{n},YV{n}] = meshgrid(x_grid_t2_LONG{n},y_grid_t1{n});
    [XP{n},YP{n}] = meshgrid(x_grid_t2{n},y_grid_t2_SHORT{n});
    for type = 1:length(interpolation_type)
        UC{n,type} = interp2(XU{n},YU{n},U{n}                        ,X{n},Y{n},interpolation_type{type});
        VC{n,type} = interp2(XV{n},YV{n},[V{n}(:,end),V{n},V{n}(:,1)],X{n},Y{n},interpolation_type{type});
        mag_u{n,type} = sqrt(UC{n,type}.^2 + VC{n,type}.^2);
    end
    PC{n} = P{n} - mean(P{n},'all');
    %% Error Calcs
    for j = 1:Cgd{n}(2)
        for i = 1:Cgd{n}(1)
            Ua{n}(i,j) = sin(2*pi*y_grid_t1{n}(j))*sin(2*pi*x_grid_t1{n}(i));            
            Va{n}(i,j) = -3.5 + cos(2*pi*x_grid_t1{n}(i))*(cos(2*pi*y_grid_t1{n}(j))-1); 
        end
    end
    for j = 1:Ugd{n}(2)-2
        for i = 1:Ugd{n}(1)
            Ua_shift{n}(i,j) = sin(2*pi*y_grid_t2_SHORT{n}(j))*sin(2*pi*x_grid_t1{n}(i));      
        end
    end
    for j = 1:Vgd{n}(2)
        for i = 1:Vgd{n}(1)
            Va_shift{n}(i,j) = -3.5 + cos(2*pi*x_grid_t2{n}(i))*(cos(2*pi*y_grid_t1{n}(j))-1);          
        end
    end
    for j = 1:Pgd{n}(2)
        for i = 1:Pgd{n}(1)
            Pa{n}(i,j) = sin(2*pi*y_grid_t2_SHORT{n}(j))*cos(2*pi*x_grid_t2{n}(i));      
        end
    end
    Ua{n} = Ua{n}';
    Va{n} = Va{n}';
    Ua_shift{n} = Ua_shift{n}';
    Va_shift{n} = Va_shift{n}';
    Pa{n} = Pa{n}';
end
%% Errors Continued
UErr_p_2 = deal(zeros(length(N),1+length(interpolation_type))); 
VErr_p_2 = deal(zeros(length(N),1+length(interpolation_type))); 
PErr_p_2 = deal(zeros(length(N),1+length(interpolation_type))); 
for n = 1:length(N)
    Utop_diff = Ua_shift{n} - U{n}(2:end-1,:); 
    Utop_diff = Utop_diff(:); 
    UErr_p_2(n,1) = norm(Utop_diff)/norm(Ua_shift{n}(:)); 
    Vtop_diff = Va_shift{n} - V{n}; 
    Vtop_diff = Vtop_diff(:); 
    VErr_p_2(n,1) = norm(Vtop_diff)/norm(Va_shift{n}(:)); 
    
    Ptop_diff = (Pa{n} - mean(Pa{n},'all')) - (P{n} - mean(P{n},'all')); 
    Ptop_diff = Ptop_diff(:); 
    PErr_p_2(n,1) = norm(Ptop_diff)/norm(Pa{n}(:)); 
    for type = 2:length(interpolation_type)+1
        Utop_diff = Ua{n} - UC{n, type-1}; 
        Utop_diff = Utop_diff(:); 
        UErr_p_2(n,type) = norm(Utop_diff)/norm(Ua{n}(:)); 
    
        Vtop_diff = Va{n} - VC{n, type-1}; 
        Vtop_diff = Vtop_diff(:); 
        VErr_p_2(n,type) = norm(Vtop_diff)/norm(Va{n}(:)); 
    end
end
%% Plotting
figure(1)
t = tiledlayout(1,2, 'TileSpacing', 'compact'); 

% Flow Field
nexttile; 
hold on
pcolor(X{end},Y{end},mag_u{end,1}) 
shading interp;
try colormap(slanCM('bubblegum')); catch, colormap('parula'); end
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$||\vec{u}||$',Interpreter='latex')
c = colorbar;
set(c, 'TickLabelInterpreter', 'latex');
stream = streamline(X{end},Y{end},UC{end,1},VC{end,1},x_grid_t2{end},1*ones(length(x_grid_t2{end}),1));
set(stream, 'Color', [1 1 1], 'LineWidth', lwlrg);
axis square
fontsize(font,"points")
set(gca, 'TickLabelInterpreter', 'latex')
box on
grid on
xticks(0:.2:1)
yticks(0:.2:1)
ax1 = gca;
ax1.LineWidth = lwsml;

% Combined Error Plot
nexttile;
loglog(N,UErr_p_2(:,1), '-o', 'LineWidth',lw,'Color','#7E2F8E')
hold on
loglog(N,VErr_p_2(:,1), '-s', 'LineWidth',lw,'Color','#EDB120')
loglog(N,PErr_p_2(:,1), '-d', 'LineWidth',lw,'Color',[0,0,0])
loglog(N, UErr_p_2(:,2), '-x', 'LineWidth',lw, 'Color', '#0072BD') 
loglog(N, UErr_p_2(:,3), '-+', 'LineWidth',lw, 'Color', '#D95319') 
loglog(N, UErr_p_2(:,4), '-*', 'LineWidth',lw, 'Color', '#77AC30') 
ref_y = UErr_p_2(1,1) * (N(1)./N).^2;
loglog(N, ref_y, 'k--', 'LineWidth', 2);
legend('$U$ Staggered','$V$ Staggered','$P$ Staggered',...
       '$U$ Linear mapped','$U$ Cubic mapped', '$U$ Spline mapped', ...
       '$\mathcal{O}(N^{-2})$', 'Interpreter', 'latex',Location='southwest')
xlabel('$N$', 'Interpreter', 'latex');
ylabel('Relative Error ($L_2$)', 'Interpreter', 'latex');
title('Convergence Rates', 'Interpreter', 'latex');
axis square
fontsize(font,"points")
grid on
box on
set(gca, 'TickLabelInterpreter', 'latex')
ax2 = gca;
ax2.LineWidth = lwsml;

set(gcf, 'Position', window_size_for_export);

%% Helper Functions
function k = lcl_i(i,j,Xgd) 
    k = i + (j-1)*Xgd(1);
end
function g_i = g_i(i,j,grid,n) 
global Ugd Vgd Pgd
    switch grid
        case 'U'
            g_i = 0                     + i + (j-1)*Ugd{n}(1);
        case 'V'
            g_i = prod(Ugd{n})          + i + (j-1)*Vgd{n}(1);
        case 'P'
            g_i = prod(Ugd{n}) + prod(Vgd{n}) + i + (j-1)*Pgd{n}(1);
    end
end