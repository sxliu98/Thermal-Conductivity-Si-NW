clear clc
set(0,'defaultAxesFontSize',12)
set(0,'defaultLineLineWidth',2)
set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder',':|-.|--|-')

h = 6.626e-34 / (2*pi); % Planck's constant divided by 2𝜋, J*s
k_B = 1.38e-23; % Boltzmann constant, J/K
a = 5.4310e-10; % Lattice constant Si, m 
eV = 1.602e-19; % The energy of a single electron, J
k_max = 2*pi/a; % Maximum wave vector, m^-1 
xX = 0:0.001:1;  % x = k/k_max
kX = xX*k_max;
TX = 10:500;
TX = unique(TX(:));
T_s = 300;

%% Model Holland
Cph = @(x) k_B*x.^2.*exp(x)./(exp(x)-1).^2; 
vT = 5860; vT1 = 2000;
vL = 8480; vL1 = 4240;
theta1 = 180; theta2 = 210;
theta3 = 570; theta4 = 350;
thetaD = 643;
Bl = 2e-24; Bt = 9.3e-13; Btu = 5.5e-18;
A = 1.32e-45;
tau1I = @(omega) A*omega.^4;  alpha = A*(k_B/h)^4;    tau1i = @(x,T) alpha*x.^4.*T.^4;
tau1LN = @(omega,T) Bl*omega.^2.*T.^3;  betaL = Bl*(k_B/h)^2;    tau1ln = @(x,T) betaL*x.^2.*T.^5;
tau1TN = @(omega,T) Bt*omega.*T.^4;  betaT = Bt*(k_B/h)^1;    tau1tn = @(x,T) betaT*x.*T.^5;
tau1TU = @(omega,T) Btu*omega.^2./sinh(hm*omega./T/k_B); betaTU = Btu*(k_B/h)^2;    tau1tu = @(x,T) betaTU*x.^2.*T.^2./sinh(x);

tauT01 = @(x,T) tau1i(x,T)+tau1tn(x,T); 
tauTU1 = @(x,T) tau1i(x,T)+tau1tu(x,T);
tauL1 = @(x,T) tau1i(x,T)+tau1ln(x,T);
Const = @(v) 1/(2*pi^2*v)*(k_B/h)^3;

for iT = 1:length(TX)    
    T1 = TX(iT);
    f = @(x) x.^2.*Cph(x)./tauT01(x,T1);
        kT0(iT) = 2/3*Const(vT)*T1^3*integral(f,0,theta1/T1);
    f = @(x) x.^2.*Cph(x)./tauTU1(x,T1);
        kTU(iT) = 2/3*Const(vT1)*T1^3*integral(f,theta1/T1,theta2/T1);
    f = @(x) x.^2.*Cph(x)./tauL1(x,T1);
        kL1(iT) = 1/3*Const(vL)*T1^3*integral(f,0,theta4/T1);
        kL2(iT) = 1/3*Const(vL1)*T1^3*integral(f,theta4/T1,theta3/T1);
    f = @(x) x.^2.*Cph(x);
        Cp(iT) = 3*(T1/thetaD)^3*integral(f,0,thetaD/T1);        
end
k_bulk1 = kT0+kTU+kL1+kL2;

%% Impurity
A_imp = 1.32e-45;
TAU_imp = @(omega) A_imp*omega.^4;

%% Model Ward-Broido
coeff = [0 6.759 -3.353 -2.288 1.702; 0 9.756 -0.371 -1.644 0];
OMEGA_bar = @(x) (coeff(:,5).*x.^4+coeff(:,4).*x.^3+coeff(:,3).*x.^2+coeff(:,2).*x+coeff(:,1))*10^13;
Vg_bar = @(x) abs((4*coeff(:,5).*x.^3 + 3*coeff(:,4).*x.^2 + 2*coeff(:,3).*x + coeff(:,2))*10^13./k_max);
Vp_bar = @(x) OMEGA_bar(x)./(x*k_max);
omegaX_bar = OMEGA_bar(xX);
vgX_bar = Vg_bar(xX);
vpX_bar = Vp_bar(xX);
temp_D = 643; 
A = [253322; 163921; 2012; 507].*[(h/eV*10^3)^2; (h/eV*10^3)^2; (h/eV*10^3)^4; (h/eV*10^3)^4];
TAU_N_WB = @(omega,T) A([1 2]).*omega.^2*T*(1-exp(-3*T/temp_D));
TAU_U_WB = @(omega,T) A([3 4]).*omega.^4*T*(1-exp(-3*T/temp_D));
TAU_WB = @(omega,T) TAU_N_WB(omega,T) + TAU_U_WB(omega,T);

%%
k_bulk_mode2 = zeros(2,length(TX));
k_bulk2 = zeros(1,length(TX));
for mode = 1:2
    omegaM2 = omegaX_bar(mode,:);
    vg = vgX_bar(mode,:);
    vp = vpX_bar(mode,:); 
    DX2 = 1/(2*pi^2)*omegaM2.^2./(vg.*vp.^2);    
    DX2(1) = DX2(2); % приближение
    for iT = 1:length(TX)        
        T = TX(iT);
        XX2 = h*omegaM2./(k_B*T);        
        CphX2 = Cph(XX2);        
        if isnan(CphX2(1))
            CphX2(1) = CphX2(2); % приближение
        end
        tauX2 = (TAU_WB(omegaM2,T)+TAU_imp(omegaM2)).^-1;
        tauXX2 = tauX2(mode,:);
        tauXX2(1) = 0;        
        y2 = CphX2.*DX2.*vg.^2.*tauXX2;        
        k_bulk_mode2(mode,iT) = trapz(omegaM2,y2);       
    end    
end
k_bulk2(:) = 1/3*(2*k_bulk_mode2(1,:)+k_bulk_mode2(2,:));

%% Model Slack
v0 = [5840; 8430];
v_s = 6400;
theta_slack = [240; 586];
theta_D = 643;
omega_D = theta_D*k_B/h;
omega_slack = [5.1 12.4]*2*pi*1e12;
omegaX_slack = [linspace(0,omega_slack(1),length(kX)); linspace(0,omega_slack(2),length(kX))];
beta = 1/omega_slack(1)^2*(omega_D*v0(1)/v_s/omega_slack(1)-1);
alpha = 1/omega_slack(2)*(omega_D*v0(2)/v_s/omega_slack(2)-1);
kX1_slack = omegaX_slack(1,:)./v0(1).*(1+beta.*omegaX_slack(1,:).^2);
kX2_slack = omegaX_slack(2,:)./v0(2).*(1+alpha.*omegaX_slack(2,:));
vpX_slack = [omegaX_slack(1,:)./kX1_slack; omegaX_slack(2,:)./kX2_slack];
%%% Herring
B_N = [7.1e-13; 2.4e-24];
TAU_NT_Her = @(omega,T) B_N(1).*omega.*T.^3.3;
TAU_NL_Her = @(omega,T) B_N(2).*omega.^2.*T.^3;
TAU_N_Her = @(omega,T) [TAU_NT_Her(omega,T); TAU_NL_Her(omega,T)];
%%% Slack
B_slack_U = [10; 5.5]*10^-20;
Theta_slack = [240; 586];
TAU_U_slack = @(omega,T) B_slack_U.*omega.^2*T.*exp(-Theta_slack./(3*T));

%%
k_bulk_mode3 = zeros(2,length(TX));
k_bulk3 = zeros(1,length(TX));
for mode = 1:2
    omegaM3 = omegaX_slack(mode,:);
    vg = vgX_bar(mode,:);
    vp = vpX_bar(mode,:); 
    DX3 = 1/(2*pi^2)*omegaM3.^2./(vg.*vp.^2);    
    DX3(1) = DX3(2); % approximation
    for iT = 1:length(TX)        
        T = TX(iT);
        XX3 = h*omegaM3./(k_B*T);        
        CphX3 = Cph(XX3);        
        if isnan(CphX3(1))
            CphX3(1) = CphX3(2); % approximation
        end
        tauX3 = (TAU_N_Her(omegaM3,T)+TAU_U_slack(omegaM3,T)+TAU_imp(omegaM3)).^-1;
        tauXX3 = tauX3(mode,:);
        tauXX3(1) = 0;        
        y3 = CphX3.*DX3.*vg.^2.*tauXX3;        
        k_bulk_mode3(mode,iT) = trapz(omegaM3,y3);       
    end    
end
k_bulk3(:) = 1/3*(2*k_bulk_mode3(1,:)+k_bulk_mode3(2,:));

%%
k_bulk_mode4 = zeros(2,length(TX));
k_bulk4 = zeros(1,length(TX));
for mode = 1:2
    omegaM4 = omegaX_slack(mode,:);
    vX4 = vpX_slack(mode,:);
    vX4(1) = v0(mode);
    DX4 = 1/(2*pi^2)*omegaM4.^2./(vX4.^3);    
    DX4(1) = DX4(2); % approximation
    for iT = 1:length(TX)        
        T = TX(iT);
        XX4 = h*omegaM4./(k_B*T);        
        CphX4 = Cph(XX4);        
        if isnan(CphX4(1))
            CphX4(1) = CphX4(2); % approximation
        end
        tauX4 = (TAU_N_Her(omegaM4,T)+TAU_U_slack(omegaM4,T)+TAU_imp(omegaM4)).^-1;
        tauXX4 = tauX4(mode,:);
        tauXX4(1) = 0;        
        y4 = CphX4.*DX4.*vX4.^2.*tauXX4;        
        k_bulk_mode4(mode,iT) = trapz(omegaM4,y4);       
    end    
end
k_bulk4(:) = 1/3*(2*k_bulk_mode4(1,:)+k_bulk_mode4(2,:));

%% Plot

T_bulk_expholl = [1.6249	1.9358	2.9751	4.9511	7.43   10.017	14.977	19.804	24.839	34.267	49.358	59.23	68.65	78.88	89.85	99.7	124.9	148.53	175.15	199.5	247.84	297.43];
k_bulk_expholl = [27.861	43.425	146.89	556.22	1418.2 1988.3	3589	4575.7	4653.4	3806.4	2616.1	2122.1	1691.9	1384.6	1153	951.9	610.5	419.75	326.02	266.81	194.96	158.15];
T_bulk_expslack =[396.08	496.42	595.7	696.5	786.4	903.6	1000];
k_bulk_expslack =[104.99	80.13	63.32	52.273	42.78	35.316	31.5];

figure('OuterPosition',[0 0 800 600])
plot(TX,k_bulk1,'k:')
hold on
plot(TX,k_bulk4,'k-')
hold on
plot(TX,k_bulk2,'k--')
hold on
plot(TX,k_bulk3,'k-.')
hold on
plot(T_bulk_expholl,k_bulk_expholl,'o','Color','r');
hold on
plot(T_bulk_expslack,k_bulk_expslack,'s','Color','r','MarkerSize',8);

set(gca,'linewidth',1.5,'FontSize',16)
xlim([50 500])
grid on
ylim([0 2000])
legend('Model Holland','Model Slack','Model Ward Broido','Experiment Holland','Experiment Slack')
legend boxoff
xlabel("$T$, K",'Interpreter','latex','FontSize',22)
ylabel('$\kappa_{bulk}$, W/(mK)','Interpreter','latex','FontSize',22)
%%
T_bulk_exp = [T_bulk_expholl(12:end) T_bulk_expslack(1:2)];
k_bulk_exp = [k_bulk_expholl(12:end) k_bulk_expslack(1:2)];

k_bulk_the1 = interp1(TX,k_bulk1,T_bulk_exp);
k_bulk_the4 = interp1(TX,k_bulk4,T_bulk_exp);
k_bulk_the2 = interp1(TX,k_bulk2,T_bulk_exp);

error1 = abs(k_bulk_exp-k_bulk_the1)./k_bulk_exp*100;
error4 = abs(k_bulk_exp-k_bulk_the4)./k_bulk_exp*100;
error2 = abs(k_bulk_exp-k_bulk_the2)./k_bulk_exp*100;

figure('OuterPosition',[0 0 800 600])
plot(T_bulk_exp,error1,':o')
hold on
plot(T_bulk_exp,error4,'-s')
hold on
plot(T_bulk_exp,error2,'--*')

xlim([50 500])
set(gca,'linewidth',1.5,'FontSize',16)
grid on
legend('Model Holland','Model Slack','Model Ward Broido','Location','best')
legend boxoff
xlabel("$T$, K",'Interpreter','latex','FontSize',22)
ylabel('$\frac{(\kappa_{bulk}^{model}-\kappa_{bulk}^{exp})}{\kappa_{bulk}^{exp}}\cdot100\%$','Interpreter','latex','FontSize',22)

%% Size Effect for wire

pX = 0:0.05:1;
Kn_1 = (1:1:10)'*(10.^(-8:8));
Kn = unique(Kn_1(:));
m = 1:50;
S = zeros(length(m),length(Kn));
F = zeros(length(pX),length(Kn));

for i = 1:length(m)
    for j = 1:length(Kn)
        f = @(r,J) sqrt(1-J.^2).*exp(-m(i)*J./(r*Kn(j))).*(-sqrt(r.^2-r.^4));
        S(i,j) = integral2(f,1,0,0,1,'RelTol',1e-20);
    end
end

for k = 1:length(pX)
    Y = S.*(m'.*pX(k).^(m'-1));
    F(k,:) = 1-(12*((1-pX(k))^2/pi))*sum(Y,1);
end

%%%%%% Plot %%%%%%
figure('OuterPosition',[0 0 800 600]);
p_plotX = [0.2 0.4 0.6 0.8];
for ip = 1:length(p_plotX)
    p_plot = p_plotX(ip);
    ind_p = find( abs(pX-p_plot)<1e-15 );
    semilogx(Kn,F(ind_p,:),'LineWidth',2)
    hold on
end
hold off
set(gca,'XDir','reverse')
set(gca,'linewidth',1.5,'FontSize',16)
grid on
xlim([1e-2 1e2])
legend('p=0.2','p=0.4','p=0.6','p=0.8','Location','best')
legend boxoff
xlabel("Kn",'Interpreter','latex','FontSize',22)
ylabel("$ F $",'Interpreter','latex','FontSize',22)

%% Roughness

omegaX = omegaX_slack;
rmsX = (1:10)'*10.^(-11:-8);
rms = unique(rmsX(:));

for i = 11:37
    rms(i) = rms(i+1);
end
rms = rms(1:37);

kX = xX*k_max;
theta = linspace(0,pi/2,100);

p = zeros(length(kX),length(theta),length(rms));
for ir = 1:length(rms)
    p(:,:,ir) = exp(-4*(kX.^2)'*rms(ir).^2*(cos(theta)).^2);
end

omegaX(:,1) = omegaX(:,2);

f0 = zeros(2,length(omegaX),length(TX));
for iT = 1:length(TX)
    f0(:,:,iT) = 1./(exp(h*omegaX/k_B/TX(iT))-1);
end

DX = zeros(4,length(kX));
for i = 1:2
    DX(i,:) = 1/(2*pi^2)*omegaX(i,:).^2./(vpX_slack(i,:).^3);
end
DX(:,1)=DX(:,2);

p_ave = zeros(2,length(TX),length(rms));
for i = 1:2
    for ir = 1:length(rms)
        for iT = 1:length(TX)
            p_ave(i,iT,ir) = sum( (f0(i,:,iT))'.*(DX(i,:))'.*p(:,:,ir).*sin(theta),[1 2] ) ./ (sum((f0(i,:,iT))'.*(DX(i,:))' * sin(theta), [1 2]));
        end
    end
end

%%%%%%%%%$$ Plot %%%%%%%%%$$
figure('OuterPosition',[0 0 800 600]);

rms_plotX = 10.^[-11 -10 -9 -8];

t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

nexttile
for irms = 1:length(rms_plotX)
    rms_plot = rms_plotX(irms);
    ind_r = find(abs(rms-rms_plot)<1e-20);
    semilogx(TX,p_ave(1,:,ind_r),"LineWidth",2);
    hold on
end
hold off
grid on
set(gca,'linewidth',1.5,'FontSize',14)
ylim([0 1])
xlim([10 500])
legend('\sigma=0.01nm','\sigma=0.1nm','\sigma=1nm','\sigma=10nm','Location','northeast','FontSize',15)
legend boxoff
xlabel('TA','FontSize',16)

nexttile
for irms = 1:length(rms_plotX)
    rms_plot = rms_plotX(irms);
    ind_r = find(abs(rms-rms_plot)<1e-20);
    semilogx(TX,p_ave(2,:,ind_r),"LineWidth",2)
    hold on
end
hold off
grid on
set(gca,'linewidth',1.5,'FontSize',14)
ylim([0 1])
xlim([10 500])
legend('\sigma=0.1nm','\sigma=1nm','\sigma=10nm','\sigma=100nm','Location','northeast','FontSize',15)
legend boxoff
xlabel('LA','FontSize',16)

ylabel(t,"$ p_{ave} $",'Interpreter','latex','FontSize',22)
xlabel(t,"$ T $, K",'Interpreter','latex','FontSize',22)
%% Conductivity of Nanowire
%% For diameters

dX = (1:0.1:10)'*[1e-8 1e-7 1e-6 1e-5];
dX = unique(dX(:));
kw_mode = zeros(2,length(dX),length(rms));
for mode = 1:2
    omegaM3 = omegaX_slack(mode,:);
    vg = vpX_slack(mode,:);
    vg(1) = vg(2);
    DX = 1/(2*pi^2)*omegaM3.^2./(vg.^3);    
    DX(1) = DX(2); % approximation        
    T_fix = 300;
    ind_T = abs(TX-T_fix)<1e-20;
    XX = h*omegaM3./(k_B*T_fix);
    CphX = Cph(XX);
    if isnan(CphX(1))
        CphX(1) = CphX(2); % approximation
    end
    tauX3 = (TAU_N_Her(omegaM3,T_fix)+TAU_U_slack(omegaM3,T_fix)+TAU_imp(omegaM3)).^-1;
    tauXX3 = tauX3(mode,:);
    tauXX3(1) = 0;
    for id = 1:length(dX)
        KnX = (tauXX3.*v0(mode))./dX(id);
        KnX(1) = KnX(2);   % approximation
        for ir = 1:length(rms)
            FX = interp2(Kn,pX,F,KnX,p_ave(mode,ind_T,ir));
%             FX1 = diag(FX);
            y = CphX.*DX.*vg.^2.*tauXX3.*FX;
            kw_mode(mode,id,ir) = trapz(omegaM3,y);
        end
    end     
end

kX_w_Ts(:,:) = 1/3*(2*kw_mode(1,:,:) + kw_mode(2,:,:));
%%

%%%%%%%%%$$ Plot %%%%%%%%%$$
figure('OuterPosition',[0 0 800 600]);
rms_plotX = [0.05 0.1 5]*10^-9;

for ir = 1:length(rms_plotX)
    rms_plot = rms_plotX(ir);
    ind_r = find(abs(rms-rms_plot)<1e-20);
    ind_r = ind_r(1);
    kX_w_plot = kX_w_Ts(:,ind_r);  
    semilogx(dX'*1e9,kX_w_plot);
    hold on
end

ind_T = find(abs(TX-T_fix)<1e-20);
k_bulk_fix = ones(1,length(dX))*k_bulk4(ind_T);
semilogx(dX'*1e9,k_bulk_fix);
hold on

d_exp = [37 56 115];
k_exp = [18 26 40];
 
semilogx(d_exp,k_exp,'s','Color','r')

set(gca,'linewidth',1.5,'FontSize',16)
legend('\sigma=0.05nm','\sigma=0.1nm','\sigma=5.0nm','bulk','Эксперимент Li','Location','east','FontSize',16)
legend boxoff
grid on
xlabel("$d$, nm",'Interpreter','latex','FontSize',22)
ylabel('$\kappa_{wire}$, W/(mK)','Interpreter','latex','FontSize',22)

axes('position',[0.55,0.2,0.3,0.3])
rms_plotX = [0.5 1 5]*10^-9;
for ir = 1:length(rms_plotX)
    rms_plot = rms_plotX(ir);
    ind_r = find(abs(rms-rms_plot)<1e-20);
    ind_r = ind_r(1);
    kX_w_plot = kX_w_Ts(:,ind_r);  
    semilogx(dX'*1e9,kX_w_plot);
    hold on
end
loglog(d_exp,k_exp,'s','Color','r')
xlim([30 120])
legend('\sigma=0.5nm','\sigma=1nm','\sigma=5nm','Location','best','FontSize',12)
legend boxoff
grid on
%%

dX_new = [37 56 115]*1e-9;
kw_mode = zeros(2,length(TX),length(dX_new),length(rms));
kX_w = zeros(length(TX),length(dX_new),length(rms));

for mode = 1:2
    omegaM3 = omegaX_slack(mode,:);
    vg = vpX_slack(mode,:);
    vg(1) = vg(2);
    DX = 1/(2*pi^2)*omegaM3.^2./(vg.^3);    
    DX(1) = DX(2); % approximation
    for iT = 1:length(TX)        
        T = TX(iT);
        XX = h*omegaM3./(k_B*T);
        CphX = Cph(XX);
        if isnan(CphX(1))
            CphX(1) = CphX(2); % approximation
        end
        tauX3 = (TAU_N_Her(omegaM3,T)+TAU_U_slack(omegaM3,T)+TAU_imp(omegaM3)).^-1;
        tauXX3 = tauX3(mode,:);
        tauXX3(1) = 0;
        for id = 1:length(dX_new)
            KnX = (tauXX3.*v0(mode))./dX_new(id);
            KnX(1) = KnX(2);   % approximation
            for ir = 1:length(rms)
                FX = interp2(Kn,pX,F,KnX,p_ave(mode,iT,ir));
                y = CphX.*DX.*vg.^2.*tauXX3.*FX;
                kw_mode(mode,iT,id,ir) = trapz(omegaM3,y);
            end
        end 
    end    
end

kX_w(:,:,:) = 1/3*(2*kw_mode(1,:,:,:) + kw_mode(2,:,:,:));
%%
%%%%%%%%%$$ Plot %%%%%%%%%$$
figure('OuterPosition',[0 0 800 600]);
rms_plotX = [2 3 6]*10^-9;

for id = 1:length(dX_new)
    rms_plot = rms_plotX(id);
    ind_r = find(abs(rms-rms_plot)<1e-20);
    ind_r = ind_r(1);
    plot(TX',kX_w(:,id,ind_r));
    hold on
end


load exp_wire_li.mat % Experiment Li
T_exp = T_exp_li;
k_exp = k_exp_li;


T_exp1 = T_exp(2,11:end);
T_exp2 = T_exp(3,8:end);
T_exp3 = T_exp(4,12:end);
k_exp1 = k_exp(2,11:end);
k_exp2 = k_exp(3,8:end);
k_exp3 = k_exp(4,12:end);
plot(T_exp1,k_exp1,'s','Color','r','MarkerSize',6)
hold on
plot(T_exp2,k_exp2,'s','Color','r','MarkerSize',6)
hold on
plot(T_exp3,k_exp3,'s','Color','r','MarkerSize',6)

hold off
set(gca,'linewidth',2,'FontSize',20)
legend('d=37nm,\sigma=2nm','d=56nm,\sigma=3nm','d=115nm,\sigma=6nm','Эксперимент Li','Location','northeast','FontSize',14)
legend boxoff

ylim([0 70])
xlim([100 350])
grid on
xlabel("$T$, K",'Interpreter','latex','FontSize',20)
ylabel('$\kappa_{wire}$, W/(mK)','Interpreter','latex','FontSize',20)

%%
dX_exp_lim = [54 57 68 70 78 84 93 99 50 71 73 78 80 54 66 69 70 78 84 92 100 61 64 69]*10^-9;
rms_exp_lim = [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 2.75 2.75 2.75 2.75 2.75 3.25 3.25 3.25 3.25 3.25 3.25 3.25 3.25...
    4 4 4]*10^-9;
error_lim = [1 1 1 1 1 1 1 1 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.5 0.5 0.5];
error_li = [0.5 0.6 1.1];
figure('OuterPosition',[0 0 800 600]);
errorbar(dX_new*1e9,rms_plotX*1e9,error_li,'rs','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on
errorbar(dX_exp_lim*1e9,rms_exp_lim*1e9,error_lim,'ro','MarkerSize',8)
hold on
fplot(@(d) 0.04*d,':','linewidth',2)
hold on
fplot(@(d) 0.07*d,'--','linewidth',2)
hold off
set(gca,'linewidth',2,'FontSize',18)
grid on
xlabel("$d$, nm",'Interpreter','latex','FontSize',22)
ylabel("$\sigma$, nm",'Interpreter','latex','FontSize',22)
legend('This work','Experiment Lim','\sigma=0.04d','\sigma=0.07d','Location','northeast','FontSize',14)
legend boxoff
xlim([0 150])
ylim([0 8])
%%
ind_r = find(abs(rms-rms_plotX(1))<1e-20);
ind_r = ind_r(1);
k_wire_1 = interp1(TX,kX_w(:,1,ind_r),T_exp1);
error1 = abs(k_exp1-k_wire_1)./k_exp1*100;
ind_r = find(abs(rms-rms_plotX(2))<1e-20);
ind_r = ind_r(1);
k_wire_2 = interp1(TX,kX_w(:,2,ind_r),T_exp2);
error2 = abs(k_exp2-k_wire_2)./k_exp2*100;
ind_r = find(abs(rms-rms_plotX(3))<1e-20);
ind_r = ind_r(1);
k_wire_3 = interp1(TX,kX_w(:,3,ind_r),T_exp3);
error3 = abs(k_exp3-k_wire_3)./k_exp3*100;

figure('OuterPosition',[0 0 800 600])
plot(T_exp1,error1,':x')
hold on
plot(T_exp2,error2,'-.o')
hold on
plot(T_exp3,error3,'--s')
hold off
xlim([100 350])
set(gca,'linewidth',1.5,'FontSize',16)
grid on
legend('d=37nm,\sigma=2nm','d=56nm,\sigma=3nm','d=115nm,\sigma=6nm','Location','northeast','FontSize',14)
legend boxoff
xlabel("$T$, K",'Interpreter','latex','FontSize',22)
ylabel('Relative error','Interpreter','latex','FontSize',22)