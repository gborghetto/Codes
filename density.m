clear all
clc
clf

c=0.5;
G = 6.674*10^(-8);
Ra = 1.1*6.955*10^10;
M = 0.8*1.989*10^33;    
M_acc = 0.5*10^(-9)*(1.988*1e33/(60*60*24*365)) ;  

l = 1;
l2 = 3;
Rs = 8.3;

phase_d = 0.75;
phase_o = 0.15;
psi_d = (1-phase_d)*2*pi;
psi_o = (1-phase_o)*2*pi;
beta_d = 5*pi/180;
beta_o = 175*pi/180;

th = pi/2;
Phi = 0:0.1:2*pi;

if rem(l,2) == 1
    N1 = (l-1)/2;
else
    N1 = l/2;
end 
if rem(l2,2)== 1
    N2 = (l2-1)/2;
else
    N2 = l2/2;
end 

Bdip = 1200;
Boct = 1600;
F = Boct./Bdip;

ds = 0.01;

rt = [2.5;8.25];


f_acc = 1.5;


figure(1)
hold on
set(gca,'TicklabelInterpreter','latex')
set(gca,'FontSize',28)
for i = 1:numel(rt)
        mo = Boct/(l2+1);
        md = Bdip/(l+1); 
        theta = pi/2;
        phi = pi;
        R = rt(i);
        vv = ((2*G*M)/(R*Ra))^0.5*(1/R-1/rt(i))^0.5;
        rloop =[R];
        theta_loop = [theta];
        v = [vv];
        Bloop = [];
        tt = 0;
        while R>1 
                mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
                mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
                mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
                mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
                mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
                mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
            
                k = 0;
                Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
                Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
                Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
            
                Pr_o = 0;
                Pt_o = 0;
                Pf_o = 0;
                for k = 0:N2
                    Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
                    Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
                    Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
                end
        
              
                Btdip = Bdip*(1/R)^(l+2)*Pt_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                Btoct = Boct*(1/R)^(l2+2)*Pt_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                Brdip =  Bdip*(1/R)^(l+2)*Pr_d*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                Broct =  Boct*(1/R)^(l2+2)*Pr_o*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                Bfdip = Bdip*(1/R)^(l+2)*Pf_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                Bfoct = Boct*(1/R)^(l2+2)*Pf_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));

                Bt = Btdip + Btoct;
                Br = Brdip + Broct;
                Bf = Bfdip + Bfoct;

                B = (Bt^2+Br^2+Bf^2)^0.5;
                dr = ds*Br/B;
                dtheta = ds*Bt/(R*B);
                R = R-dr;
                theta = theta - dtheta;
                rloop = [rloop;R];
                theta_loop = [theta_loop;theta];
                vv = ((2*G*M)/(R*Ra))^0.5*(1/R-1/rt(i))^0.5;
                v = [v;vv];
                Bloop = [Bloop;B];
        end
        A = f_acc.*(4*pi*Ra^2);
        rho1 = (Bloop./B).*(M_acc./(v(1:(end-1)).*A));
        rho1 = rho1/(1e-24);
        plot(rloop(1:(end-1)),(log10(rho1)),"LineWidth",2)
end

ylabel('$log (\rho_N) [cm^{-3}]$','Interpreter','latex')
xlabel('$r [R_{\star}]$','Interpreter','latex')
xlim([1;8.5])
hold off


