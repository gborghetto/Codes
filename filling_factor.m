clear 
clc
clf

c=0.5;
G = 6.674*10^(-8);
Ra = 2*6.955*10^10;
M = 1.989*10^33;
l = 1;
l2 = 3;   
M_acc = 1e-8 * (1.988*1e33/(60*60*24*365)) ;  

F1 = 3000/1000;
ds = 0.01;

Rt1 = [4.1;5.1];
Rs = 6;
ra1 = Rt1.*Ra/c;
th = pi/2;
Phi = 0:2*pi/100:2*pi;

phase_d = 0.0; 
phase_o = 0.0;
psi_d = (1-phase_d)*2*pi;
psi_o = (1-phase_o)*2*pi;
beta_d = 0*pi/180;
tilt_o = (0:1:90)*pi/180;
Bdip = 1000;
Boct = 3000;

% t = 0:pi/62:pi;
% p = 0:2*pi/62:2*pi;
% 
% tstep = fliplr(logspace(log10(0.001),log10(0.01),10));
% pstep = fliplr(logspace(log10(0.003),log10(0.006),10));

% ar = 2*pi.*(cos(t(1:(numel(t)-1)))-cos(t(2:end)));
% ar2 = ar(2:numel(ar))-ar(1:numel(ar)-1);
% arr = ar2(1:numel(ar2)-1)-ar2(2:numel(ar2));
% pstep = 0.1;
% tstep = 0.04;
% for i = 1:30
%     tstep(i+1) = tstep(i)-arr(i);
%     pstep(i+1) = pstep(i)-arr(i);
% end

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

% tstep = 0.006:-0.005/90:0.001;
% pstep = 0.06:-0.057/90:0.003;
tstep = 0.006:-0.0045/90:0.0015;
pstep = 0.06:-0.059/90:0.001;


figure(1)
hold on
set(gca,'FontSize',26)

f_acc = zeros(1,numel(tilt_o));
for m = 1:numel(tilt_o)
    t1 = [];
    t2 = [];
    p1 = [];
    p2 = [];
    beta_o = tilt_o(m);
    for p = 1:numel(Phi)
        for j = 1:numel(Rt1)

                theta = th;
                phi = Phi(p);
                R = Rt1(j);
                theta_loop = [theta];
                phi_loop = [phi];
                rloop =[R];
                md = Bdip/(l+1);
                mo = Boct/(l2+1);
                iterations = 0;
                while R>1  && R<=Rt1(j)
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
                        % for pure dipole comment out this part
                        %
                        for k = 0:N2
                            Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
                            Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
                            Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
                        end
                        %
                        if R<=Rs
                            Btdip = Bdip*(1/R)^(l+2)*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                            Btoct = Boct*(1/R)^(l2+2)*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                            Bfdip = Bdip*(1/R)^(l+2)*Pf_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                            Bfoct = Boct*(1/R)^(l2+2)*Pf_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                        elseif R>Rs
                            Btdip =0;
                            Btoct=0;
                            Bfdip = 0;
                            Bfoct=0;
                        end
%                         Btdip = Bdip*(1/R)^(l+2)*Pt_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%                         Btoct = Boct*(1/R)^(l2+2)*Pt_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                        Brdip =  Bdip*(1/R)^(l+2)*Pr_d*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                        Broct =  Boct*(1/R)^(l2+2)*Pr_o*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%                         Bfdip = Bdip*(1/R)^(l+2)*Pf_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%                         Bfoct = Boct*(1/R)^(l2+2)*Pf_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
        
                        Bt = Btdip + Btoct;
                        Br = Brdip + Broct;
                        Bf = Bfdip + Bfoct;
    
                        B = (Bt^2+Br^2+Bf^2)^0.5;
    
                        dr = ds*Br/B;
                        dtheta = ds*Bt/(R*B);
                        dphi = ds*Bf/(B*R*sin(theta));
                        R = R-dr;
                        theta = theta - dtheta;
                        phi = phi - dphi;
    
                        rloop = [rloop;R];
                        theta_loop = [theta_loop;theta];
                        phi_loop = [phi_loop;phi]; 
    
                        iterations = iterations+1;
                end
                if iterations > 100
                    t1 = [t1;theta];
                    p1 = [p1;phi];
                end
                Rx = rloop.*sin(theta_loop).*cos(phi_loop);
                Rz = rloop.*cos(theta_loop);
                Ry = rloop.*sin(theta_loop).*sin(phi_loop);
                
        
        end
    end
    
    for p = 1:numel(Phi)
        for j = 1:numel(Rt1)

                theta = th;
                phi = Phi(p);
                R = Rt1(j);
                theta_loop = [theta];
                rloop =[R];
                phi_loop = [phi];
                md = Bdip/(l+1);
                mo = Boct/(l2+1);
                iterations = 0;
                while R>1 && R<=Rt1(j)
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
                        % for pure dipole comment out this part
                        %
                        for k = 0:N2
                            Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
                            Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
                            Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
                        end
                        %
        
                        if R<=Rs
                            Btdip = Bdip*(1/R)^(l+2)*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                            Btoct = Boct*(1/R)^(l2+2)*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                            Bfdip = Bdip*(1/R)^(l+2)*Pf_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                            Bfoct = Boct*(1/R)^(l2+2)*Pf_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                        elseif R>Rs
                            Btdip =0;
                            Btoct=0;
                            Bfdip = 0;
                            Bfoct=0;
                        end
%                         Btdip = Bdip*(1/R)^(l+2)*Pt_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%                         Btoct = Boct*(1/R)^(l2+2)*Pt_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                        Brdip =  Bdip*(1/R)^(l+2)*Pr_d*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                        Broct =  Boct*(1/R)^(l2+2)*Pr_o*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%                         Bfdip = Bdip*(1/R)^(l+2)*Pf_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%                         Bfoct = Boct*(1/R)^(l2+2)*Pf_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
        
                        Bt = Btdip + Btoct;
                        Br = Brdip + Broct;
                        Bf = Bfdip + Bfoct;
    
                        B = (Bt^2+Br^2+Bf^2)^0.5;
    
                        dr = ds*Br/B;
                        dtheta = ds*Bt/(R*B);
                        dphi = ds*Bf/(B*R*sin(theta));
                        R = R+dr;
                        theta = theta + dtheta;
                        phi = phi + dphi;
    
                        rloop = [rloop;R];
                        theta_loop = [theta_loop;theta];
                        phi_loop = [phi_loop;phi];
    
                        iterations = iterations+1;
                end
                if iterations > 100
                    t2 = [t2;theta];
                    p2 = [p2;phi];
                end
                Rx = rloop.*sin(theta_loop).*cos(phi_loop);
                Rz = rloop.*cos(theta_loop);
                Ry = rloop.*sin(theta_loop).*sin(phi_loop);
                     
        end
    end
    
    footpoints_t = [t1;t2];
    footpoints_p = [p1;p2]; 

    for i = 1:numel(footpoints_t)
        if footpoints_t(i)<0
            footpoints_t(i)=footpoints_t(i)+(2*pi);
        elseif footpoints_t(i)>(pi)
            footpoints_t(i)=footpoints_t(i)-(2*pi);
        end
        if footpoints_p(i)<0
            footpoints_p(i)=footpoints_p(i)+(2*pi);
        elseif footpoints_p(i)>(2*pi)
            footpoints_p(i)=footpoints_p(i)-(2*pi);
        end
    end
    
%     tt1 = abs(mean(footpoints_t(1:length(t1)-1)-footpoints_t(2:length(t1))));
%     pp1 = abs(mean(footpoints_p(1:length(p1)-1)-footpoints_p(2:length(p1))));
%     tstep = mean(abs(footpoints_t(2:end)-footpoints_t(1:(end-1))));
%     pstep =  mean(abs(footpoints_p(3:end)-footpoints_p(1:(end-2))));
    ttt = pi/tstep(m);
    ppp = 2*pi/pstep(m);
    tt = 0:pi/(floor(ttt)):pi;
    pp = 0:2*pi/(floor(ppp)):2*pi;

    area = zeros(numel(footpoints_t),1);
    
    
    for i = 1:numel(footpoints_t)
        a = 0;
        a1 = 0;
        for t = 1:(numel(tt)-1)
            if footpoints_t(i)>=tt(t) && footpoints_t(i)<tt(t+1)
                a = cos(tt(t))-cos(tt(t+1));
            elseif footpoints_t(i)==tt(end)
                a = cos(tt(end))-cos(tt(end-1));
            end
        end
        for t = 1:(numel(pp)-1)
            if footpoints_p(i)>=pp(t) && footpoints_p(i)<pp(t+1)
                a1 = pp(t+1)-pp(t);
            elseif footpoints_p(i)==pp(end)
                a1 = pp(end)-pp(end-1);
            end
        end
        area(i) = a*a1;
    end
    
    
    
    f = sum(area);
    f_acc(1,m) = 100*f/(4*pi);
end


figure(1)
hold on
set(gca,'FontSize',26)
set(gca,'TicklabelInterpreter','latex')

scatter(tilt_o*180/pi,f_acc,50, 'filled')
ylim([0,0.5])
xlim([0,90])
xlabel('$\beta_{oct}$ (degrees)','Interpreter','latex')
ylabel('Filling factor ($\%$)','Interpreter','latex')