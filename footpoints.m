clear all
clc
clf
set(gca,'TicklabelInterpreter','latex')
set(gca,'FontSize',26)
l = 1;
l2 = 3;
c=1/2;
G = 6.674e-8;
Rsun =6.955e10;
Ra = 1.1*Rsun;     
Msun = 1.989e33;
M = 0.8*Msun;
M_acc = 0.5e-9 * (Msun/(60*60*24*365)) ; 
T = 3.56*24*60*60;       
v = 2*pi*Ra/T;
w = v/Ra;
R_co = ((G*M/(w^2))^(1/3))/Ra;
Bdip = 500;      
Boct = 2800;

Rs = R_co;
Rt1 = [7.2;Rs];  
% ;2.9:3.1;3.3;3.5;3.7;3.9;4.1
th = pi/2;
% th = 0:pi/6:pi;
Phi = 0:2*pi/100:2*pi;
ds = 0.01;

phase_d = 0.75;
phase_o = 0.15;
psi_d = (1-phase_d)*2*pi;
psi_o = (1-phase_o)*2*pi;
beta_d = 5*pi/180;
beta_o = 175*pi/180;
md = Bdip/(l+1);
mo = Boct/(l2+1);
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
t1 = [];
t2 =[];
p1 = [];
p2 = [];

% figure(1)
% hold on
% set(gca,'FontSize',26)
% axis equal
for p = 1:numel(Phi)
    for t=1:numel(th)
        Rmax = Rt1;
        for i = 1:numel(Rmax)
            theta = th(t);
            R = Rmax(i);
            rloop =[R];
            phi = Phi(p);
            theta_loop = [theta];
            phi_loop = [phi];
            while R>1 && R<=Rmax(i)
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
                        
    %                     Btdip = Bdip*(1/R)^(l+2)*Pt_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
    %                     Btoct = Boct*(1/R)^(l2+2)*Pt_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                        Brdip =  Bdip*(1/R)^(l+2)*Pr_d*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                        Broct =  Boct*(1/R)^(l2+2)*Pr_o*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
    %                     Bfdip = Bdip*(1/R)^(l+2)*Pf_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
    %                     Bfoct = Boct*(1/R)^(l2+2)*Pf_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
        
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
            end
            Rx = rloop.*sin(theta_loop).*cos(phi_loop);
            Rz = rloop.*cos(theta_loop);
            Ry = rloop.*sin(theta_loop).*sin(phi_loop);

            if round(R,1)==1.0
                t1 = [t1;theta];
                p1 = [p1;phi];
            end
%             plot3(Rx,Ry,Rz,'black','linewidth',2)


        end 
    
    end
end

for p = 1:numel(Phi)
    for t=1:numel(th)
        Rmax = Rt1;
        for i = 1:numel(Rmax)
            theta = th(t);
            R = Rmax(i);
            rloop =[R];
            phi = Phi(p);
            theta_loop = [theta];
            phi_loop = [phi];
            while R>1 && R<=Rmax(i)
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
                        
    %                     Btdip = Bdip*(1/R)^(l+2)*Pt_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
    %                     Btoct = Boct*(1/R)^(l2+2)*Pt_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
                        Brdip =  Bdip*(1/R)^(l+2)*Pr_d*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                        Broct =  Boct*(1/R)^(l2+2)*Pr_o*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
    %                     Bfdip = Bdip*(1/R)^(l+2)*Pf_d*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
    %                     Bfoct = Boct*(1/R)^(l2+2)*Pf_o*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
        
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
            end
            Rx = rloop.*sin(theta_loop).*cos(phi_loop);
            Rz = rloop.*cos(theta_loop);
            Ry = rloop.*sin(theta_loop).*sin(phi_loop);
            if round(R,1)==1.0
                t2 = [t2;theta];
                p2 = [p2;phi];
            end
%             plot3(Rx,Ry,Rz,'black','linewidth',2)        

        end
            
    end
end

footpoints_t = [t1;t2];
footpoints_p = [p1;p2];

theta = 0:pi/100:pi;
p = 0:2*pi/100:2*pi;
B = zeros(numel(theta),numel(p));
R = 1;

md = Bdip/(l+1);
mo = Boct/(l2+1);


for i = 1 : numel(p)
    phi = p(i);
    mr_d = md.*sin(theta).*cos(phi).*cos(psi_d).*sin(beta_d)+md.*sin(theta).*sin(phi).*sin(psi_d).*sin(beta_d)+md.*cos(theta).*cos(beta_d);
    mt_d = md.*cos(theta).*cos(phi).*cos(psi_d).*sin(beta_d)+md.*cos(theta).*sin(phi).*sin(psi_d).*sin(beta_d)-md.*sin(theta).*cos(beta_d);
    mf_d = -md.*sin(phi).*cos(psi_d).*sin(beta_d)+md.*cos(phi).*sin(psi_d).*sin(beta_d);
    mr_o = mo.*sin(theta).*cos(phi).*cos(psi_o).*sin(beta_o)+mo.*sin(theta).*sin(phi).*sin(psi_o).*sin(beta_o)+mo.*cos(theta).*cos(beta_o);
    mt_o = mo.*cos(theta).*cos(phi).*cos(psi_o).*sin(beta_o)+mo.*cos(theta).*sin(phi).*sin(psi_o).*sin(beta_o)-mo.*sin(theta).*cos(beta_o);
    mf_o = -mo.*sin(phi).*cos(psi_o).*sin(beta_o)+mo.*cos(phi).*sin(psi_o).*sin(beta_o);

    k = 0;
    Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k)).*(mr_d./md).^(l-2*k);
    Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1)).*(mr_d./md).^(l-2*k-1).*mt_d./md;
    Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1)).*(mr_d./md).^(l-2*k-1).*mf_d./md;

    Pr_o = 0;
    Pt_o = 0;
    Pf_o = 0;
    for k = 0:N2
      Pr_o = Pr_o+((-1)^k.*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k)).*(mr_o./mo).^(l2-2*k);
      Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1)).*(mr_o./mo).^(l2-2*k-1).*mt_o./mo;
      Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1)).*(mr_o./mo).^(l2-2*k-1).*mf_o./mo;
    end

    Brdip =  Bdip.*(1/R)^(l+2).*Pr_d.*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
    Broct =  Boct.*(1/R)^(l2+2).*Pr_o.*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
    Bt = 0;
    Bf = 0;
    B(:,i) = Brdip + Broct;

end
% 
% [x,y,z] = sphere(100);
% s1 = surf(x,y,z,flipud(B)/1000,'EdgeColor','none');
% colormap("turbo")
% colorbar;
% % xlabel('Rx[ $R_{\star}$]','Interpreter','latex')
% % zlabel('Rz [$R_{\star}$]','Interpreter','latex')
% % ylabel('Rz [$R_{\star}$]','Interpreter','latex')
% a = colorbar;
% a.TickLabelInterpreter = 'latex';
% a.Label.String = 'kG';
% set(a.XLabel,{'String','Rotation','Position'},{'kG',0,[0.9 3.6]},'Interpreter','latex')
% xlabel('x [$R_{\star}$]','Interpreter','latex')
% zlabel('z [$R_{\star}$]','Interpreter','latex')
% ylabel('y [$R_{\star}$]','Interpreter','latex')
% 
% % xlim([-3,3])
% % zlim([-2,2])
% 
% hold off


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


step1 = pi/0.004;
step2 = 2*pi/0.06;

tt = 0:pi/(floor(step1)):pi;
pp = 0:2*pi/(floor(step2)):2*pi;


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

scatter(footpoints_p,footpoints_t)

f = sum(area);
f_acc = 100*f/(4*pi);
% end
numel(footpoints_p(3:end))
numel(footpoints_p(1:(end-2)))
