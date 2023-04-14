clear all
clc
clf
ax = gca;
set(gca,'TicklabelInterpreter','latex')
set(gca,'FontSize',28)
c=1/2;
G = 6.674e-8;
Rsun =6.955e10;
R = 2*Rsun;     
Msun = 1.989e33;
M = 1*Msun;
M_acc = 1*10^(-9) * (Msun/(60*60*24*365)) ;   
T = 3.56*24*60*60;   
v = 2*pi*R/T;
w = v/R;
R_co = ((G*M/(w^2))^(1/3))/R;

ds = 0.01;

l = 1;
l2 = 3;
Rs = 1000;
phase_d = 0.0;
phase_o = 0.50;
psi_d = (1-phase_d)*2*pi;
psi_o = (1-phase_o)*2*pi;
beta_d = 0*pi/180;
beta_o = 30*pi/180;

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
t2 = [];
p1 = [];
p2 = [];
Bdip = 3000;
Boct = 1000;
mo = Boct/(l2+1);
md = Bdip/(l+1);   
theta = pi/2;
p = 0:2*pi/30:2*pi;
r = (1:0.001:30);
rr = r*R;
n = numel(r);
Rt = zeros(numel(p),1);
BB = zeros(1,numel(p));
for i = 1:numel(p)
            phi = p(i);
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
    
            n = numel(r);

            Btdip = Bdip.*(1./r).^(l+2).*Pt_d.*((-r.^(2*l+1)+Rs^(2*l+1))./(l+(l+1)*Rs^(2*l+1)));
            Btoct = Boct.*(1./r).^(l2+2).*Pt_o.*((-r.^(2*l2+1)+Rs^(2*l2+1))./(l2+(l2+1)*Rs^(2*l2+1)));
            Bfdip = Bdip.*(1./r).^(l+2).*Pf_d.*((-r.^(2*l+1)+Rs^(2*l+1))./(l+(l+1)*Rs^(2*l+1)));
            Bfoct = Boct.*(1./r).^(l2+2).*Pf_o.*((-r.^(2*l2+1)+Rs^(2*l2+1))./(l2+(l2+1)*Rs^(2*l2+1)));


            Brdip =  Bdip.*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
            Broct =  Boct.*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
            
%             Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
%             Btoct = Boct(j)/(l2+1).*(1./r).^(l2+2)*Pt_o;
%             Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
%             Broct =  Boct(j).*(1./r).^(l2+2)*Pr_o;
%             Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
%             Bfoct = Boct(j)/(l+1).*(1./r).^(l2+2)*Pf_o;
            
    
            Bt = Btdip + Btoct;
            Br = Brdip + Broct;
            Bf = Bfdip + Bfoct;
        
            B = (Bt.^2+Br.^2+Bf.^2).^0.5;
            
            p_mag = 1/(8*pi).*B.^2;
            p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
            P = p_mag-p_ram;
            while P(n)<0
                n = n-1;
            end
            Rt(i) = c*r(n);
            BB(i) = B(n)^2;
end
b = mean(BB);
p_mag = 1/(8*pi)*b;
a = p_mag-p_ram;
idx = find(abs(a)==min(abs(a)));
rt = c*r(idx);
    

figure(1)
hold on
set(gca,'TicklabelInterpreter','latex')
set(gca,'FontSize',28)
axis equal
Rs = 10;
Phi = p;
th = pi/2;
RRt = zeros(numel(Rt),1);
RRt(:) = R_co-0.01;
Rt1 = [Rt RRt];
for s =1:2
    RR = Rt1(:,s);
    for i = 1:numel(RR)
        
            Rmax = RR;
        
            theta = th;
            R = Rmax(i);
            rloop =[R];
            phi = Phi(i);
            theta_loop = [theta];
            phi_loop = [phi];
            Bloop = [];
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
                        Bloop=[Bloop;B];
            end
            Rx = rloop.*sin(theta_loop).*cos(phi_loop);
            Rz = rloop.*cos(theta_loop);
            Ry = rloop.*sin(theta_loop).*sin(phi_loop);
            
            if round(R,1)==1.0
                t1 = [t1;theta];
                p1 = [p1;phi];
            end
            plot3(ax,Rx,Ry,Rz,'black','linewidth',2)
    end       
end

for s =1:2
    RR = Rt1(:,s);
    for i = 1:numel(RR)
        
            Rmax = RR;
        
            theta = th;
            R = Rmax(i);
            rloop =[R];
            phi = Phi(i);
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
            plot3(ax,Rx,Ry,Rz,'black','linewidth',2)        

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



theta = 0:0.01:pi;
p = (pi:0.01:3*pi);
B = zeros(numel(theta),numel(p));
% Bdip = 1000;
% Boct = 3000;
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
    Br = Brdip + Broct;

    B(:,i) = Br; 

end
theta = 0:pi/100:pi;
p = (pi:2*pi/100:3*pi);
BB = zeros(numel(theta),numel(p));
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
    BB(:,i) = Brdip + Broct;

end
% [x,y,z] = sphere;
% s1 = surf(x,y,z,"FaceColor",[1 1 0],'EdgeColor','none');

%  
% footpoints_p2 = footpoints_p+pi;
% for f = 1:numel(footpoints_p2)
%     if footpoints_p2(f)>2*pi
%         footpoints_p2(f) = footpoints_p2(f)-2*pi;
%     end
% end
[x,y,z] = sphere(100);
s1 = surf(ax,x,y,z,flipud(BB)/1000,'EdgeColor','none');
colormap("turbo")
colorbar;
% xlabel('Rx[ $R_{\star}$]','Interpreter','latex')
% zlabel('Rz [$R_{\star}$]','Interpreter','latex')
% ylabel('Rz [$R_{\star}$]','Interpreter','latex')
a = colorbar;
a.TickLabelInterpreter = 'latex';
a.Label.String = 'kG';
set(a.XLabel,{'String','Rotation','Position'},{'kG',0,[0.9 3.6]},'Interpreter','latex')
xlabel('x [$R_{\star}$]','Interpreter','latex')
zlabel('z [$R_{\star}$]','Interpreter','latex')
ylabel('y [$R_{\star}$]','Interpreter','latex')
% xlim([-3,3])
% zlim([-2,2])
hold off
figure(2)
hold on
set(gca,'TicklabelInterpreter','latex')
set(gca,'FontSize',28)
tt = [-90 90];
pp = [0 360];
image(pp,tt,flipud(B/1000),'CDataMapping','scaled')
a = colorbar;
a.Label.String = 'kG';
set(a.XLabel,{'String','Rotation','Position'},{'kG',0,[0 0]})
colormap("turbo")
scatter(footpoints_p*180/pi,(footpoints_t*180/pi)-90,120,"w","filled")
xlim([0;360])
ylim([-90;90])
xticks([0,90,180,270,360])
yticks([-90,-45,0,45,90])

xlabel("Longitude (degrees)",'Interpreter','latex')
ylabel("Latitude (degrees)",'Interpreter','latex')
hold off