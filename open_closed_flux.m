clear 
clc
clf

c=0.5;
G = 6.674*10^(-8);
Ra = 2*6.955*10^10;
M = 1.989*10^33;
l = 1;
l2 = 3;   

ds = 0.01;

Rsurf = 1:0.1:8;
Phi = 0:2*pi/250:2*pi;
Theta = 0:pi/200:pi;


Bdip = 1000;      
Boct = 3000;

phase_d = 0.65;
phase_o = 0.5;
psi_d = (1-phase_d)*2*pi;
psi_o = (1-phase_o)*2*pi;
beta_d = 25*pi/180;
beta_o = 50*pi/180;
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

figure(1)
hold on
set(gca,'FontSize',26)
open_flux = zeros(1,numel(Rsurf));
closed_flux = zeros(1,numel(Rsurf));
for s = 1:numel(Rsurf)
    Rs = Rsurf(s);
    flux_open = 0;
    flux_closed = 0;
    for p = 1:numel(Phi)-1
        phi1= Phi(p);
        phi2 = Phi(p+1);
        for t=1 :numel(Theta)-1
            R = 1;
            theta1 = Theta(t);
            theta2 = Theta(t+1);
            phi_cell = phi1:(phi2-phi1)/10:phi2;
            theta_cell = theta1:(theta2-theta1)/10:theta2;
            for pp = 1:numel(phi_cell)
                phi = phi_cell(pp);
                theta = theta_cell;
                mr_d = md.*sin(theta).*cos(phi).*cos(psi_d).*sin(beta_d)+md.*sin(theta).*sin(phi).*sin(psi_d).*sin(beta_d)+md.*cos(theta).*cos(beta_d);
                mr_o = mo.*sin(theta).*cos(phi).*cos(psi_o).*sin(beta_o)+mo.*sin(theta).*sin(phi).*sin(psi_o).*sin(beta_o)+mo.*cos(theta).*cos(beta_o);
                
                k = 0;
                Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k)).*(mr_d./md).^(l-2*k);
              
                Pr_o = 0;
             
                for k = 0:N2
                  Pr_o = Pr_o+((-1)^k.*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k)).*(mr_o./mo).^(l2-2*k);
                end

                Brdip =  Bdip.*(1/R)^(l+2).*Pr_d.*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                Broct =  Boct.*(1/R)^(l2+2).*Pr_o.*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
     
                Br = Brdip + Broct;
   
            end
            Br_av = mean(abs(Br)); 
            Br_cell = Ra^2*Br_av*(phi2-phi1)*(cos(theta1)-cos(theta2));
            R = 1;
            theta = mean([theta1,theta2]);
            phi = mean([phi1,phi2]);

            while R>=1 && R<=Rs
                        mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
                        mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
                        mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
                        mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
                        mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
                        mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
        
    %                     Pr_d = 0;
    %                     Pt_d = 0;
    %                     Pf_d = 0;
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
                        R = R - dr;
                        theta = theta - dtheta;
                        phi = phi - dphi;
        
            end
            if round(R,1)==1.0 && round(theta,6)~=round(mean([theta1,theta2]),6)
                flux_closed = flux_closed + Br_cell;
            elseif round(R,1)==round(Rs,1)
                flux_open = flux_open + Br_cell;
            elseif round(R,1)==1.0 && round(theta,6)==round(mean([theta1,theta2]),6)
                R = 1;
                theta = mean([theta1,theta2]);
                phi = mean([phi1,phi2]);
    
                while R>=1 && R <= Rs
                            mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
                            mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
                            mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
                            mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
                            mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
                            mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
        
    %                         Pr_d = 0;
    %                         Pt_d = 0;
    %                         Pf_d = 0;
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
                            
                            Brdip =  Bdip*(1/R)^(l+2)*Pr_d*((l*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
                            Broct =  Boct*(1/R)^(l2+2)*Pr_o*((l2*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
        %                  
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
                end
                if round(R,1)==1.0 && round(theta,6)~=round(mean([theta1,theta2]),6)
                    flux_closed = flux_closed + Br_cell;
                elseif round(R,1)==round(Rs,1)
                    flux_open = flux_open + Br_cell;
                end
   
            
            end

        end
    end
    open_flux(s) = flux_open;
    closed_flux(s) = flux_closed;
end
figure(1)
hold on
set(gca,'FontSize',26)
plot(Rsurf,open_flux,'LineWidth',3,'Color','blue')
plot(Rsurf,closed_flux,'LineWidth',3,'Color','red')

xlabel("$R_s/R_{\star}$",'Interpreter','latex')
ylabel("Flux [$G$ $cm^2$]",'Interpreter','latex')
xlim([1,8])
% legend("Open flux", "Closed flux",'Interpreter','latex')

