clear all
clc
clf
set(gca,'TicklabelInterpreter','latex')
set(gca,'FontSize',28)
c=1/2;
G = 6.674e-8;
Rsun =6.955e10;
Ra = 1.95*Rsun;     
Msun = 1.989e33;
M = 0.65*Msun;
M_acc = 3*10^(-8) * (Msun/(60*60*24*365)) ;  

l = 1;
l2 = 3;
Rs = 1000;
phi = 0;
psi = 0;
phase_d = 0.65;
phase_o = 0.15;
psi_d = (1-phase_d)*2*pi;
psi_o = (1-phase_o)*2*pi;
beta_d = 20*pi/180;
beta_o = 10*pi/180;
theta = pi/2;
p = 0:0.1:2*pi;

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

% 
% Bdip = 1000;
% Boct = [1000;2000;3000;4000;5000];
% F = Boct./Bdip;
% % figure(1)
% % hold on
% % set(gca,'TicklabelInterpreter','latex')
% % set(gca,'FontSize',28)
% % for j = 1:numel(F)
% %     mo = Boct(j)/(l2+1);
% %     md = Bdip/(l+1);   
% %     r = (1:0.001:30);
% %     rr = r*R;
% %     n = numel(r);
% %     Rt = zeros(1,numel(p));
% %     for i = 1:numel(p)
% %         phi = p(i);
% %         mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
% %         mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
% %         mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
% %         mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
% %         mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
% %         mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
% %     
% %         k = 0;
% %         Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
% %         Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
% %         Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
% %     
% %         Pr_o = 0;
% %         Pt_o = 0;
% %         Pf_o = 0;
% %         for k = 0:N2
% %             Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
% %             Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
% %             Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
% %         end
% % 
% %         n = numel(r);
% % %         Btdip = [];
% % %         Btoct = [];
% % %         for k = 1:n
% % %             if r(k)<=Rs
% % %                 R = r(k);
% % %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% % %             elseif r(k)>Rs
% % %                 Btdip(k) =0;
% % %                 Btoct(k)=0;
% % %             end
% % %         end
% % % 
% % %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% %         
% %         Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
% %         Btoct = Boct(j)/(l2+1).*(1./r).^(l2+2)*Pt_o;
% %         Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
% %         Broct =  Boct(j).*(1./r).^(l2+2)*Pr_o;
% %         Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
% %         Bfoct = Boct(j)/(l+1).*(1./r).^(l2+2)*Pf_o;
% %         
% % 
% %         Bt = Btdip + Btoct;
% %         Br = Brdip + Broct;
% %         Bf = Bfdip + Bfoct;
% %     
% %         B = (Bt.^2+Br.^2+Bf.^2).^0.5;
% %         
% %         p_mag = 1/(8*pi).*B.^2;
% %         p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
% %         P = p_mag-p_ram;
% %         while P(n)<0
% %             n = n-1;
% %         end
% %         Rt(i) = c*r(n);
% %     end
% %     plot(p*180/pi,Rt,'LineWidth',3)
% % end
% % 
% % xlabel('$\phi$','Interpreter','latex')
% % ylabel('$r_t [R_{\star}]$','Interpreter','latex')
% % legend('$B_{oct}/B_{dip}=1$','$B_{oct}/B_{dip}=2$','$B_{oct}/B_{dip}=3$','$B_{oct}/B_{dip}=4$','$B_{oct}/B_{dip}=5$','Interpreter','latex','Location','NorthEastOutside')
% % title('$R_t$ vs. $\phi$ with $\beta_{dip}=5$ and $\beta_{oct}=45$','Interpreter','latex')
% 
% % Bdip = 1000;
% % Boct = [1000;2000;3000;4000;5000];
% % F = Boct./Bdip;
% % figure(1)
% % hold on
% % set(gca,'TicklabelInterpreter','latex')
% % set(gca,'FontSize',28)
% % for j = 1:numel(F)
% %     mo = Boct(j)/(l2+1);
% %     md = Bdip/(l+1);   
% %     r = (1:0.001:30);
% %     rr = r*R;
% %     n = numel(r);
% %     Rt = zeros(1,numel(p));
% %     for i = 1:numel(p)
% %         phi = p(i);
% %         mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
% %         mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
% %         mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
% %         mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
% %         mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
% %         mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
% %     
% %         k = 0;
% %         Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
% %         Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
% %         Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
% %     
% %         Pr_o = 0;
% %         Pt_o = 0;
% %         Pf_o = 0;
% %         for k = 0:N2
% %             Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
% %             Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
% %             Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
% %         end
% % 
% %         n = numel(r);
% % %         Btdip = [];
% % %         Btoct = [];
% % %         for k = 1:n
% % %             if r(k)<=Rs
% % %                 R = r(k);
% % %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% % %             elseif r(k)>Rs
% % %                 Btdip(k) =0;
% % %                 Btoct(k)=0;
% % %             end
% % %         end
% % % 
% % %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% %         
% %         Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
% %         Btoct = Boct(j)/(l2+1).*(1./r).^(l2+2)*Pt_o;
% %         Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
% %         Broct =  Boct(j).*(1./r).^(l2+2)*Pr_o;
% %         Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
% %         Bfoct = Boct(j)/(l+1).*(1./r).^(l2+2)*Pf_o;
% %         
% % 
% %         Bt = Btdip + Btoct;
% %         Br = Brdip + Broct;
% %         Bf = Bfdip + Bfoct;
% %     
% %         B = (Bt.^2+Br.^2+Bf.^2).^0.5;
% %         
% %         p_mag = 1/(8*pi).*B.^2;
% %         p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
% %         P = p_mag-p_ram;
% %         while P(n)<0
% %             n = n-1;
% %         end
% %         Rt(i) = c*r(n);
% %     end
% %     plot(p,Rt,'LineWidth',3)
% % end
% % 
% % xlabel('$\phi$','Interpreter','latex')
% % ylabel('$r_t [R_{\star}]$','Interpreter','latex')
% % legend('$B_{oct}/B_{dip}=1$','$B_{oct}/B_{dip}=2$','$B_{oct}/B_{dip}=3$','$B_{oct}/B_{dip}=4$','$B_{oct}/B_{dip}=5$','Interpreter','latex','Location','NorthEastOutside')
% % % title('$R_t$ vs. $B_{oct}/B_{dip}$ with $\beta_{dip}=5$ and $\beta_{oct}=15$','Interpreter','latex')
% 
% beta_d = 5*pi/180;
% tilt_o = (0:5:90)*pi/180;
% Bdip = 1000;
% Boct = [1000;1500;2000;2500;3000];
% F = Boct./Bdip;
% figure(1)
% hold on
% set(gca,'TicklabelInterpreter','latex')
% set(gca,'FontSize',28)
% for j = 1:numel(F)
%     mo = Boct(j)/(l2+1);
%     md = Bdip/(l+1);   
%     r = (1:0.001:30);
%     rr = r*R;
%     n = numel(r);
%     rt = zeros(1,numel(tilt_o));
%     for t = 1:numel(tilt_o)
%         beta_o = tilt_o(t);
%         BB = zeros(1,numel(p));
%         for i = 1:numel(p)
%             phi = p(i);
%             mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
%             mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
%             mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
%             mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
%             mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
%             mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
%         
%             k = 0;
%             Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
%             Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
%             Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
%         
%             Pr_o = 0;
%             Pt_o = 0;
%             Pf_o = 0;
%             for k = 0:N2
%                 Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
%                 Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
%                 Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
%             end
%     
%             n = numel(r);
%     %         Btdip = [];
%     %         Btoct = [];
%     %         for k = 1:n
%     %             if r(k)<=Rs
%     %                 R = r(k);
%     %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%     %             elseif r(k)>Rs
%     %                 Btdip(k) =0;
%     %                 Btoct(k)=0;
%     %             end
%     %         end
%     % 
%     %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%             
%             Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
%             Btoct = Boct(j)/(l2+1).*(1./r).^(l2+2)*Pt_o;
%             Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
%             Broct =  Boct(j).*(1./r).^(l2+2)*Pr_o;
%             Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
%             Bfoct = Boct(j)/(l+1).*(1./r).^(l2+2)*Pf_o;
%             
%     
%             Bt = Btdip + Btoct;
%             Br = Brdip + Broct;
%             Bf = Bfdip + Bfoct;
%         
%             B = (Bt.^2+Br.^2+Bf.^2).^0.5;
%             
%             p_mag = 1/(8*pi).*B.^2;
%             p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
%             P = p_mag-p_ram;
%             while P(n)<0
%                 n = n-1;
%             end
%             BB(i) = B(n)^2;
%         end
%         b = mean(BB);
%         p_mag = 1/(8*pi)*b;
%         a = p_mag-p_ram;
%         idx = find(abs(a)==min(abs(a)));
%         rt(t) = c*r(idx);
%     end
%     
%     plot(rt,tilt_o*180/pi,'LineWidth',3)
% end
% 
% xlabel('$r_t [R_{\star}]$','Interpreter','latex')
% ylabel('$\beta_{oct} [degrees]$','Interpreter','latex')
% legend('$B_{oct}/B_{dip}=1$','$B_{oct}/B_{dip}=2$','$B_{oct}/B_{dip}=3$','$B_{oct}/B_{dip}=4$','$B_{oct}/B_{dip}=5$','Interpreter','latex','Location','NorthEastOutside')
% ylim([0,90])
% 
% % title('$R_t$ vs. $B_{oct}/B_{dip}$ with $\beta_{dip}=5$ and $\beta_{oct}=15$','Interpreter','latex')
% 
% 
% tilt_d = (0:5:90)*pi/180;
% beta_o = 5*pi/180;
% Bdip = 1000;
% Boct = [1000;1500;2000;2500;3000];
% F = Boct./Bdip;
% figure(2)
% hold on
% set(gca,'TicklabelInterpreter','latex')
% set(gca,'FontSize',28)
% for j = 1:numel(F)
%     mo = Boct(j)/(l2+1);
%     md = Bdip/(l+1);   
%     r = (1:0.001:30);
%     rr = r*R;
%     n = numel(r);
%     rt = zeros(1,numel(tilt_d));
%     for t = 1:numel(tilt_d)
%         beta_d = tilt_d(t);
%         BB = zeros(1,numel(p));
%         for i = 1:numel(p)
%             phi = p(i);
%             mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
%             mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
%             mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
%             mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
%             mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
%             mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
%         
%             k = 0;
%             Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
%             Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
%             Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
%         
%             Pr_o = 0;
%             Pt_o = 0;
%             Pf_o = 0;
%             for k = 0:N2
%                 Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
%                 Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
%                 Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
%             end
%     
%             n = numel(r);
%     %         Btdip = [];
%     %         Btoct = [];
%     %         for k = 1:n
%     %             if r(k)<=Rs
%     %                 R = r(k);
%     %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%     %             elseif r(k)>Rs
%     %                 Btdip(k) =0;
%     %                 Btoct(k)=0;
%     %             end
%     %         end
%     % 
%     %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%             
%             Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
%             Btoct = Boct(j)/(l2+1).*(1./r).^(l2+2)*Pt_o;
%             Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
%             Broct =  Boct(j).*(1./r).^(l2+2)*Pr_o;
%             Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
%             Bfoct = Boct(j)/(l+1).*(1./r).^(l2+2)*Pf_o;
%             
%     
%             Bt = Btdip + Btoct;
%             Br = Brdip + Broct;
%             Bf = Bfdip + Bfoct;
%         
%             B = (Bt.^2+Br.^2+Bf.^2).^0.5;
%             
%             p_mag = 1/(8*pi).*B.^2;
%             p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
%             P = p_mag-p_ram;
%             while P(n)<0
%                 n = n-1;
%             end
%             BB(i) = B(n)^2;
%         end
%         b = mean(BB);
%         p_mag = 1/(8*pi)*b;
%         a = p_mag-p_ram;
%         idx = find(abs(a)==min(abs(a)));
%         rt(t) = c*r(idx);
%     end
%     
%     plot(rt,tilt_o*180/pi,'LineWidth',3)
% end
% 
% xlabel('$r_t [R_{\star}]$','Interpreter','latex')
% ylabel('$\beta_{dip} [degrees]$','Interpreter','latex')
% legend('$B_{oct}/B_{dip}=1$','$B_{oct}/B_{dip}=2$','$B_{oct}/B_{dip}=3$','$B_{oct}/B_{dip}=4$','$B_{oct}/B_{dip}=5$','Interpreter','latex','Location','NorthEastOutside')
% ylim([0,90])
% 
% % title('$R_t$ vs. $B_{oct}/B_{dip}$ with $\beta_{dip}=5$ and $\beta_{oct}=15$','Interpreter','latex')
% 
% tilt_d = (0:18:90)*pi/180;
% tilt_o = (0:5:90)*pi/180;
% Bdip = 1200;
% Boct = 1600;
% F = Boct/Bdip;
% figure(3)
% hold on
% set(gca,'TicklabelInterpreter','latex')
% set(gca,'FontSize',28)
% for j = 1:numel(tilt_d)
%     beta_d = tilt_d(j);
%     mo = Boct/(l2+1);
%     md = Bdip/(l+1);   
%     r = (1:0.001:30);
%     rr = r*R;
%     n = numel(r);
%     rt = zeros(1,numel(tilt_o));
%     for t = 1:numel(tilt_o)
%         beta_o = tilt_o(t);
%         BB = zeros(1,numel(p));
%         for i = 1:numel(p)
%             phi = p(i);
%             mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
%             mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
%             mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
%             mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
%             mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
%             mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
%         
%             k = 0;
%             Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
%             Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
%             Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
%         
%             Pr_o = 0;
%             Pt_o = 0;
%             Pf_o = 0;
%             for k = 0:N2
%                 Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
%                 Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
%                 Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
%             end
%     
%             n = numel(r);
%     %         Btdip = [];
%     %         Btoct = [];
%     %         for k = 1:n
%     %             if r(k)<=Rs
%     %                 R = r(k);
%     %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%     %             elseif r(k)>Rs
%     %                 Btdip(k) =0;
%     %                 Btoct(k)=0;
%     %             end
%     %         end
%     % 
%     %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%             
%             Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
%             Btoct = Boct/(l2+1).*(1./r).^(l2+2)*Pt_o;
%             Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
%             Broct =  Boct.*(1./r).^(l2+2)*Pr_o;
%             Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
%             Bfoct = Boct/(l+1).*(1./r).^(l2+2)*Pf_o;
%             
%     
%             Bt = Btdip + Btoct;
%             Br = Brdip + Broct;
%             Bf = Bfdip + Bfoct;
%         
%             B = (Bt.^2+Br.^2+Bf.^2).^0.5;
%             
%             p_mag = 1/(8*pi).*B.^2;
%             p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
%             P = p_mag-p_ram;
%             while P(n)<0
%                 n = n-1;
%             end
%             BB(i) = B(n)^2;
%         end
%         b = mean(BB);
%         p_mag = 1/(8*pi)*b;
%         a = p_mag-p_ram;
%         idx = find(abs(a)==min(abs(a)));
%         rt(t) = c*r(idx);
%     end
%     
%     plot(rt,tilt_o*180/pi,'LineWidth',3)
% end
% 
% xlabel('$r_t [R_{\star}]$','Interpreter','latex')
% ylabel('$\beta_{oct} [degrees]$','Interpreter','latex')
% legend('$\beta_{dip}=0^{\circ}$','$\beta_{dip}=18^{\circ}$','$\beta_{dip}=36^{\circ}$','$\beta_{dip}=54^{\circ}$','$\beta_{dip}=72^{\circ}$','$\beta_{dip}=90^{\circ}$','Interpreter','latex','Location','NorthEastOutside')
% ylim([0,90])
% 
% % title('$R_t$ vs. $B_{oct}/B_{dip}$ with $\beta_{dip}=5$ and $\beta_{oct}=15$','Interpreter','latex')
% 
% 
% tilt_d = (0:5:90)*pi/180;
% tilt_o = (0:18:90)*pi/180;
% Bdip = 1200;
% Boct = 1600;
% F = Boct/Bdip;
% figure(4)
% hold on
% set(gca,'TicklabelInterpreter','latex')
% set(gca,'FontSize',28)
% for j = 1:numel(tilt_o)
%     beta_o = tilt_o(j);
%     mo = Boct/(l2+1);
%     md = Bdip/(l+1);   
%     r = (1:0.001:30);
%     rr = r*R;
%     n = numel(r);
%     rt = zeros(1,numel(tilt_d));
%     for t = 1:numel(tilt_d)
%         beta_d = tilt_d(t);
%         BB = zeros(1,numel(p));
%         for i = 1:numel(p)
%             phi = p(i);
%             mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
%             mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
%             mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
%             mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
%             mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
%             mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
%         
%             k = 0;
%             Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
%             Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
%             Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
%         
%             Pr_o = 0;
%             Pt_o = 0;
%             Pf_o = 0;
%             for k = 0:N2
%                 Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
%                 Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
%                 Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
%             end
%     
%             n = numel(r);
%     %         Btdip = [];
%     %         Btoct = [];
%     %         for k = 1:n
%     %             if r(k)<=Rs
%     %                 R = r(k);
%     %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%     %             elseif r(k)>Rs
%     %                 Btdip(k) =0;
%     %                 Btoct(k)=0;
%     %             end
%     %         end
%     % 
%     %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
%     %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
%             
%             Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
%             Btoct = Boct/(l2+1).*(1./r).^(l2+2)*Pt_o;
%             Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
%             Broct =  Boct.*(1./r).^(l2+2)*Pr_o;
%             Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
%             Bfoct = Boct/(l+1).*(1./r).^(l2+2)*Pf_o;
%             
%     
%             Bt = Btdip + Btoct;
%             Br = Brdip + Broct;
%             Bf = Bfdip + Bfoct;
%         
%             B = (Bt.^2+Br.^2+Bf.^2).^0.5;
%             
%             p_mag = 1/(8*pi).*B.^2;
%             p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
%             P = p_mag-p_ram;
%             while P(n)<0
%                 n = n-1;
%             end
%             BB(i) = B(n)^2;
%         end
%         b = mean(BB);
%         p_mag = 1/(8*pi)*b;
%         a = p_mag-p_ram;
%         idx = find(abs(a)==min(abs(a)));
%         rt(t) = c*r(idx);
%     end
%     
%     plot(rt,tilt_d*180/pi,'LineWidth',3)
% end
% 
% xlabel('$r_t [R_{\star}]$','Interpreter','latex')
% ylabel('$\beta_{dip} [degrees]$','Interpreter','latex')
% ylim([0,90])
% legend('$\beta_{oct}=0^{\circ}$','$\beta_{oct}=18^{\circ}$','$\beta_{oct}=36^{\circ}$','$\beta_{oct}=54^{\circ}$','$\beta_{oct}=72^{\circ}$','$\beta_{oct}=90^{\circ}$','Interpreter','latex','Location','NorthEastOutside')
% % 
% % Bdip = 1000;
% % Boct = 2000;
% % mo = Boct/(l2+1);
% % md = Bdip/(l+1); 
% % beta_d = 5*pi/180;
% % tilt_o = (0:18:90)*pi/180;
% % figure(5)
% % hold on
% % set(gca,'TicklabelInterpreter','latex')
% % set(gca,'FontSize',28)
% % for j = 1:numel(tilt_o)
% %     beta_o = tilt_o(j);
% %       
% %     r = (1:0.001:30);
% %     rr = r*R;
% %     n = numel(r);
% %     Rt = zeros(1,numel(p));
% %     for i = 1:numel(p)
% %         phi = p(i);
% %         mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
% %         mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
% %         mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
% %         mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
% %         mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
% %         mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
% %     
% %         k = 0;
% %         Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
% %         Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
% %         Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
% %     
% %         Pr_o = 0;
% %         Pt_o = 0;
% %         Pf_o = 0;
% %         for k = 0:N2
% %             Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
% %             Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
% %             Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
% %         end
% % 
% %         n = numel(r);
% % %         Btdip = [];
% % %         Btoct = [];
% % %         for k = 1:n
% % %             if r(k)<=Rs
% % %                 R = r(k);
% % %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% % %             elseif r(k)>Rs
% % %                 Btdip(k) =0;
% % %                 Btoct(k)=0;
% % %             end
% % %         end
% % % 
% % %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% %         
% %         Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
% %         Btoct = Boct/(l2+1).*(1./r).^(l2+2)*Pt_o;
% %         Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
% %         Broct =  Boct.*(1./r).^(l2+2)*Pr_o;
% %         Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
% %         Bfoct = Boct/(l+1).*(1./r).^(l2+2)*Pf_o;
% %         
% % 
% %         Bt = Btdip + Btoct;
% %         Br = Brdip + Broct;
% %         Bf = Bfdip + Bfoct;
% %     
% %         B = (Bt.^2+Br.^2+Bf.^2).^0.5;
% %         
% %         p_mag = 1/(8*pi).*B.^2;
% %         p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
% %         P = p_mag-p_ram;
% %         while P(n)<0
% %             n = n-1;
% %         end
% %         Rt(i) = c*r(n);
% %     end
% %     plot(p*180/pi,Rt,'LineWidth',3)
% % end
% % 
% % xlabel('$\phi$','Interpreter','latex')
% % ylabel('$r_t [R_{\star}]$','Interpreter','latex')
% % legend('$\beta_{oct}=0$','$\beta_{oct}=18$','$\beta_{oct}=36$','$\beta_{oct}=54$','$\beta_{oct}=72$','$\beta_{oct}=90$','Interpreter','latex','Location','NorthEastOutside')
% % % title('$R_t$ vs. $B_{oct}/B_{dip}$ with $\beta_{dip}=5$ and $\beta_{oct}=15$','Interpreter','latex')
% 
% % Bdip = 1000;
% % Boct = 2000;
% % mo = Boct/(l2+1);
% % md = Bdip/(l+1); 
% % tilt_d = (0:18:90)*pi/180;
% % beta_o = 5*pi/180;
% % figure(6)
% % hold on
% % set(gca,'TicklabelInterpreter','latex')
% % set(gca,'FontSize',28)
% % for j = 1:numel(tilt_d)
% %     beta_d = tilt_d(j);
% %       
% %     r = (1:0.001:30);
% %     rr = r*R;
% %     n = numel(r);
% %     Rt = zeros(1,numel(p));
% %     for i = 1:numel(p)
% %         phi = p(i);
% %         mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
% %         mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
% %         mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
% %         mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
% %         mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
% %         mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
% %     
% %         k = 0;
% %         Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
% %         Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
% %         Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
% %     
% %         Pr_o = 0;
% %         Pt_o = 0;
% %         Pf_o = 0;
% %         for k = 0:N2
% %             Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
% %             Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
% %             Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
% %         end
% % 
% %         n = numel(r);
% % %         Btdip = [];
% % %         Btoct = [];
% % %         for k = 1:n
% % %             if r(k)<=Rs
% % %                 R = r(k);
% % %                 Btdip(k) = Bdip(i).*(1/R).^(l+2).*Pt_d.*((-R^(2*l+1)+Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %                 Btoct(k) = Boct(j).*(1/R).^(l2+2).*Pt_o.*((-R^(2*l2+1)+Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% % %             elseif r(k)>Rs
% % %                 Btdip(k) =0;
% % %                 Btoct(k)=0;
% % %             end
% % %         end
% % % 
% % %         Brdip =  Bdip(i).*(1./r).^(l+2).*Pr_d.*((l.*r.^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
% % %         Broct =  Boct(j).*(1./r).^(l2+2).*Pr_o.*((l2.*r.^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
% %         
% %         Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
% %         Btoct = Boct/(l2+1).*(1./r).^(l2+2)*Pt_o;
% %         Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
% %         Broct =  Boct.*(1./r).^(l2+2)*Pr_o;
% %         Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
% %         Bfoct = Boct/(l+1).*(1./r).^(l2+2)*Pf_o;
% %         
% % 
% %         Bt = Btdip + Btoct;
% %         Br = Brdip + Broct;
% %         Bf = Bfdip + Bfoct;
% %     
% %         B = (Bt.^2+Br.^2+Bf.^2).^0.5;
% %         
% %         p_mag = 1/(8*pi).*B.^2;
% %         p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
% %         P = p_mag-p_ram;
% %         while P(n)<0
% %             n = n-1;
% %         end
% %         Rt(i) = c*r(n);
% %     end
% %     plot(p*180/pi,Rt,'LineWidth',3)
% % end
% % 
% % xlabel('$\phi$','Interpreter','latex')
% % ylabel('$r_t [R_{\star}]$','Interpreter','latex')
% % legend('$\beta_{dip}=0$','$\beta_{dip}=18$','$\beta_{dip}=36$','$\beta_{dip}=54$','$\beta_{dip}=72$','$\beta_{dip}=90$','Interpreter','latex','Location','NorthEastOutside')
% % % title('$R_t$ vs. $B_{oct}/B_{dip}$ with $\beta_{dip}=5$ and $\beta_{oct}=15$','Interpreter','latex')
% 
% % 
beta_d = 20*pi/180;
beta_o = 10*pi/180;
Bdip = 1200;
Boct = 1600;
theta = pi/2;
p = 0:0.1:2*pi;

mo = Boct/(l2+1);
md = Bdip/(l+1);   
r = (1:0.001:30);
rr = r*Ra;
n = numel(r);

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
            
            BB(i) = B(n)^2;
end
b = mean(BB);
p_mag = 1/(8*pi)*b;
a = p_mag-p_ram;
idx = find(abs(a)==min(abs(a)));
rt = c*r(idx);


% 
% phase_d = 0.00;
% phase_o = 0.0;
% psi_d = (1-phase_d)*2*pi;
% psi_o = (1-phase_o)*2*pi;
% tilt_d = 5*pi/180;
% tilt_o = 20*pi/180;
% Bd = [1000;1500;2000;2500;3000];
% theta = pi/2;
% p = 0:0.1:2*pi;
% 
% figure(4)
% hold on
% set(gca,'TicklabelInterpreter','latex')
% set(gca,'FontSize',28)
% for j = 1:numel(Bd)
%     beta_o = tilt_o;
%     Bdip = Bd(j);
%     Bo = Bdip:100:Bdip*3;
%     F = Bo./Bdip;
%     md = Bdip/(l+1);   
%     r = (1:0.001:30);
%     rr = r*R;
%     n = numel(r);
%     rt = zeros(1,numel(F));
%     for t = 1:numel(F)
%         Boct = Bo(t);
%         mo = Boct/(l2+1);
%         beta_d = tilt_d;
%         BB = zeros(1,numel(p));
%         for i = 1:numel(p)
%             phi = p(i);
%             mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
%             mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
%             mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
%             mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
%             mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
%             mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
%         
%             k = 0;
%             Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
%             Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
%             Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
%         
%             Pr_o = 0;
%             Pt_o = 0;
%             Pf_o = 0;
%             for k = 0:N2
%                 Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
%                 Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
%                 Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
%             end
%     
%             n = numel(r);
%     %      
%             Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
%             Btoct = Boct/(l2+1).*(1./r).^(l2+2)*Pt_o;
%             Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
%             Broct =  Boct.*(1./r).^(l2+2)*Pr_o;
%             Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
%             Bfoct = Boct/(l+1).*(1./r).^(l2+2)*Pf_o;
%             
%     
%             Bt = Btdip + Btoct;
%             Br = Brdip + Broct;
%             Bf = Bfdip + Bfoct;
%         
%             B = (Bt.^2+Br.^2+Bf.^2).^0.5;
%             
%             p_mag = 1/(8*pi).*B.^2;
%             p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
%             P = p_mag-p_ram;
%             while P(n)<0
%                 n = n-1;
%             end
%             BB(i) = B(n)^2;
%         end
%         b = mean(BB);
%         p_mag = 1/(8*pi)*b;
%         a = p_mag-p_ram;
%         idx = find(abs(a)==min(abs(a)));
%         rt(t) = c*r(idx);
%     end
%     
%     plot(rt,F,'LineWidth',3)
% end
% 
% xlabel('$r_t [R_{\star}]$','Interpreter','latex')
% ylabel('$B_{oct}/B_{dip}$','Interpreter','latex')
% % legend('$\beta_{oct}=0^{\circ}$','$\beta_{oct}=18^{\circ}$','$\beta_{oct}=36^{\circ}$','$\beta_{oct}=54^{\circ}$','$\beta_{oct}=72^{\circ}$','$\beta_{oct}=90^{\circ}$','Interpreter','latex','Location','NorthEastOutside')


% tilt_d = 5*pi/180;
% beta_o = 20*pi/180;
% Bdip = 1200;
% Boct = 1600;
% phase_dd = [0.0;0.25;0.5;0.75];
% phase_oo = 0:0.1:1;
% 
% F = Boct./Bdip;
% figure(2)
% hold on
% set(gca,'TicklabelInterpreter','latex')
% set(gca,'FontSize',28)
% for j = 1:numel(phase_dd)
%     phase_d = phase_dd(j);
%     psi_d = (1-phase_d)*2*pi;
%     mo = Boct/(l2+1);
%     md = Bdip/(l+1);   
%     r = (1:0.001:30);
%     rr = r*R;
%     n = numel(r);
%     rt = zeros(1,numel(phase_oo));
%     for t = 1:numel(phase_oo)
%         phase_o = phase_oo(t);
%         psi_o = (1-phase_o)*2*pi;
%         BB = zeros(1,numel(p));
%         for i = 1:numel(p)
%             phi = p(i);
%             mr_d = md*sin(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*sin(theta)*sin(phi)*sin(psi_d)*sin(beta_d)+md*cos(theta)*cos(beta_d);
%             mt_d = md*cos(theta)*cos(phi)*cos(psi_d)*sin(beta_d)+md*cos(theta)*sin(phi)*sin(psi_d)*sin(beta_d)-md*sin(theta)*cos(beta_d);
%             mf_d = -md*sin(phi)*cos(psi_d)*sin(beta_d)+md*cos(phi)*sin(psi_d)*sin(beta_d);
%             mr_o = mo*sin(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*sin(theta)*sin(phi)*sin(psi_o)*sin(beta_o)+mo*cos(theta)*cos(beta_o);
%             mt_o = mo*cos(theta)*cos(phi)*cos(psi_o)*sin(beta_o)+mo*cos(theta)*sin(phi)*sin(psi_o)*sin(beta_o)-mo*sin(theta)*cos(beta_o);
%             mf_o = -mo*sin(phi)*cos(psi_o)*sin(beta_o)+mo*cos(phi)*sin(psi_o)*sin(beta_o);
%         
%             k = 0;
%             Pr_d = ((-1)^k*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k))*(mr_d/md)^(l-2*k);
%             Pt_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mt_d/md;
%             Pf_d = ((-1)^(k+1)*factorial(2*l-2*k))/(2^l*factorial(k)*factorial(l-k)*factorial(l-2*k-1))*(mr_d/md)^(l-2*k-1)*mf_d/md;
%         
%             Pr_o = 0;
%             Pt_o = 0;
%             Pf_o = 0;
%             for k = 0:N2
%                 Pr_o = Pr_o+((-1)^k*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k))*(mr_o/mo)^(l2-2*k);
%                 Pt_o = Pt_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mt_o/mo;
%                 Pf_o = Pf_o+((-1)^(k+1)*factorial(2*l2-2*k))/(2^l2*factorial(k)*factorial(l2-k)*factorial(l2-2*k-1))*(mr_o/mo)^(l2-2*k-1)*mf_o/mo;
%             end
%     
%            
%             Btdip = Bdip/(l+1).*(1./r).^(l+2)*Pt_d;
%             Btoct = Boct/(l2+1).*(1./r).^(l2+2)*Pt_o;
%             Brdip =  Bdip.*(1./r).^(l+2)*Pr_d;
%             Broct =  Boct.*(1./r).^(l2+2)*Pr_o;
%             Bfdip = Bdip/(l+1).*(1./r).^(l+2)*Pf_d;
%             Bfoct = Boct/(l+1).*(1./r).^(l2+2)*Pf_o;
%             
%     
%             Bt = Btdip + Btoct;
%             Br = Brdip + Broct;
%             Bf = Bfdip + Bfoct;
%         
%             B = (Bt.^2+Br.^2+Bf.^2).^0.5;
%             
%             p_mag = 1/(8*pi).*B.^2;
%             p_ram = 1/(8*pi).*(2*G*M)^0.5.*M_acc.*rr.^(-5/2);
%             P = p_mag-p_ram;
%             while P(n)<0
%                 n = n-1;
%             end
%             BB(i) = B(n)^2;
%         end
%         b = mean(BB);
%         p_mag = 1/(8*pi)*b;
%         a = p_mag-p_ram;
%         idx = find(abs(a)==min(abs(a)));
%         rt(t) = c*r(idx);
%     end
%     
%     plot(rt,(1-phase_oo)*360,'LineWidth',3)
% end
% 
% xlabel('$r_t [R_{\star}]$','Interpreter','latex')
% ylabel('$\psi_{dip} [degrees]$','Interpreter','latex')
% ylim([0,360])
