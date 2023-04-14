clear all
clc
clf
set(gca,'TicklabelInterpreter','latex')
set(gca,'FontSize',26)
l = 1;
l2 = 3;
Bdip = 1700;      
Boct = 500;

Rs =7.6;
R_max = [1.4];  
% ;2.9:3.1;3.3;3.5;3.7;3.9;4.1

th = 0:pi/20:pi;
Phi = 0:2*pi/20:2*pi;
ds = 0.01;

phase_d = 0.65;
phase_o = 0.5;
psi_d = (1-phase_d)*2*pi;
psi_o = (1-phase_o)*2*pi;
beta_d = 170*pi/180;
beta_o = 10*pi/180;
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
ax = gca;
fig = gcf;
hold on

set(ax,'FontSize',26)
axis equal
for p = 1:numel(Phi)
    for t=1:numel(th)
        Rmax = R_max;
        for i = 1:numel(Rmax)
            theta = th(t);
            R = Rmax(i);
            rloop =[R];
            phi = Phi(p);
            theta_loop = [theta];
            phi_loop = [phi];
            while R>1 && R<=Rs
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
            plot3(ax,Rx,Ry,Rz,'black','linewidth',2)    
        end
    end
end

for p = 1:numel(Phi)
    for t=1:numel(th)
        Rmax = R_max;
        for i = 1:numel(Rmax)
            theta = th(t);
            R = Rmax(i);
            rloop =[R];
            phi = Phi(p);
            theta_loop = [theta];
            phi_loop = [phi];
            while R>1 && R<=Rs
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
            plot3(ax,Rx,Ry,Rz,'black','linewidth',2)    
        end
    end
end


theta = 0:pi/100:pi;
p = pi:2*pi/100:3*pi;
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
% [x,y,z] = sphere;
% s1 = surf(x,y,z,"FaceColor",[1 1 0]);
[x,y,z] = sphere(100);
s1 = surf(ax,x,y,z,flipud(B)/1000,'EdgeColor','none');
colormap("turbo")
colorbar;
xlabel('Rx[ $R_{\star}$]','Interpreter','latex')
zlabel('Rz [$R_{\star}$]','Interpreter','latex')
ylabel('Rz [$R_{\star}$]','Interpreter','latex')
a = colorbar;
a.TickLabelInterpreter = 'latex';
a.Label.String = 'kG';
set(a.XLabel,{'String','Rotation','Position'},{'kG',0,[0.9 3.6]},'Interpreter','latex')
xlabel('x [$R_{\star}$]','Interpreter','latex')
zlabel('z [$R_{\star}$]','Interpreter','latex')
ylabel('y [$R_{\star}$]','Interpreter','latex')
xlim([-3,3])
zlim([-2,2])
hold off

% fig.Position = [100 100 1000 1000];
% fig = gcf;
% plotly_fig = fig2plotly(fig, 'offline', false,'strip',false);
% plotly_url = 'https://chart-studio.plotly.com';
% signin('gborghetto', 'vNzIpavF8Gl8rHaL0hVf', 'PlotlyURL', plotly_url);
% export_fig('myplot.u3d','-u3d')
% fig = get(groot,'CurrentFigure');
% verLessThan('matlab','8.4.0')
% vrml(gcf,'try') 

% fig2u3d(ax, 'test', '-pdf')
% copyfile('test.u3d', '..\tex\personal\3dheart\img\test.u3d')

% fig = gcf;
% plotly_fig = fig2plotly(fig);
% html_str = plotly_fig.url;
% json_str = savejson(html_str);

theta = (0:120/500:120)*pi/180;
p = 0:2*pi/500:2*pi;
B1 = zeros(numel(theta),numel(p));
B2 = zeros(numel(theta),numel(p));
B3 = zeros(numel(theta),numel(p));

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
    Btdip = Bdip/(l+1)*(1/R)^(l+2)*Pt_d*((-(l+1)*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
    Btoct = Boct/(l2+1)*(1/R)^(l2+2)*Pt_o*((-(l2+1)*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
    Bfdip = Bdip/(l+1)*(1/R)^(l+2)*Pf_d*((-(l+1)*R^(2*l+1)+(l+1)*Rs^(2*l+1))/(l+(l+1)*Rs^(2*l+1)));
    Bfoct = Boct/(l2+1)*(1/R)^(l2+2)*Pf_o*((-(l2+1)*R^(2*l2+1)+(l2+1)*Rs^(2*l2+1))/(l2+(l2+1)*Rs^(2*l2+1)));
    
    B1(:,i) = Brdip + Broct;
    B2(:,i) = Btdip + Btoct;
    B3(:,i) = Bfdip + Bfoct;

end
pos = [0 30 60 90 120];
% Rticks = {'0^{\circ}','30^{\circ}','60^{\circ}','90^{\circ}','100^{\circ}','Interpreter','latex'};
RR = linspace(0,120,501); % (distance in km)
Az = linspace(0,360,501); % in degrees
P = B1/1000;
P1 = rescale(B2/1000,min(P,[],'all'),max(P,[],'all'));
figure
set(gca,'FontSize',26)
subplot(1,3,1)
[~,c]=polarPcolor(RR,Az,P,'Nspokes',5,'circlesPos',pos,'ncolor',10);
caxis([-2.2,2.2])
% caxis([-max(P,[],'all'),max(P,[],'all')])

subplot(1,3,3)
[~,b]=polarPcolor(RR,Az,B2/1000,'Nspokes',5,'circlesPos',pos,'ncolor',10);
% caxis([-max(P,[],'all'),max(P,[],'all')])
caxis([-2.2,2.2])

subplot(1,3,2)
[~,a]=polarPcolor(RR,Az,B3/1000,'Nspokes',5,'circlesPos',pos,'ncolor',10);
% ylabel(c,'kG','Interpreter','latex','FontSize',26);
colormap("turbo")
% caxis([-max(P,[],'all'),max(P,[],'all')])
caxis([-2.2,2.2])


c.TickLabelInterpreter = 'latex';
% set(c,'FontSize',26);

