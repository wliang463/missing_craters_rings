function topof = diffusion(topo)

%topo: topography of the Moon

siz1 = size(topo,1);
siz2 = size(topo,2);

topo0 = topo;

steps = 1e3;
num = 1e3;
lats = flipud(repmat(linspace(-89.9,89.9,siz1)',1,siz2));
theta = repmat(linspace(0.1,179.9,siz1)',1,siz2);


lats = flipud(repmat(linspace(-90,90,siz1)',1,siz2));
theta = repmat(linspace(0,180,siz1)',1,siz2);

theta = theta(8:1500,:);
topo = topo(8:1500,:);

D = 1;

dtheta = pi/siz1;
dphi = 2*pi/siz2;

t0 = 1e6;

dt = 1e2;
tic

for i=1:1e7
    down_1 = circshift(topo,-1,1);down1 = circshift(topo,1,1);
    side_1 = circshift(topo,-1,2);side1 = circshift(topo,1,2);
    
    dz_dtheta = (down1-2*topo+down_1)/dtheta;%(topo-down1)/dtheta;
    dz2_dtheta2 = (down1-2*topo+down_1)/(dtheta)^2;

    dz_dtheta(1,:) = (-topo(1,:)+down_1(1,:))/dtheta;
    dz_dtheta(end,:) = (down1(end,:)-topo(end,:))/dtheta;
    dz2_dtheta2(1,:) = (-topo(1,:)+down_1(1,:))/(dtheta)^2;
    dz2_dtheta2(end,:) = (down1(end,:)-topo(end,:))/(dtheta)^2;
    
    dz2_dphi2 = (side1-2*topo+side_1)/(dphi)^2;  
    
    theta_portion = ((topo).^2.*sind(theta)).^(-1).*((cosd(theta).*dz_dtheta)+sind(theta).*dz2_dtheta2);
    phi_portion = (topo.*sind(theta)).^(-2).*dz2_dphi2;
    r_portion = 2*(topo).^(-1);
    
    dz_dt = (theta_portion + phi_portion + r_portion)*D;
    
    topo = topo + dz_dt.*dt;    
    
    if 0%rem(i,1000) == 0%sum(sum(isnan(topo))) > 0
        disp(i)
        toc
    end
    
    if rem(i,1e5) == 0
        topo0(8:1500,:) = topo;
        topo0(1:7,:) = 1737e3;
        topo0(1501:end,:) = 1737e3;
        save(['topof_' num2str(i) '.mat'],'topo');
        %disp('saved')
        disp(i)
        toc
    end
    
    if sum(num2str(i)-'0') == 1
        disp(i)
        toc
        save(['topof_' num2str(i) '_check.mat'],'topo');
    end
    
    
end

topo0(236:1770,:) = topo;
topo0(1:235,:) = 1737e3;
topo0(1771:end,:) = 1737e3;


topof = topo0;

save topof.mat topof
