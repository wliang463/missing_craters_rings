function [tempf,depth]= tneal()

%50m resolution
%1km lava, not held at constant temp
%boundary is held at 250K
%top is all lava

dz = 50;%1e3;
dx = 50;%1.98e3;

load tycho.mat
z1 = -(z0-(max(z0)+1e3))/1e3;%-(z0/1e3-3);
t_bound = [repmat(z1(1),1,5) z1 repmat(z1(end),1,5)];
t1 = interp1((0:70)*1.98e3,t_bound,(50:50:140000));
t1(isnan(t1)) = [];

num = 1e3/dz;

sz2 = size(t1,2);
sz1 = num*10;%400;

midp = 1307;

grid0 = ones(sz1,sz2);%sz2 is diameter in pixel of tycho
x = repmat((1:sz2)-midp,sz1,1)*dx;%1.98e3;%repmat(-round(sz2/2)+1:round(sz2/2),sz1,1)*1.98e3;

top_ind = round((t1-1)*num+1);
temp = grid0;

above_row = [];
above_col = [];

hla = 1;

depths = zeros(1,sz2);

for i=1:sz2
    add = (2:top_ind(i)+num*hla-1);
    if isempty(add)
        continue
    end
    above_row = [above_row add];
    above_col = [above_col repmat(i,1,numel(add))];
end
lava = sub2ind(size(grid0),above_row,above_col);

temp(lava) = 1423;
%tind3 = sub2ind(size(grid0),t_bound+2,1:sz2);temp(tind3) = 1423;

temp(temp~=1423) = 250;

temp_bound = 250;
tind = sub2ind(size(grid0),ones(1,sz2),1:sz2);temp(tind) = temp_bound;
%tind = sub2ind(size(grid0),top_ind,1:sz2);temp(tind) = temp_bound;


%boundary = [sub2ind(size(grid0),repmat((1:600),1,2),[repmat(1,1,600) repmat(sz2,1,600)]) ...
%    sub2ind(size(grid0),[repmat(1,1,sz2) repmat(600,1,sz2)],repmat((1:sz2),1,2))];

boundary = [sub2ind(size(grid0),[1 2 1 2],repmat(sz2,1,4)) ...
    sub2ind(size(grid0),[repmat(1,1,sz2) repmat(sz1,1,sz2)],repmat((1:sz2),1,2))];

above_row = [];
above_col = [];
for i=1:sz2
    add = (1:top_ind(i)-1);
    if isempty(add)
        continue
    end
    above_row = [above_row add];
    above_col = [above_col repmat(i,1,add(end))];
end

above = sub2ind(size(grid0),above_row,above_col);

no_change = unique([boundary tind]);%above
no_change(no_change==554202) = [];
no_change(no_change==554203) = [];


q = 0.0916;%0.1;
restt = sz1-(min(top_ind)+num*hla)+1;
temp_down = 250+q/2.3*dz*(1:restt);
temp(min(top_ind)+num*hla:end,1:end) = repmat(temp_down',1,sz2);
temp(lava) = 1423;

load('poro_dep.mat')
poro_down = interp1(depth0,poro0,dz*(1:restt)/1e3,'spline');
poro = grid0;
poro(min(top_ind)+num*hla:end,1:end) = repmat(poro_down',1,sz2);
poro(lava) = 0.15;
poro(tind) = 0.15;%boundary


rho = 2800;
mol = 0.27;%molar mass
c = 1100;%270/mol;%original cp in j/(mol*K)


t0 = 1000000;

alpha = 2.3/(rho*c);
dt = 3.154e7;%(dx)^2/alpha*1e-5*0.66;
%dt = 0.01;%0.5*(1737e3*dx*cosd(lats)).^2/(2*D);
%dt = 0.5*(1737e3*dx*cosd(lats)).^2/(2*D);
tic

temp0 = temp;

%{
for i=1:sz2
    ind = find(temp0(:,i)>923);
    depths0(i) = max(ind)*dz;%-z1(1)*1e3;
end
%}

inds = find(temp(:)>923);
[ii,jj] = ind2sub(size(temp),inds);
max0 = max(ii(ii<140));

annealed = temp>923;

%k = grid0*2.3;%*3.15e7;
%k = 100./(7.4+0.05*temp);
%k(temp>500) = k(temp>500) + 0.0023*(temp(temp>500)-500);

m = -log(10)/0.2;
%k = m*poro+k;

lava_grid = grid0-1;
lava_grid(lava) = 1;

for i=1:t0%linspace(1,num,steps)
    
    k = 100./(7.4+0.05*temp);
    k(temp>500) = k(temp>500) + 0.0023*(temp(temp>500)-500);
    
    poro(~lava_grid & temp>923) = 0;
    k = exp(m*poro).*k;
    %kmin = min(k(:))
%    k = grid0*0.13;%2.3;
%    k(lava) = 1;
    if i==222
%        figure;imagesc(k);colorbar
    end

    down_1 = circshift(temp,-1,1);down1 = circshift(temp,1,1);
    side_1 = circshift(temp,-1,2);side1 = circshift(temp,1,2);

    kdown_1 = circshift(k,-1,1);kdown1 = circshift(k,1,1);
    kside_1 = circshift(k,-1,2);kside1 = circshift(k,1,2);

    dk_dz = (kdown_1-k*2+kdown1)/dz;%(temp-down1)/dz;
    dk_dx = (kside_1-k*2+kside1)/dx;%(temp-down1)/dz;
    
    dtemp_dz = (down1-2*temp+down_1)/dz;%(temp-down1)/dz;
    dtemp2_dz2 = (down1-2*temp+down_1)/(dz)^2;
    %dtemp2_dz2(1,:) = 0;dtemp2_dz2(end,:) = 0;dtemp_dz(end,:) = 0;
    
    dtemp_dx = (side1-2*temp+side_1)/dx;%(temp-down1)/dz;
    dtemp2_dx2 = (side1-2*temp+side_1)/(dx)^2;  

    grad_temp = dtemp_dx + dtemp_dz;
    lap_temp = dtemp2_dx2 + dtemp2_dz2 + 1./x.*dtemp_dx;
    lap_temp(:,midp) = dtemp2_dz2(:,midp);

    grad_k = dk_dx + dk_dz;

    prod1 = 0;%dtemp_dx.*dk_dx + dtemp_dz.*dk_dz;%0;%dtemp_dx.*dk_dx + 
    
    prod2 = k.*(lap_temp);
    q_change = 0;
    
    dtemp_dt = (prod1 + prod2)/(rho*c);

    dtemp_dt(no_change) = 0;
    %(dz_dx.*dx.^(-1) + dz_dy.*dy.^(-1))*D;x`
    
   % dt = 1e3;
    %dt = 0.9*(1737e3*dx*cosd(lats)).^2/(2*D);
    
    temp = temp + dtemp_dt.*dt;  

%    disp([i temp(94,1303)])

    intermed = annealed+(temp>923);
    annealed(intermed>=1) = 1;
    %annealed = annealed || (temp>923);
    
    if rem(i,1000) == 0%sum(sum(isnan(temp))) > 0
        disp(i)
        toc
        inds = find(temp(:)>923);
        [ii,jj] = ind2sub(size(temp),inds);
        maxf = max(ii(ii<140));
        depth = (maxf-max0)*dz;
        disp(depth)

    end
    
    if rem(i,1e4) == 0
        save(['tempf_' num2str(i) '.mat'],'temp');
        toc
    end
    
    if 0%sum(num2str(i)-'0') == 1
        disp(i)
        toc
        save(['tempf_' num2str(i) '.mat'],'temp');
    end
    %{
    imagesc(theta_portion);axis equal;colorbar
    saveas(gcf,['diffusion' num2str(i) '.png'])
    imagesc(phi_portion);axis equal;colorbar
    saveas(gcf,['diffusion' num2str(i) '_phi.png'])
    
    %}

    
    if 0%rem(i,1e4) == 0
        depths = [];
        if max(temp(:)) < 923
            disp('brrrrrr')
            depths = 0;
        else
            for j=1:sz2
                ind = find(temp(:,j)>923);
                ind(ind>329) = [];
                if ~isempty(ind)% && ind > depths(j)
                    depths(j) = max(ind)*dz;%-z1(1)*1e3;
                end
            end
        end
        depth = max(depths)-max(depths0);%median(depths);
        disp(depth)

    end
end

tempf = temp;


inds = find(tempf(:)>923);
[ii,jj] = ind2sub(size(tempf),inds);
maxf = max(ii(ii<140));

if max(tempf(:)) < 923
    disp('brrrrrr')
    depths = 0;
else
    for i=1:sz2
        ind = find(tempf(:,i)>923);
        ind(ind>329) = [];
        if ~isempty(ind)
            depths(i) = max(ind)*dz;%-z1(1)*1e3;
        end
    end
end
%}


%depth = max(depths)-max(depths0);%median(depths);
depth = (maxf-max0)*dz;
%}
%depth = 0;
xtt = arrayfun(@num2str, dx*0.5*(1:5), 'UniformOutput', 0);
ytt = arrayfun(@num2str, dz*0.02*(1:20), 'UniformOutput', 0);%dz*0.02

figure;imagesc(temp);colorbar
xlabel('Distance (km)','FontSize',40)
xticklabels(xtt)
ylabel('Depth (km)','FontSize',40);
yticklabels(ytt)
title('Temperature (K) after 30000 yr')
ax = gca;
ax.FontSize = 24; 
set(gcf,'color','w')
set(ax, 'YTick', 20:20:200)
caxis([251.9913 1423])


save depths.mat depths
save annealed.mat annealed

%figure;plot((1:sz2),-depths0);hold on;plot((1:sz2),-depths)




save tempf.mat tempf
