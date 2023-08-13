function [counts,rads] = crater_search(excel,eig,grav)

%excel: crater catalog
%eig: Bouguer gravity gradient
%grav: Bouguer gravity

%counts: number of detected craters
%rads: radius of detected craters

eig(1,:) = [];
eig(:,1) = []; 

global siz2 siz1

siz2 = size(eig,2);
siz1 = size(eig,1);

side = 'far';

eig = fliplr(eig);
grav = fliplr(grav);

siz = size(excel);%locs = round(locs);

excel(:,2) = excel(:,2)-0.1795;
excel(:,1) = excel(:,1)+0.1795;

edges = 20:20:200;
Y = discretize(excel(:,3),edges);
Y(isnan(Y)) = numel(edges);

total = zeros(numel(edges),1);
good_bin = zeros(numel(edges),1);


e_arr = (0:0.005:0.4);

radius = [];
bad = 0;
skipped = 0;

y00 = repmat((1:siz1)',1,siz2);
x00 = repmat((1:siz2),siz1,1);

tic

data = [];
all_c = [];

load('masks_erasure.mat');
bw_f = circshift(bw_f,-501,2);
bw_n = circshift(bw_n,-501,2);

data_f = [];
for i=1:siz(1)
    
    if rem(i,1000)==1
        toc
    end

    switch side
        case 'far'
            cond = abs(excel(i,1)) < 90;
        case 'near'
            cond = abs(excel(i,1)) > 90;
    end
    
    if cond
        rads(i) = nan;
        skipped = skipped + 1;
        continue
    end

    data_f = [data_f excel(i,3)];
    
    load('mare_mask.mat');

    
    x0 = excel(i,1);y0 = excel(i,2);
    
    if x0 > 0 && abs(x0) > 90
        x0 = round((x0-90)/360*siz2);
    else
        x0 = round((x0+270)/360*siz2);
    end
    
    y0 = round((90-y0)/180*siz1);
    
    if x0 == 0
        x0 = 1;
    elseif y0 == 0
        y0 = 1;
    end
    
    ind1 = sub2ind(size(eig),y0,x0);
    
    switch side
        case 'far'
            cond = ~bw_f(ind1);
        case 'near'
            cond = ~bw_n(ind1);
    end
    if cond
        rads(i) = nan;
        skipped = skipped + 1;
        continue
    end
     %} 
    
    all_c = [all_c excel(i,3)];
    
    a = linspace(excel(i,3)/2*0.7,excel(i,3)/2*1.3,50);
    a = a(randperm(numel(a)));
    %a = excel(i,3)/2;
    vals = [];
    
    

    for j=1:numel(a)

        [xx,yy] = get_xy(excel(i,1),excel(i,2),a(j)*cosd(excel(i,2)));
        if isnan(sum(xx)) || isnan(sum(yy))
            continue
        end
        inds = sub2ind(size(eig),yy,xx);
        %vals(j,:)= eig(inds);
        vals = eig(inds);
        
        eig_in = eigin(xx,yy,eig);
        rim_f = rimf(vals);

        [xx,yy] = get_xy(excel(i,1),excel(i,2),a(j)*cosd(excel(i,2))*0.7);
        grav_in = eigin(xx,yy,grav);

        [xx,yy] = get_xy(excel(i,1),excel(i,2),a(j)*cosd(excel(i,2))*1.3);
        inds2 = sub2ind(size(grav),yy,xx);
        vals2 = grav(inds2);

        
        
        if x0 < 1250 && x0 > 1230 && y0 < 390 && y0 > 375
                disp(std(grav_in)*1e5)
        end
        
        if (max(rim_f)<12 || quantile(eig_in,0.5) > 0) && j==numel(a) && (abs(mean(grav_in)-mean(vals2)) < 25/1e5 || std(grav_in) > 15/1e5) 
            if x0 == 363 && y0 == 154
                aaa = 1;
            end
            bad = bad+1;
        elseif (abs(mean(grav_in)-mean(vals2)) > 25/1e5 && std(grav_in) < 15/1e5) || (max(rim_f)>12 && quantile(eig_in,0.5) < 0) %mean(grav_in) > mean(vals2)+25/1e5
            rads(i) = a(j);
            data = [data excel(i,3)];

            if x0 < 525 && x0 > 521 && y0 < 464 && y0 > 460
                %523.1049  462.2861
                %disp(std(grav_in)*1e5)
            end
        
            
            break
        elseif ((abs(mean(grav_in)-mean(vals2)) > 15/1e5 && std(grav_in) < 35/1e5) || (max(rim_f)>12 && quantile(eig_in,0.5) < 0)) && excel(i,3) > 140
            
            rads(i) = a(j);
            data = [data excel(i,3)];
            
            break
        else

            continue
        
        end
       

    end

end

[counts, bins] = histcounts(data,edges);
[counts_a, bins] = histcounts(all_c,edges);

figure;plot(edges(1:9),counts./counts_a)

r_p = counts./counts_a;

cdf = cumsum(counts,'reverse');
cdf = cdf/max(cdf);

good = numel(data);

save crater_search_far_flip.mat good cdf r_p counts rads data

end


function [xx,yy] = get_xy(x0,y0,a)

    global siz2 siz1

    xx = [];
    yy = [];
    lat = y0;
  
    dlong = a/(1737*cosd(lat));
    r = dlong/(2*pi)*siz2;
    
    if r < 1
        xx = nan;
        yy = nan;
        return;
    end
    
    if x0 > 0 && abs(x0) > 90
        x0 = round((x0-90)/360*siz2);
    else
        x0 = round((x0+270)/360*siz2);
    end
    
    y0 = round((90-y0)/180*siz1);

    if x0 < 1250 && x0 > 1230 && y0 < 390 && y0 > 375
            x0 = x0-1;
            y0 = y0-2;
    end

    if x0 < 1214 && x0 > 1206 && y0 < 428 && y0 > 418
            x0 = x0-1;
            y0 = y0-2;
    end
    
    for theta = 1:360
        xx = [xx x0+r*cosd(theta)/abs(cosd(lat))];
        yy = [yy y0+r*sind(theta)];
    end
    
    xx = round(xx); yy = round(yy);
    
    yy0 = yy;
    xx0 = xx;
    
    xx(xx0>-siz2 & xx0<1) = xx(xx0>-siz2 & xx0<1)+siz2;
    xx(xx0<=-siz2) = 1;

    yy(yy0>siz1) = siz1;
    yy(yy0<1) = 1;
    
    xx0 = xx;

    xx(xx0>siz2 & xx0<2*siz2) = xx(xx0>siz2 & xx0<2*siz2)-siz2;
    xx(xx0>=2*siz2) = siz2;

end

function eig_in = eigin(xx,yy,eig)
    y_list = [];
    x_list = [];
    for j=min(xx):max(xx)
        yy_bound = sort(unique(yy(xx==j)));
        
        if isempty(yy_bound)
            continue
        end
        
        y_list = [y_list yy_bound(1)+1:yy_bound(end)-1];
        x_list = [x_list repmat(j,numel(yy_bound(1)+1:yy_bound(end)-1),1)'];
        
    end
    in = sub2ind(size(eig),y_list,x_list);
    
    in0 = sub2ind(size(eig),yy,xx);
    
    
    eig_in = unique(eig(setdiff(in,in0)));
end

function rim_f = rimf(rim)

    rimm = [rim rim];
    grid = round(numel(rim)*0.3);%continuous measurement of 0.3 of the rim
    rim_a = movmin(rimm,grid);
    rim_f = rim_a(grid-1:numel(rim)+grid-2);
end

