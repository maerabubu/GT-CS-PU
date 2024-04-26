function [ind,xy0,xy,grid_ij,nzix,ph0,ph] = CtSent_downsample4PU(ph,xy_input,pix_size,poly)

fprintf('   Resampling phase to grid...\n')

if nargin<2
    error('not enough arguments')
end
if nargin<3
    pix_size = 250;
end
if nargin<4
    poly = [];
end


if ~isreal(ph) && sum(ph(:)==0)>0
    error('Containing Zero Phase Values!')
end
ind = inpolygon(xy_input(:,1),xy_input(:,2),poly(:,1),poly(:,2));
xy0 = xy_input(ind,:);
ph0 = ph(ind,:);
xy = xy_input(~ind,:);
phsbas = ph(~ind,:);
[n_ps,n_ifg]=size(phsbas);

fprintf('   Number of interferograms  : %d\n',n_ifg)
fprintf('   Number of points per ifg  : %d\n',n_ps)
fprintf('   Number of non-resampled points: %d\n',size(xy0,1))
% Define the resampled grid size
if pix_size==0
    grid_x_min = 1;
    grid_y_min = 1;
    n_i = max(xy(:,2));
    n_j = max(xy(:,1));
    grid_ij = [xy(:,2),xy(:,1)];
else
    grid_x_min = min(xy(:,1));
    grid_y_min = min(xy(:,2));

    grid_ij(:,1)=ceil((xy(:,2)-grid_y_min+1e-3)/pix_size);
    grid_ij(grid_ij(:,1)==max(grid_ij(:,1)),1)=max(grid_ij(:,1))-1;
    grid_ij(:,2)=ceil((xy(:,1)-grid_x_min+1e-3)/pix_size);
    grid_ij(grid_ij(:,2)==max(grid_ij(:,2)),2)=max(grid_ij(:,2))-1;

    n_i = max(grid_ij(:,1));
    n_j = max(grid_ij(:,2));
end

ph_grid = zeros(n_i,n_j,'single');
 
for i1=1:n_ifg
    if isreal(phsbas)
        ph_this = exp(1i*phsbas(:,i1));
    else
        ph_this = phsbas(:,i1);
    end   
    ph_grid(:) = 0;
    if pix_size==0
        ph_grid((xy(:,1)-1)*n_i+xy(:,2)) = ph_this;
    else
        for i=1:n_ps     
            ph_grid(grid_ij(i,1),grid_ij(i,2)) = ph_grid(grid_ij(i,1),grid_ij(i,2))+ph_this(i);
        end
    end
  
    if i1==1
        nzix=ph_grid~=0;
        n_ps_grid=sum(nzix(:));
        ph=zeros(n_ps_grid,n_ifg,'single');
    end
        ph(:,i1)=ph_grid(nzix);
end

n_ps = n_ps_grid;

fprintf('   Number of resampled points: %d\n',n_ps)

[nz_i,nz_j]=find(ph_grid~=0);
if pix_size==0
    xy = xy;
else
    xy=[grid_x_min+(nz_j-0.5)*pix_size,grid_y_min+(nz_i-0.5)*pix_size];
end
gridix=zeros(size(nzix));
gridix(nzix)=[1:n_ps];
% save('phwgrid.mat','ph','ph_lowpass','xy','ij','nzix','grid_x_min',...
%     'grid_y_min','n_i','n_j','n_ifg','n_ps','grid_ij','pix_size','-v7.3');
end
