function phuw = CtSent_Interp4PU(ind,nzix,grid_ij,phw,phu)

fprintf('Unwrapping from grid...\n')
phwgrid = phw(~ind,:);
nfull = sum(ind);
phufull = phu(1:nfull,:);
phugrid = phu(nfull+1:end,:);
[n_ps,n_ifg]=size(phugrid);
n_psfull = size(phwgrid,1);
gridix=zeros(size(nzix));
gridix(nzix)=[1:n_ps];
    
ph_uwgrid=zeros(n_psfull,n_ifg,'single');
for i=1:n_psfull
    ix=gridix(grid_ij(i,1),grid_ij(i,2));
    if ix==0
        ph_uwgrid(i,:)=nan; % wrapped phase values were zero
    else
        ph_uw_pix=phugrid(ix,:);
        if isreal(phugrid)
          ph_uwgrid(i,:)=ph_uw_pix+angle(exp(1i*(phwgrid(i,:)-ph_uw_pix)));
        else
          ph_uwgrid(i,:)=ph_uw_pix+angle(phwgrid(i,:).*exp(-1i*ph_uw_pix));
        end
    end
end
phuw = zeros(size(phw));
phuw(ind,:) = phufull;
phuw(~ind,:) = ph_uwgrid;
end

