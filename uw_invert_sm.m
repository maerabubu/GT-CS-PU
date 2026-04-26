function [ph_uw,res,res1] = uw_invert_sm(ifgday_ix,phU,n_image)
n_ifg = size(phU,2);
G=zeros(n_ifg,n_image);
for i=1:n_ifg
    G(i,ifgday_ix(i,1))=1;
    G(i,ifgday_ix(i,2))=-1;
end
G(:,1)=[]; % take out master as ref by setting to zero
% ph_uw = lscov(G,double(phU'))';
ph_uw = (pinv(G'*G)*G'*phU')';
ph_res=single(phU - (G*ph_uw')');
ph_uw = [ph_uw(:,1)*0,ph_uw];
res = sqrt(sum(abs(ph_res).^2,2)./(n_ifg-1));
res1 = abs(mean(ph_res,2));
% ph_model = (G*ph_uw')';
end


