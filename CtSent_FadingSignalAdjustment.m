function [fading_signal] = CtSent_FadingSignalAdjustment(lonlat,edgs,closure,max_rm_fraction)
n_edge = size(edgs,1);
n_points = size(lonlat,1);
n_tri = size(closure,2);
if nargin < 4
    max_rm_fraction = 0.005; % maximum fraction of outlier edges to remove at one time
%     smooth_threshold = 0.5;
end
dph_space_uw = double(angle(exp(1j*(closure(edgs(:,1),:) - closure(edgs(:,2),:)))));
res_tol=1e-2;         % tolerance for dph residual
ref_ix = find(sum(abs(angle(exp(1j*closure))),2)==min(sum(abs(angle(exp(1j*closure))),2)));
ref_ix=ref_ix(1);
idx = [(1:n_edge)';(1:n_edge)';];
idy = [edgs(:,1);edgs(:,2);];
value = [ones(size(edgs(:,1)));-ones(size(edgs(:,2)));];
A = sparse(idx,idy,value);
A=A(:,[1:ref_ix-1,ref_ix+1:n_points]);
A = double(A);
A_unweighted=A;
fading_signal=zeros(n_points,n_tri);
progressBar = CommandLineProgressBar(n_tri);
progressBar.message = 'Estimating Fading Signal: ';
progressBar.barLength = 42;
parfor i = 1:n_tri
fading = inversionFading(edgs,A_unweighted,dph_space_uw(:,i),max_rm_fraction,res_tol,ref_ix,n_edge,n_points);
intcycle = round(mean(fading - closure(:,i))./(2*pi))*2*pi;
fading = fading - intcycle;
fading_signal(:,i) = fading;
progressBar.increment;
end
end
function fading_signal = inversionFading(edgs,A_unweighted,dph_space_uw,max_rm_fraction,res_tol,ref_ix,n_edge,n_points)
warning off;
    dph_abs_exp = exp(abs(dph_space_uw));
    weighting_use = 1./dph_abs_exp;
    A = spdiags(weighting_use,0, n_edge, n_edge)*A_unweighted;
    dph_use = dph_space_uw.*weighting_use;
    exit_flag = 0;
    edges_use = edgs;
    n_dph = length(dph_use);
    n_bad = ceil(length(dph_use)*max_rm_fraction);
    while exit_flag == 0
        if sprank(A) >= size(A,2)
            fading_temp = A\dph_use;
            A_save = A;
            edges_save = edges_use;
            weighting_save = weighting_use;
            dph_save = dph_use; 
            n_dph_save = n_dph;
            dph_hat = A_save*fading_temp;
            res = (dph_save-dph_hat);            
        else
            fading_temp=A_save\dph_save;
            n_bad=ceil(n_bad/10);
        end
        if abs(res) < res_tol
            exit_flag=1;
        else
            [~,sort_ix] = sort(abs(res));
            arc_in = true(n_dph_save,1);   
            arc_out = zeros(n_points,1);
            for i2=length(sort_ix):-1:length(sort_ix)-n_bad+1
                    bad_ix=(sort_ix(i2));
                    if arc_out(edges_save(bad_ix,1:2))==0
                        arc_in(bad_ix)=0;
                        arc_out(edges_save(bad_ix,1:2))=1;
                    end
            end
            edges_use=edges_save(arc_in,:);
            weighting_use=weighting_save(arc_in);
            dph_use=dph_save(arc_in,:);
            A=A_save(arc_in,:);
            n_dph=length(dph_use);
        end
    end
    fading_signal = [fading_temp(1:ref_ix-1);0;fading_temp(ref_ix:end)]; 
end
function smooth = ctSent_PhaseSmooth(ph,X,Y,K,complexflag)
    Indice = knnsearch([X,Y],[X,Y],'K',K);
    tempdiff = zeros(length(X),1);
    for i = 1:K-1
        edgeIndice(:,1) = Indice(:,1);
        edgeIndice(:,2) = Indice(:,i+1);
        temp = temporaldiff(ph,edgeIndice,complexflag);
        tempdiff = tempdiff + temp;
    end
  smooth = sqrt(tempdiff/(K-1)); 
end
function temp = temporaldiff(ph,edgs,complexflag)
if complexflag
temp  = abs(angle(exp(1j*(ph(edgs(:,2)) - ph(edgs(:,1)))))).^2;
else
temp  = abs((ph(edgs(:,2)) - ph(edgs(:,1)))).^2;
end
end
