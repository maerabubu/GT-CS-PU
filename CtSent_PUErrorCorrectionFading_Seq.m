function  [PhU,fading] = CtSent_PUErrorCorrectionFading_Seq(ac_date,lonlat,ph,order,methodflag,fadingflag,networktype)
% Phase Unwrapping Error Correction and Fading Signal Correction Based on
% Basis Pursuit Denoising and Closure-based Stacking
% Use Small-Baseline (Sequential Network) InSAR Data
% Created by Zhang-Feng Ma on 06/12/2022
% Earth Observatory of Singapore
% Email: jspcmazhangfeng@gmail.com

%step 1. automatic reference point selection
[ph,robust_ind] = setref_auto(ph,lonlat);
%step 2. point-wise unwrapping error correction & fading signal estimation
[PhU,fading] = CtSent_PUErrorCorrectionPoint(ac_date,order,ph,lonlat,robust_ind,methodflag,fadingflag,networktype,1);
end
function [PhU,fading] = CtSent_PUErrorCorrectionPoint(ac_date,order,Ph,lonlat,robust_ind,methodflag,fadingflag,networktype,interval)
% Phase Unwrapping Error Correction and Fading Signal Correction Based on
% Basis Pursuit Denoising and Closure-based Stacking
% Use Small-Baseline (Sequential Network) InSAR Data
% Created by Zhang-Feng Ma on 06/12/2022
% Earth Observatory of Singapore
% Email: jspcmazhangfeng@gmail.com

%CtSent_PUErrorCorrectionFading Summary of this function goes here
%Phase Unwrapping Error Correction Based on BPDN
%n_image: number of SLC
%nconsec:sequence of the sequential network
%lonlat: longitude and latitude matrix, Mx2
%ph: unwrapped phase matrix, MxN

fprintf('#######################CtSent v2.1####################### \n');
fprintf('######################################################### \n');
fprintf('#############  PU Error Correction (Point) ############## \n');
fprintf('######################################################### \n');
fprintf('######Zhangfeng Ma, Earth Observatory of Singapore####### \n');
fprintf('###################   21,Nov,2022    #################### \n');

%replace NaN to zero
Ph(isnan(Ph)) = 0;
n_image = length(ac_date);
fprintf('Calculating Residuals ... \n');
% [~,incMat,order,eleidx] = TriangleNetwork(ac_date,n_image,nconsec);
[incMat,eleidx] = findtri(order);
[eleidx_full,incMat1,order_full,~] = TriangleNetwork(ac_date,n_image,n_image);
[ia,ib] = ismember(order_full,order,'rows');
incMat1(:,find(ib==0)) = [];
incMat1(sum(abs(incMat1),2)<3,:) = [];
unwerr = round((Ph*incMat1')./(2*pi));
noerrorNum_before = sum(sum(abs(unwerr),2)==0);
clear unwerr;
[addCharg,fading] = SovePointAmb(ac_date,order,order_full,Ph,eleidx,incMat,incMat1,lonlat,robust_ind,methodflag,fadingflag,networktype,interval);
PhU = Ph + addCharg * 2*pi;
unwerr = round((PhU*incMat')./(2*pi));
noerrorNum_after = sum(sum(abs(unwerr),2)==0);
fprintf('Phase Unwrapping Error Correction Done! \n');
fprintf('%d Points with Unw Errors are corrected \n',noerrorNum_after - noerrorNum_before);
fprintf('%d Points with Unw Errors still cannot be corrected \n',size(PhU,1) - noerrorNum_after);
end
function [intAmbAll,fading] = SovePointAmb(ac_date,order,order_full,ph,eleidx,incMat,incMat1,lonlat,robust_ind,methodflag,fadingflag,networktype,interval)
gcycle = round(incMat*ph'./(-2*pi))';
zeronumber = sum(gcycle==0,2);
correct_ind = zeronumber==size(gcycle,2);
idx = find(correct_ind==0);
% n_corr = sum(correct_ind==0);
intAmbAll = zeros(size(ph));
incMat = zeros(size(eleidx,1),size(ph,2));
for i = 1:size(eleidx,1)
    incMat(i,eleidx(i,1)) = 1;incMat(i,eleidx(i,2)) = 1;incMat(i,eleidx(i,3)) = -1;
end
closure = (incMat1*ph')';
closure0 = (incMat*ph')';
if fadingflag
fading = ctSent_fading_stacking(ac_date,order,eleidx,lonlat,closure0,networktype);
fading_loop = fading*incMat1';
closure = closure - fading_loop;
else
    fading = 0;
end
nifg = size(ph,2);
closure = closure(correct_ind==0,:);
closure = closure(1:interval:end,:);
idx = idx(1:interval:end);
n_corr = length(idx);
lambda = sqrt(sum(abs(angle(exp(1j*(closure)))).^2,2)./size(closure,2));
intAmb = zeros(n_corr,size(ph,2));
progressBar = CommandLineProgressBar(n_corr);
progressBar.message = 'Solving Point Integer Ambiguity: ';
progressBar.barLength = 42;
parfor count = 1:n_corr
weight = weightgoodifg(closure(count,:)',eleidx,nifg);
if strcmp(methodflag,'ilp')
temp = Lpsolver(incMat1,closure(count,:)',weight);
intAmb(count,:) = temp;
elseif strcmp(methodflag,'bpdn')
temp = BPDNsolver(double(incMat1),double(closure(count,:)'./(-2*pi)),double(lambda(count)./(2*pi)));
temp = round(temp);
intAmb(count,:) = temp;
elseif strcmp(methodflag,'lasso')
temp = LASSOsolver(incMat1,closure(count,:)',0.01);    
intAmb(count,:) = temp;
else
    warning('Please choose intlinprog or bpdn method !!!\n');
end
progressBar.increment;
end
intAmbAll(idx,:) =  intAmb;
end
function [A,eleidx] = findtri(Intflist)
%findtri: find all triangle loops in SBAS graph
%Intflist: two-columns ifg id matrix [reference secondary]
Intf = unique(Intflist(:));
Intf = sort(Intf,'ascend');
alltri = nchoosek(Intf,3);
ntri = size(alltri,1);
eleidx = zeros(ntri,3);
progressBar = CommandLineProgressBar(ntri);
progressBar.message = 'Search All Possible Triangle Loops: ';
progressBar.barLength = 42;
parfor i = 1:ntri
    triedgs = [alltri(i,1),alltri(i,2);alltri(i,2),alltri(i,3);alltri(i,1),alltri(i,3)];
    [temp,ix] = ismember(triedgs,Intflist,'rows');
    if sum(temp)==3
        eleidx(i,:) = ix;
    end  
    progressBar.increment;
end
eleidx(sum(eleidx,2)==0,:) = [];
ntri = size(eleidx,1);
progressBar = CommandLineProgressBar(ntri);
progressBar.message = 'Generating Incidence Matrix: ';
progressBar.barLength = 42;
parfor i = 1:ntri
   tempA = zeros(1,size(Intflist,1));
   tempA(eleidx(i,1)) = 1;tempA(eleidx(i,2)) = 1;tempA(eleidx(i,3)) = -1;
   A(i,:) = tempA;
   progressBar.increment;
end
end
function weight = weightgoodifg(closure,eleidx,nifg)
goodedge = eleidx(find(round(closure./(2*pi))==0),:);
statistic = tabulate(goodedge(:));
count = [(1:nifg)',zeros(nifg,1)];
count(statistic(:,1),2) = statistic(:,2);
weight = ones(nifg,1)*1000;
weight(count(:,2)==0) = 0.1;
end
function [PhU,robust_ind] = setref_auto(PhU,lonlat)
[ref_lonlat,robust_ind] = DBScanForReliable(lonlat(1:5:end,:),PhU(1:5:end,:));
t = linspace(0,2*pi,114); % 114 pts
r=200;                   % change radius if you want to include more points
x = r*cos(t);
y = r*sin(t);
%
% SELECT POINTS based on LONLAT
xy=llh2local(lonlat',[ref_lonlat(1),ref_lonlat(2)])'*1000;
ind = inpolygon(xy(:,1),xy(:,2),x,y);
PhU = (PhU-repmat(nanmean(PhU(ind,:),1),size(lonlat,1),1));
end
function [ref_lonlat,bound] = DBScanForReliable(lonlat,ph)
tempcoh = edgesetTempcoh(ph,lonlat(:,1),lonlat(:,2));
[grid_lon,grid_lat,grid_value,~]=DownSamp2Grid(tempcoh,ph,lonlat,1000);
[Data_Wh, ~, ~, ~] = whiten([grid_value,grid_lon,grid_lat],0.01);
label  =  dbscan(Data_Wh,0.2,100);
labeluni = unique(label);
for i = 1:length(labeluni)
    meancoh(i) = mean(grid_value(label == labeluni(i)));
end
idx = labeluni(meancoh == max(meancoh));
ind = label==idx;
lon = grid_lon(ind);
lat = grid_lat(ind);
value = grid_value(ind);
k = boundary(double(lon),double(lat));
ref_lonlat = [lon(value == max(value)),lat(value == max(value))];
bound = [lon(k),lat(k)];
% subplot(1,2,1);scatter(grid_lon,grid_lat,5,grid_value,'s','filled');colormap(jet);set(gca,'clim',[0.8 1]);
% hold on;scatter(ref_lonlat(1),ref_lonlat(2),100,'o','filled');plot(lon(k),lat(k));
% subplot(1,2,2);scatter(lon,lat,5,value,'s','filled');colormap(jet);set(gca,'clim',[0.8 1]);
end
function tempcoh = edgesetTempcoh(ph,X,Y)
    Indice = knnsearch([X,Y],[X,Y],'K',21);
    tempcoh = zeros(length(X),1);
    for i = 1:20
        edgeIndice(:,1) = Indice(:,1);
        edgeIndice(:,2) = Indice(:,i+1);
        temp = temporalcoherence(ph,edgeIndice);
        tempcoh = tempcoh + temp;
    end
  tempcoh = tempcoh/20; 
end
function coh = temporalcoherence(ph,edgs)
nimage = size(ph,2);
temp = zeros(size(edgs,1),1);
ph = exp(1j*ph);
for i = 1:nimage
    temp = temp + (ph(edgs(:,2),i).*conj(ph(edgs(:,1),i)));
end
coh = abs(temp)/nimage;
end
function [grid_lon,grid_lat,grid_value,grid_value1] = DownSamp2Grid(tpc,ph,lonlat,pix_size)
%this function is used to regularly downsample the InSAR time series
if nargin < 3
    pix_size = 500;
end
%read data from path 1
min_lon_common = min(lonlat(:,1));
max_lon_common = max(lonlat(:,1));
min_lat_common = min(lonlat(:,2));
max_lat_common = max(lonlat(:,2));

pix_size = km2deg(pix_size/1000);

y = min_lat_common:pix_size:max_lat_common;
x = min_lon_common:pix_size:max_lon_common;
grid_value = zeros(length(y),length(x),'single');
grid_value1 = zeros(length(y),length(x),'single');
numInGrid1 = grid_value;
% grid_la1 = grid_value;
grid_lon = grid_value;
grid_lat = grid_value;
for ii=1:length(y)-1
    for jj=1:length(x)-1
        ix1 = lonlat(:,1)>=x(jj)&lonlat(:,1)<=x(jj+1)&lonlat(:,2)>=y(ii)&lonlat(:,2)<=y(ii+1);
        if sum(ix1)~=0
            grid_value(ii,jj) = mean(tpc(ix1));
            grid_value1(ii,jj) = mean(std(ph(ix1,:)));
%             grid_la1(ii,jj) = mean(inc(ix1));
            numInGrid1(ii,jj) = sum(ix1);
        end
        grid_lon(ii,jj) = (x(jj)+x(jj+1))*0.5;
        grid_lat(ii,jj) = (y(ii)+y(ii+1))*0.5;
    end
end
grid_value=grid_value(:);
grid_value1=grid_value1(:);
ind = grid_value==0|isnan(grid_value);
grid_lon=grid_lon(:);
grid_lat=grid_lat(:);
grid_value(ind)=[];
grid_value1(ind)=[];
grid_lon(ind)=[];
grid_lat(ind)=[];
end
function [Xwh, mu, invMat, whMat] = whiten(X,epsilon)


if ~exist('epsilon','var')
    epsilon = 0.0001;
end

mu = mean(X); 
X = bsxfun(@minus, X, mu);
A = X'*X;
[V,D,notused] = svd(A);
whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
Xwh = X*whMat;  
invMat = pinv(whMat);

end 
function fading_out = FadingNeigh(closure,X,Y,networktype)
if strcmp(networktype,'delaunay')
tr = delaunayTriangulation(double(X),double(Y));
edg = edges(tr);
elseif strcmp(networktype,'apsp')
edg = APSPnet(closure,double(X),double(Y));
else
    warning('Please define the netwrok type [delaunay or apsp]');
end
delete(gcp('nocreate'));
% parpool('local',8);
% fading = CtSent_NetworkAdjustment([X Y],edg,closure);
fading_out = CtSent_FadingSignalAdjustment([X Y],edg,closure);
% delete(gcp('nocreate'));
end
function Modula2PI = Lpsolver(incMat,Ph,weight)
n_edge = size(incMat,2);
incMatrix = [incMat,-incMat];
Aeq = incMatrix;
Beq = round(Ph./(-2*pi));
weight = [weight;weight];
L = zeros(n_edge*2,1);%lower boundary
U = Inf(n_edge*2,1);%upper boundary
options = optimoptions('intlinprog','Display','off');
[sol,~] = intlinprog(weight,(1:size(L))',[],[],Aeq,double(Beq),L,U,options);
Hp = sol(1:n_edge,1);Hm = sol(n_edge+1:end,1);%plus charge and minus charge
Modula2PI = (Hp-Hm);
end
function Modula2PI = LASSOsolver(incMat,Ph,lambda)
Aeq = incMat;
Beq = round(Ph/(-2*pi));
warning off;
U = lasso(Aeq,Beq,'Lambda',lambda);
Modula2PI = round(U);
end
function edgs = APSPnet(ph,X,Y)
 lengthPh = size(ph,1);lengthimg = size(ph,2);
    K = 21;edgeIndice = zeros(lengthPh*(K-1),2,'single');
%     fprintf("Generating Delaunay Triangulation Network...\n");
    tr = delaunayTriangulation(double(X(:)),double(Y(:)));
    Indice = knnsearch([X,Y],[X,Y],'K',K);
    for i = 1:K-1
        edgeIndice(lengthPh*(i-1)+1:lengthPh*i,1) = Indice(:,1);
        edgeIndice(lengthPh*(i-1)+1:lengthPh*i,2) = Indice(:,i+1);
    end
    fulledgs = sort(edgeIndice',1);
    fulledgs = fulledgs';
    tredgs = edges(tr);
    edgs = unique([fulledgs;tredgs],'rows');
%     phc = exp(1j*ph);
   temp_res = temporalresidual(ph,edgs);
% G=graph(edgs(:,1),edgs(:,2),temp_coh,259361);
% p = plot(G,'XData',X,'YData',Y,'EdgeCData',temp_coh,'EdgeColor','flat','Marker','none');colormap jet
    edgs = APSPTri(tr.Points,tredgs,edgs,temp_res.^2);
end
function ft = APSPTri(Points,iniedgs,xyfull,dismat)
%iniedgs: initial edges
%dismat : distance matrix of full edges , size: n*(n-1)/2
% iniedgs = double(iniedgs);xyfull = double(xyfull);dismat = double(dismat);Points = double(Points);
G=graph(xyfull(:,1),xyfull(:,2),dismat);
% G = sparse(double(xyfull(:,1)),double(xyfull(:,2)),double(dismat),size(Points,1),size(Points,1));
max_edge_num = 2000;
path1 = zeros(size(iniedgs,1),max_edge_num,'int32');
% parfor_progress(size(iniedgs,1));
fprintf('Constructing APSP Network... \n');
progressBar = CommandLineProgressBar(size(iniedgs,1));
progressBar.message = 'Searching: ';
progressBar.barLength = 42;
    parfor i=1:size(iniedgs,1)
        pathtemp = shortestpath(G,iniedgs(i,1),iniedgs(i,2));
        nedge = length(pathtemp);
        tempedge = zeros(max_edge_num,1);
        tempedge(1:nedge) = int32(pathtemp);
        path1(i,:) = tempedge;
        progressBar.increment;
    end
    startpath = path1(:,1:end-1);
    endpath = path1(:,2:end);
    clear path1;
    startpath = startpath(:);
    endpath = endpath(:);
    id = startpath==0 | endpath==0;
    startpath(id) = [];endpath(id) = [];
    qfrom    = min(startpath,endpath);% index in X_REF
    qto      = max(startpath,endpath);% index in X_REF
    ft       = [qfrom, qto];
    ft       = double(unique(ft,'rows'));
end
function err = temporalresidual(ph,edgs)
nimage = size(ph,2);
temp = zeros(size(edgs,1),1);
for i = 1:nimage
    temp = temp + (angle(exp(1j*(ph(edgs(:,2),i) - ph(edgs(:,1),i))))).^2;
end
err = sqrt(temp/nimage);
end
function  [order,A,edgs,eleidx] = TriangleNetwork(ac_date,stacksize,N)
%N: degree
%stack size
day = datenum(num2str(ac_date),'yyyymmdd');
m=1;k=1;
for i=1:stacksize-1
    for j=i+1:i+N
        if j>=stacksize
            continue;
        end
   order(m,1)=i;order(m,2)=i+1; order(m,3)=j+1;
   edgs(k,1) = i;edgs(k,2) = i+1;edgs(k+1,1) = i+1;edgs(k+1,2) = j+1;edgs(k+2,1) = i;edgs(k+2,2) = j+1;
   k = k+3;
   m=m+1;
    end
end
edgs = unique(edgs,'rows');
t_day = day(edgs(:,2)) - day(edgs(:,1));
for i = 1:size(order,1)
       temp = ismember(edgs,[order(i,1) order(i,2)],'rows');
    eleidx(i,1) = (find(temp == 1));
        temp = ismember(edgs,[order(i,2) order(i,3)],'rows');
    eleidx(i,2) = (find(temp == 1));
        temp = ismember(edgs,[order(i,1) order(i,3)],'rows');
    eleidx(i,3) = (find(temp == 1));
end
% day_tri = t_day(eleidx);
% [sorted_ele,~] = sortrows([day_tri,eleidx],[1 2 3]);
% eleidx = sorted_ele(:,4:6);
ntri = size(eleidx,1);
A = zeros(ntri,size(edgs,1));
progressBar = CommandLineProgressBar(ntri);
progressBar.message = 'Generating Incidence Matrix: ';
progressBar.barLength = 42;
parfor i = 1:ntri
   tempA = zeros(1,size(edgs,1));
   tempA(eleidx(i,1)) = 1;tempA(eleidx(i,2)) = 1;tempA(eleidx(i,3)) = -1;
   A(i,:) = tempA;
   progressBar.increment;
end
end
function xy=llh2local(llh,origin)
%llh2local     xy=llh2local(llh,origin)
%Set ellipsoid constants (WGS84)

   a=6378137.0;
   e=0.08209443794970;

%Convert to radians
%   llh=llh*pi/180;
%   origin=origin*pi/180;
   llh=double(llh)*pi/180;
   origin=double(origin)*pi/180;

%Do the projection

   z=llh(2,:)~=0;

   dlambda=llh(1,z)-origin(1);

   M=a*((1-e^2/4-3*e^4/64-5*e^6/256)*llh(2,z) - ...
        (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*llh(2,z)) + ...
        (15*e^4/256 +45*e^6/1024)*sin(4*llh(2,z)) - ...
        (35*e^6/3072)*sin(6*llh(2,z)));

   M0=a*((1-e^2/4-3*e^4/64-5*e^6/256)*origin(2) - ...
        (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*origin(2)) + ...
        (15*e^4/256 +45*e^6/1024)*sin(4*origin(2)) - ...
        (35*e^6/3072)*sin(6*origin(2)));
   
   N=a./sqrt(1-e^2*sin(llh(2,z)).^2);
   E=dlambda.*sin(llh(2,z));

   xy(1,z)=N.*cot(llh(2,z)).*sin(E);
   xy(2,z)=M-M0+N.*cot(llh(2,z)).*(1-cos(E));

%Handle special case of latitude = 0

   dlambda=llh(1,~z)-origin(1);
   xy(1,~z)=a*dlambda;
   xy(2,~z)=-M0;

%Convert to km
   
   xy=xy/1000;
end
function x_out = BPDNsolver(A,y,lambda)
numCols = size(A,2);
lb = ones(numCols,1)*-10;ub = ones(numCols,1)*10;
cvx_begin('quiet')
    % cvx_precision('best');
    cvx_solver gurobi
    variable x(numCols) 
%     minimize(lambda*norm(x,1)+0.5*power(norm(A*x-y,2),2))
%     subject to
%         lb <= x <= ub
    minimize(norm(x,1))
    subject to
    norm(A*x-y,2) <= lambda
    lb <= x <= ub
cvx_end
% x= round(x);
x_out = x;
end
function fading_sig = ctSent_fading_stacking(ac_date,order,eleidx,lonlat,fading,networktype)
% [~,eleidx,~,~,fading] = findtri(ac_date,order,lonlat,ph,networktype);
nifg = size(order,1);
ntri = size(eleidx,1);
Inc = zeros(ntri,nifg);
fading = angle(exp(1j*fading));
for i = 1:nifg
    [temp,~] = ismember(eleidx(:,1),i,'rows');
    [temp1,~] = ismember(eleidx(:,2),i,'rows');
    [temp2,~] = ismember(eleidx(:,3),i,'rows');
    tempColumb = single(temp) + single(temp1) - single(temp2);
    Inc(:,i) = tempColumb;
end
normalizeN = sum(abs(Inc),1);
normalizeNmat = repmat(normalizeN,ntri,1);
Inc = Inc./normalizeNmat;
fading_out = FadingNeigh(fading,lonlat(:,1),lonlat(:,2),networktype);
fading_sig = fading_out * Inc;
end
