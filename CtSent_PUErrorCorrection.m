function Ph1 = CtSent_PUErrorCorrection(n_image,nconsec,lonlat,Ph,flag)
% Phase Unwrapping Error Correction Based on Compressed Sensing
% Use Small-Baseline (Sequential Network) InSAR Data
% Created by Zhang-Feng Ma on 30/09/2020
% School of Earth Sciences and Engineering,Hohai University
% Email: jspcmazhangfeng@hhu.edu.cn

%CtSent_PUErrorCorrection Summary of this function goes here
%Phase Unwrapping Error Correction Based on Integer Linear Programming
%n_image: number of SLC
%nconsec:sequence of the sequential network
%lonlat: longitude and latitude matrix, Mx2
%ph: unwrapped phase matrix, MxN
%flag:1 (ILP), 2 (LS) , 1.5 (LASSO)

fprintf('#######################CtSent v1.1####################### \n');
fprintf('######################################################### \n');
fprintf('#################  PU Error Correction  ################# \n');
fprintf('######################################################### \n');
fprintf('########     Zhang-Feng Ma, Hohai University      ####### \n');
fprintf('###################  30,July,2020    #################### \n');

%replace NaN to zero
Ph(isnan(Ph)) = 0;
[npoints,nedges] = size(Ph);
fprintf('Calculating Residues & PU Errors... \n');
[~,~,edgs,eleidx] = TriangleNetwork(n_image,nconsec);
incMat = designIncMat(eleidx,nedges);
Ph1 = Ph;
% ambiguityNow = abs(Ph *(incMat'));
% Ph = setref(incMat,Ph,lonlat,ambiguityNow);
Ph = setref_auto(Ph,lonlat);
[~,~,residual2pi,~,resiual2piele] = ClosingLoops(eleidx,edgs,Ph,n_image);
fprintf('%d Triangle Loops in Total... \n',size(eleidx,1));
for i = 1:size(eleidx,1)
    fprintf('%d th Triangle Loop Contains %d 2-pi ambiguities \n',i,resiual2piele(i));
end
% incMat = designIncMat(eleidx,nedges);
idx = (1:npoints)';
idx(residual2pi==0) = [];
npointPUc = length(idx);
correctPh = zeros(npointPUc,nedges);
fprintf(' %d points contain PU errors... \n',npointPUc);
% parfor_progress(npointPUc);
fprintf('Starting Correcting PU errors... \n');
all_step = floor(npointPUc/100); p = 1;
for i = 1:npointPUc
    if flag == 1
    Modula2PI = LPsolver(Ph(idx(i),:)',incMat);
    elseif flag==1.5
    Modula2PI = ConstrucCLasso(eleidx,Ph(idx(i),:)');
    else
    Modula2PI = ConstrucC(eleidx,Ph(idx(i),:)');
    end
%     correctPh(i,:) = Ph(idx(i),:) + 2*pi*Modula2PI;
    correctPh(i,:) = 2*pi*Modula2PI;
    if i == all_step * p
       disp(['progress: ', num2str(1*p),'%']);
       p = p+1;
    end
%    ppm.increment();
% parfor_progress;
end
% parfor_progress(0);
Ph1(idx,:) = Ph1(idx,:) + correctPh;
end
function Modula2PI = LPsolver(Ph,incMat)
n_edge = size(Ph,1);
Aeq = [incMat,-incMat];
Beq = round(incMat*Ph./(-2*pi));
L = zeros(n_edge*2,1);%lower boundary
U = Inf(n_edge*2,1);%upper boundary
% fprintf('Starting to slove L1-Norm by integer linear programming(ILP)...\n');
options = optimoptions('intlinprog','Display','off');
% sol = linprog(double(ones(n_edge*2,1)),[],[],Aeq,double(Beq),L,U,options);
[sol,~] = intlinprog(double(ones(n_edge*2,1)),(1:size(L))',[],[],Aeq,double(Beq),L,U,options);
% fprintf('Ambiguities fixed! \n');
Hp = sol(1:n_edge,1);Hm = sol(n_edge+1:end,1);%plus charge and minus charge
Modula2PI = round(Hp-Hm)';
end
function incMat = designIncMat(eleidx,nedges)
ntri = size(eleidx,1);
fprintf('Constrcuting %d triangle loops... \n',ntri);
incMat = zeros(ntri,nedges);
fprintf('Designing incidence matrix... \n');
for i = 1:size(eleidx,1)
    incMat(i,eleidx(i,1)) = 1;incMat(i,eleidx(i,2)) = 1;incMat(i,eleidx(i,3)) = -1;
end
end
function Modula2PI = ConstrucC(eleidx,Ph)
n_edge = size(Ph,1);
ntri = size(eleidx,1);
incMat = zeros(size(eleidx,1),size(Ph,1));
for i = 1:size(eleidx,1)
    incMat(i,eleidx(i,1)) = 1;incMat(i,eleidx(i,2)) = 1;incMat(i,eleidx(i,3)) = -1;
end
Cmnk = incMat*Ph - angle(exp(1j*(incMat*Ph)));
goodedge = eleidx(find(abs(Cmnk)<0.1),:);
statistic = tabulate(goodedge(:));
count = statistic(:,2);
edgidx = statistic(:,1);
goodedgs = edgidx(count>=1);
Aadd = zeros(length(goodedgs),n_edge);
Aadd(sub2ind(size(Aadd),(1:length(goodedgs))',goodedgs)) = ones(size(goodedgs));
Aeq = [incMat;Aadd]*(-2*pi);
Beq = [incMat*Ph;zeros(size(goodedgs))];
U = pinv(Aeq'*Aeq)*(Aeq'*Beq);
Modula2PI = round(U');
end
function Modula2PI = ConstrucCLasso(eleidx,Ph)
n_edge = size(Ph,1);
ntri = size(eleidx,1);
incMat = zeros(size(eleidx,1),size(Ph,1));
for i = 1:size(eleidx,1)
    incMat(i,eleidx(i,1)) = 1;incMat(i,eleidx(i,2)) = 1;incMat(i,eleidx(i,3)) = -1;
end
Cmnk = incMat*Ph - angle(exp(1j*(incMat*Ph)));
Aeq = incMat*(-2*pi);
Beq = sum(Cmnk,2);
U = lasso(Aeq,Beq,'Lambda',0.1,'MaxIter',50);
Modula2PI = U';
end
function [looprob,tempCoh,resiual2pi,ph_sm,resiual2piele] = ClosingLoops(eleidx,edgs,Ph,n_image)
Ph(isnan(Ph)) = 0;
ntri = size(eleidx,1);
incMat = zeros(size(eleidx,1),size(Ph,2));
num=1;
p=1;
all = ntri;
all_step = floor(all/100);
for i = 1:ntri
    incMat(i,eleidx(i,1)) = 1;incMat(i,eleidx(i,2)) = 1;incMat(i,eleidx(i,3)) = -1;
    if num == all_step * p
       disp(['progress: ', num2str(100*p),'%']);
       p = p+1;
    end
end
residual = (incMat * Ph')';
resiual2pi = sum(abs(round(residual./(2*pi))),2);
resiual2piele = sum(abs(round(residual./(2*pi))),1);
% tempCoh = abs(sum(exp(1j*(residual)),2))./ntri;
Cmnk = residual - angle(exp(1j*(residual)));
count = sum(abs(Cmnk)>0.1,2);
looprob = count/size(eleidx,1)*100;
[ph_sm,ph_res] = uw_invert_sm(edgs,Ph,n_image);
tempCoh = abs(sum(exp(1j*(ph_res)),2))./size(ph_res,2);
end
function PhU = setref(A,PhU,lonlat,ambiguityNow)
inp = 0;
while inp == 0
figure;scatter(lonlat(:,1),lonlat(:,2),5,sum(ambiguityNow,2),'filled');box on;
grid on;set(gca,'tickdir','out');
set(gca,'gridlinestyle','--','GridColor',[1 1 1],'linewidth',1);
set(gca,'color',[0.6235    0.7137    0.8039]);axis tight;
cptcmap('GMT_haxby.cpt');colorbar('Location','southoutside');cbarrow;
%colorbar('Location','southoutside');
title('Setting reference...');
polycircle = drawcircle;
poly = polycircle.Vertices;
close;
ind = inpolygon(lonlat(:,1),lonlat(:,2),poly(:,1),poly(:,2));
% intcycles = round((PhU - angle(exp(1j.*PhU)))./(2*pi));
intcycles = (PhU-repmat(nanmean(PhU(ind,:),1),size(lonlat,1),1))*(A');
figure;subplot(2,1,1);imagesc(intcycles(ind,:));cptcmap('GMT_haxby.cpt');
subplot(2,1,2);plot(std(intcycles(ind,:)),'ko');
prompt = 'Are there ambiguities dominating the selected region ? [Yes:1 or No:0] ';
inp = input(prompt);
end
PhU = (PhU-repmat(nanmean(PhU(ind,:),1),size(lonlat,1),1));
end
function PhU = setref_auto(PhU,lonlat)
ref_lonlat = DBScanForReliable(lonlat(1:20:end,:),PhU(1:20:end,:));
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
function  [order,incMat,edgs,eleidx] = TriangleNetwork(stacksize,N)
%N: degree
%stack size
m=1;k=1;
for i=1:stacksize-1
    for j=i+1:i+N
        if j>=stacksize
            continue;
        end
   order(m,1)=i;order(m,2)=i+1; order(m,3)=j+1;
   incMat(m,i) = 1;incMat(m,i+1) = 1;incMat(m,j+1) = -1;
   edgs(k,1) = i;edgs(k,2) = i+1;edgs(k+1,1) = i+1;edgs(k+1,2) = j+1;edgs(k+2,1) = i;edgs(k+2,2) = j+1;
   k = k+3;
   m=m+1;
    end
end
edgs = unique(edgs,'rows');
for i = 1:size(order,1)
       temp = ismember(edgs,[order(i,1) order(i,2)],'rows');
    eleidx(i,1) = (find(temp == 1));
        temp = ismember(edgs,[order(i,2) order(i,3)],'rows');
    eleidx(i,2) = (find(temp == 1));
        temp = ismember(edgs,[order(i,1) order(i,3)],'rows');
    eleidx(i,3) = (find(temp == 1));
end
end
function ref_lonlat = DBScanForReliable(lonlat,ph)
tempcoh = edgesetTempcoh(ph,lonlat(:,1),lonlat(:,2));
[grid_lon,grid_lat,grid_value,~]=DownSamp2Grid(tempcoh,ph,lonlat,500);
[Data_Wh, ~, ~, ~] = whiten([grid_value,grid_lon,grid_lat],0.01);
label  =  dbscan(Data_Wh,0.2,200);
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
subplot(1,2,1);scatter(grid_lon,grid_lat,5,grid_value,'s','filled');colormap(jet);set(gca,'clim',[0.8 1]);
hold on;scatter(ref_lonlat(1),ref_lonlat(2),100,'o','filled');plot(lon(k),lat(k));
subplot(1,2,2);scatter(lon,lat,5,value,'s','filled');colormap(jet);set(gca,'clim',[0.8 1]);
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
function [label, mu, energy] = kmeans1(X, m)
% Input:
%   X: d x n data matrix
%   m: initialization parameter
% Output:
%   label: 1 x n sample labels
%   mu: d x k center of clusters
label = init(X, m);
n = numel(label);
idx = 1:n;
last = zeros(1,n);
while any(label ~= last)
    [~,~,last(:)] = unique(label);                  % remove empty clusters
    mu = X*normalize1(sparse(idx,last,1),1);         % compute cluster centers 
    [val,label] = min(dot(mu,mu,1)'/2-mu'*X,[],1);  % assign sample labels
end
energy = dot(X(:),X(:),1)+2*sum(val);
end
function label = init(X, m)
[d,n] = size(X);
if numel(m) == 1                           % random initialization
    mu = X(:,randperm(n,m));
    [~,label] = min(dot(mu,mu,1)'/2-mu'*X,[],1); 
elseif all(size(m) == [1,n])               % init with labels
    label = m;
elseif size(m,1) == d                      % init with seeds (centers)
    [~,label] = min(dot(m,m,1)'/2-m'*X,[],1); 
end
end
function Y = normalize1(X, dim)
if nargin == 1
    % Determine which dimension sum will use
    dim = find(size(X)~=1,1);
    if isempty(dim), dim = 1; end
end
Y = X./sum(X,dim);
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
function [ph_uw,ph_res] = uw_invert_sm(ifgday_ix,phU,n_image)
n_ifg = size(phU,2);
G=zeros(n_ifg,n_image);
for i=1:n_ifg
    G(i,ifgday_ix(i,1))=-1;
    G(i,ifgday_ix(i,2))=1;
end
G(:,1)=[]; % take out master as ref by setting to zero
ph_uw = lscov(G,double(phU'))';
ph_res=single(phU - (G*ph_uw')');
end
