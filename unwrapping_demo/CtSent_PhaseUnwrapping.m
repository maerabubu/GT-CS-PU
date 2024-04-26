function [PhU,edgs,temp_coh,msd] = CtSent_PhaseUnwrapping(X,Y,ph,flag)
%CtSent_PhaseUnwrapping Summary of this function goes here
%Sparse Phase Unwrapping Based on All-Pair-Shortest-Path algorithm
%X:longitude or UTM local X, Mx1 vector
%Y:latitude or UTM local Y, Mx1 vector
%ph: phase matrix, MxN matrix
%flag: 1- APSP 2-Delaunay(original MCF)

fprintf('#######################CtSent v1.1####################### \n');
fprintf('######################################################### \n');
fprintf('#################    Phase Unwrapping   ################# \n');
fprintf('######################################################### \n');
fprintf('########     Zhang-Feng Ma, Hohai University      ####### \n');
fprintf('###################  30,July,2019    #################### \n');
addpath /home/zhangfeng.ma/gurobi1002/linux64/matlab
% load("pu_input0.mat");
% X = ll(:,1);Y = ll(:,2);flag = 1;ph = ph_ll;

    if flag == 1
    lengthPh = size(ph,1);lengthimg = size(ph,2);
    if exist("apspnetwork.mat","file")
        load("apspnetwork.mat","edgs","temp_coh");
    else
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
    temp_coh = temporalcoherence(ph,edgs);
    % G=graph(edgs(:,1),edgs(:,2),temp_coh,259361);
    % p = plot(G,'XData',X,'YData',Y,'EdgeCData',temp_coh,'EdgeColor','flat','Marker','none');colormap jet
    temp_coh0 = temporalcoherence(ph,tredgs);
    ind = temp_coh0 <0.985;% Risky Edges to be replaced    
    edgs0 = APSPTri(tr.Points,tredgs(ind,:),edgs,1./temp_coh.^20);
    edgs = [tredgs(~ind,:);edgs0];
    edgs = unique(edgs,'rows');
    temp_coh = temporalcoherence(ph,edgs);
    save('apspnetwork.mat','edgs','temp_coh','X','Y','-v7.3');
    end
%     diffph = phc(edgs(:,2),:).*conj(phc(edgs(:,1),:));
%     clear phc;
    weight = 1./temp_coh;
%     weight = 1./(abs(sum(diffph,2))./lengthimg);
    else
    lengthPh = size(ph,1);lengthimg = size(ph,2);
%     fprintf("Generating Delaunay Triangulation Network...\n");
    tr = delaunayTriangulation(double(X(:)),double(Y(:)));
    edgs = edges(tr);
    temp_coh = temporalcoherence(ph,edgs);
    weight = 1./temp_coh;
    end
    PhU = zeros(size(ph));
%     progressBar = CommandLineProgressBar(lengthimg);
%     progressBar.message = 'Mutiple Phase Unwrapping: ';
%     progressBar.barLength = 42;
    for i = 1:lengthimg
    tic;
    [ambig,msd(i)] = MinimumCostFlow(edgs,ph(:,i),weight);
    timecost = toc;
    PhU(:,i) = ph(:,i) + (ambig - ambig(1))*2*pi;
    temp = PhU(:,i);
%     save('PhU_temp.mat','temp');
%     progressBar.increment;
    fprintf('%d / %d ifg completed in %d s ...\n',i,lengthimg,round(timecost));
    end 
%     save('PhU_apsp.mat','PhU','msd','-v7.3');
end
function [Ambiguity,fval] = MinimumCostFlow(edg,Ph,weight)
npoints = size(Ph,1);
n_edge = size(edg,1);
edg = double([(1:n_edge)',edg]);
% weight = ones(n_edge,1);
idx = [(1:n_edge)';(1:n_edge)';(1:n_edge)';(1:n_edge)'];
idy = [edg(:,3);edg(:,2);((1:n_edge)+ npoints)';(((n_edge+1):(2*n_edge))+ npoints)'];
value = [ones(size(edg(:,2)));-ones(size(edg(:,3)));ones(size(edg(:,2)));-ones(size(edg(:,3)))];
Aeq=sparse(idx,idy,value);
%constructing Residues(charge -1,0,+1)
    charg = Residues(edg,Ph);
    weighting = [zeros(npoints,1);weight;weight];
    L = zeros(npoints+2*n_edge,1);%lower boundary
    U = ones(size(L))*100;%upper boundary
%     sol = linprog(double(weighting),[],[],Aeq,double(charg),L,U);
%     sol = round(sol);
options = optimoptions('intlinprog','Display','iter');
    [sol,fval] = intlinprog(double(weighting),(1:size(L))',[],[],Aeq,double(charg),L,U,[],options);
    Ambiguity = sol(1:npoints,1);%plus charge and minus charge
end
function chrg = Residues(edg,Ph)
%calculated residues
chrg = round((Ph(edg(:,2)) - Ph(edg(:,3))) / (2*pi));
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
function ft = APSPTri(Points,iniedgs,xyfull,dismat)
%iniedgs: initial edges
%dismat : distance matrix of full edges , size: n*(n-1)/2
% iniedgs = double(iniedgs);xyfull = double(xyfull);dismat = double(dismat);Points = double(Points);
G=graph(xyfull(:,1),xyfull(:,2),dismat);
% G = sparse(double(xyfull(:,1)),double(xyfull(:,2)),double(dismat),size(Points,1),size(Points,1));
    startpath=[];
    endpath=[];
%     all = size(iniedgs,1);
path1 = cell(size(iniedgs,1),1);
% parfor_progress(size(iniedgs,1));
fprintf('Starting Parallel Pool & Searching APSP... \n');
tic;
progressBar = CommandLineProgressBar(size(iniedgs,1));
progressBar.message = 'Generating APSP Network: ';
progressBar.barLength = 42;
    parfor i=1:size(iniedgs,1)
        [path1{i}] = shortestpath(G,iniedgs(i,1),iniedgs(i,2));
        progressBar.increment;
    end
progressBar = CommandLineProgressBar(size(iniedgs,1));
progressBar.message = 'Assemblying APSP Network: ';
progressBar.barLength = 42;
parfor i = 1:size(iniedgs,1)
        temp = path1{i};
        startpath = [startpath,temp(1:end-1)];
        endpath = [endpath,temp(2:end)];
progressBar.increment;
end
    qfrom    = min(startpath,endpath);% index in X_REF
    qto      = max(startpath,endpath);% index in X_REF
    ft       = [qfrom.', qto.'];
    ft       = unique(ft,'rows');
t = toc;
fprintf('Finish Searching APSP in %d s... \n',t);
%     triplot(Tri);
end
function [arcTriClosureIndex] = arcTriangleDetect(IDX_temp)
%IDX_temp   = [IDX_from,IDX_to];

arc_kp_num = size(IDX_temp,1);
IDX_point  = unique(IDX_temp(:));
IDX_temp1  = IDX_temp;
[~,IDX_temp1(:,1)] = ismember(IDX_temp(:,1),IDX_point);
[~,IDX_temp1(:,2)] = ismember(IDX_temp(:,2),IDX_point);

IDX_temp1  = sort(IDX_temp1,2);
IDX_index  = [IDX_temp1(:,1),IDX_temp1(:,2),(1:arc_kp_num)'];

frompNum=max(max(IDX_index(:,1:2)));
ADJ_matrix=sparse(IDX_index(:,1),IDX_index(:,2),1,frompNum,frompNum);
pointClosure_index_cell={};

for ii=1:frompNum
    if (mod(ii,100000)==0 && ii>1)
        disp(['LS solver:: currently at points ', num2str(ii), ' of ', num2str(frompNum)]);
    end
    toPoints=find(ADJ_matrix(ii,:));
    if length(toPoints)>=2
        arcSelect2=nchoosek(toPoints,2);
        tirdPoints=[];
        for jj=1:size(arcSelect2,1)
            tirdPoints(jj,1)=ADJ_matrix(arcSelect2(jj,1),arcSelect2(jj,2));
        end
        iy=find(tirdPoints);
        if ~isempty(iy)
            arcSelect3=[ones(length(iy),1)*ii,arcSelect2(iy,:)];
            pointClosure_index_cell{ii,1}=arcSelect3;
        end
    end
end


pointClosure_index_cell(cellfun(@isempty,pointClosure_index_cell))=[];
pointClosure=cell2mat(pointClosure_index_cell);
pointClosure=sort(pointClosure,2);
clear pointClosure_index_cell
[~,ix1]=ismember(pointClosure(:,1:2),IDX_index(:,1:2),'rows');
[~,ix2]=ismember(pointClosure(:,2:3),IDX_index(:,1:2),'rows');
[~,ix3]=ismember(pointClosure(:,[1 3]),IDX_index(:,1:2),'rows');
arcTriClosureIndex=[IDX_index(ix1,3) IDX_index(ix3,3) IDX_index(ix2,3)];
end
function ft = addTriangle(edgeInParade,ft,xyfull,dismat)
    for i = 1:size(edgeInParade,1)
        nodefrom = edgeInParade(i,1);nodeto = edgeInParade(i,2);
        %find all connected lines with this line
        connectedLines = find(ft(:,1) == nodefrom | ft(:,1) == nodeto ...
            | ft(:,2) == nodefrom | ft(:,2) == nodeto);
        connectedLines = ft(connectedLines,:);
        %delete itself
        connectedLines(ismember(connectedLines,edgeInParade(i,:),'rows')) = [];
        NodesInLines = unique(connectedLines(:),'rows');
        NodesInLines(NodesInLines == nodefrom | NodesInLines == nodeto) = [];
        connectivityToDo1 = sort([nodefrom*ones(size(NodesInLines)),NodesInLines],2);
        connectivityToDo2 = sort([nodeto*ones(size(NodesInLines)),NodesInLines],2);
        if isempty(connectivityToDo1)
            continue;
        end
%          fprintf('Find Local %d Triangles... \n',length(connectivityToDo1));
        disToDo1 = dismat(ismember(xyfull,connectivityToDo1,'rows'));
        disToDo2 = dismat(ismember(xyfull,connectivityToDo2,'rows'));
        disTotal = disToDo1 + disToDo2;
        locateToDo = find(disTotal == min(disTotal));
        addedge1 = connectivityToDo1(locateToDo(1),:);
        addedge2 = connectivityToDo2(locateToDo(1),:);
        ft = [ft;addedge1;addedge2];
    end
end
