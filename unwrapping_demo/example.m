%% example of interferogram stack unwrapping
load('test_data.mat')
CtSent_PhaseUnwrappingDenseGrid(lonlat,ph,pos,1);% 1 for graph-theory based unwrapping, 0 for ordinary MCF unwrapping
load('PhU_apsp.mat','phuw');
subplot(2,1,1);scatter(lonlat(:,1),lonlat(:,2),5,ph(:,1),'filled');title('input phase');
subplot(2,1,2);scatter(lonlat(:,1),lonlat(:,2),5,phuw(:,1),'filled');title('Graph-Theory unwrapped phase');
