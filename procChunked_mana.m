clear
close all

%spath is location of this script, the getPhiStructNP script, and the
%compfiles directory
spath = pwd;
%dpath is location of the data_d1 file
dpath = pwd;

addpath(genpath([spath '/phi_toolbox_Feb2014/']));
%%
SI = 'L110519';
rls = {'pul','teo','lip','v4'};
conds = {'LFP_cueexRFp4'
         'LFP_cueinRFp2'
         'LFP_arrexRFp4'
         'LFP_arrinRFp2'};

load([dpath '/data_d1.mat']);     
%%
figure(1)
for cond = 1:4 
    cnames{cond} = conds{cond}(5:end);
    ts{cond} = 1:size(datmat{cond},2);
    
    subplot(2,2,cond)
    plot(ts{cond},mean(datmat{cond},3))
    legend(rls)
    title(cnames{cond})
    ylim(.2*[-1 1])
end;
%%

nChan = 4;
eSets = geteSets(nChan);

lls = {'Pu','Te','Li','V4'};
for x = 1:length(eSets)
    labs{x} = [lls{eSets(x,:)}];
end;
%%
Perms = size(eSets,1);
Psizes = sum(eSets,2);

PartStruct = struct([]);
for Pind = Psizes(1):Psizes(end)
    load(strcat([spath '/CompFiles/ArrayPartitions/arrayPartitions_',num2str(Pind),'.mat']));

    if(length(partitions)>1)
        partitions = sortrows([partitions max(partitions,[],2)],Pind+1);
    end;
    PartStruct(Pind).Partitions = partitions(:,1:Pind);
end;
%%
%%
ssize = 100;
taus = 8;
sstep = 10;
%maxpart = 2;%check partitions up to maxpart partitions (2 = check only bipartitions)

for cond = 1:4

sps{cond} = 200:sstep:ts{cond}(end)-ssize-200;

od(cond).phis = [];
od(cond).atomphis =[];
od(cond).MIs = [];
od(cond).Hs = [];
for S = 1:length(sps{cond})
    sp = sps{cond}(S);
    tvec = (1:ssize) + sp;
    data = datmat{cond}(:,tvec,:);
    
    shufdata = [];
    for x = 1:4
        shufdata(x,:,:) = data(x,:,randperm(size(data,3)));
    end;

    minPhi = zeros(length(eSets),length(taus));
    atomPhi = zeros(length(eSets),length(taus));
    MI = zeros(length(eSets),length(taus));
    Hcond = zeros(length(eSets),length(taus));
    for tind = 1:length(taus)
    tau = taus(tind);
    [minPhi(:,tind),minPart,MI(:,tind),Hcond(:,tind),atomPhi(:,tind),mPlist] = ...
        getPhiStructNP(data,eSets,PartStruct,tau);
    end;
    H = Hcond+MI;
    minPhi(1:nChan,:) = NaN;

    od(cond).phis(:,S,:) = minPhi;
    od(cond).atomphis(:,S,:) = atomPhi;
    od(cond).MIs(:,S,:) = MI;
    od(cond).Hs(:,S,:) = H;
end;
end;
%%

for cond = 1:4
ccs = repmat(sum(eSets,2),[1 length(sps{cond})]);
for tau = taus
figure(10)
subplot(2,2,cond)
imagesc(sps{cond}-sontimes(cond),1:length(eSets),od(cond).phis(:,:,taus==tau))
hold on
plot(sontimes(cond)+[0 0]-ssize/4-sontimes(cond),[0 1+length(eSets)],'w--')
hold off
axis xy
%caxis([0 .08])
set(gca,'YTick',1:length(labs),'YTickLabel',labs)
%colorbar
title(cnames{cond})
end;
end;
%%
%tau = 7;
allphs = [];
ps = [];
for cond = 1:4
    ccs = repmat(sum(eSets,2),[1 length(sps{cond})]);
    allphs = [allphs od(cond).phis(:,:,taus==tau)./ccs];
    ps = [ps cond+zeros(1,size(od(cond).phis,2))];
end;
dmat = pdist(allphs(1+nChan:15,:)','cityblock');
pxys = cmdscale(dmat);

cols = [.8 .4 .2;
        1 0 0;
        .2 .4 .8;
        0 0 1];
lts = {'.--','.-','.--','.-'};
figure(26)
for cond = 1:4
   pxcs{cond} = pxys(ps==cond,:);   
   
   plot(pxcs{cond}(:,1),...
        pxcs{cond}(:,2),lts{cond},'Color',cols(cond,:),'LineWidth',2)
    text(pxcs{cond}(:,1),...
        pxcs{cond}(:,2),num2str(sps{cond}'-sontimes(cond)+ssize/2),...
        'FontSize',10,'Color',cols(cond,:))
    if(cond==1)
    hold on;
    end;
end;
hold off

legend(cnames,'FontSize',16,'Location','NorthWest')

figure(27)
for cond = 1:4
    tls{cond} = sps{cond}'-sontimes(cond)+ssize/2;
    
    subplot(2,2,cond)
   pxcs{cond} = pxys(ps==cond,:);   
   
   plot(pxcs{cond}(:,1),...
        pxcs{cond}(:,2),lts{cond},'Color',cols(cond,:),'LineWidth',2)
    text(pxcs{cond}(:,1),...
        pxcs{cond}(:,2),num2str(sps{cond}'-sontimes(cond)+ssize/2),...
        'FontSize',10,'Color',cols(cond,:))
    title(cnames{cond})
    axis(1.1*[minmax(pxys(:,1)') minmax(pxys(:,2)')])
end;
%%
pois = [0 90 600;
        0 90 600;
        -470 0 100;
        -470 0 100];
    
    ms = [.025 .025 .035 .035;
          .06 .06 .025 .025;
          .05 .05 .1 .1];

    figure(25)
    for p = 1:size(pois,2)
        for cond = 1:4
            pind = find(tls{cond}==pois(cond,p));
            
            phiList = od(cond).phis(:,pind)./Psizes;
            
            subplot(size(pois,2),4,cond + 4*(p-1))
            qualePlot6d(phiList,eSets,0,max(phiList),lls);
            view([0 0])
            title({cnames{cond},['at ' num2str(pois(cond,p)) 'SOA']},'FontSize',16)
            zlim([0 1.2*max(phiList)])
            xlim(2.5*[-1 1])
            ylim(2.5*[-1 1])
        end;
    end;
%%
    figure(25)
    for p = 1:size(pois,2)
        for CX = [1 2]
            pind = find(tls{2*CX}==pois(2*CX,p));
            
            phiList = od(2*CX).phis(:,pind)./sum(eSets,2) - ...
                      od(2*CX - 1).phis(:,pind)./sum(eSets,2);
            
            subplot(size(pois,2),2,CX + 2*(p-1))
            qualePlot6d(phiList,eSets,0,max(phiList),lls);
            view([0 0])
            title({cnames{2*CX},['at ' num2str(pois(2*CX,p)) 'SOA']},'FontSize',16)
            zlim(1.2*max(max(abs(phiList)),.015)*[-1 1])
            xlim(2.5*[-1 1])
            ylim(2.5*[-1 1])
        end;
    end;
%%
cs = [2 1;4 3];
figure(12)
for c = 1:2
subplot(3,2,c)
imagesc(sps{cs(c,1)},1:length(eSets),...
    od(cs(c,1)).phis(:,:,taus==tau)-od(cs(c,2)).phis(:,:,taus==tau))
%caxis(.1*[-1 1])
axis xy
set(gca,'YTick',1:length(labs),'YTickLabel',labs,'FontSize',12)
hold on
plot(sontimes(cs(c,1))+[0 0]-ssize/4,[0 1+length(eSets)],'w--')
hold off

subplot(3,2,c+2)
imagesc(sps{cs(c,1)},1:length(eSets),...
    od(cs(c,1)).MIs(:,:,taus==tau)-od(cs(c,2)).MIs(:,:,taus==tau))
%caxis(.5*[-1 1])
axis xy
set(gca,'YTick',1:length(labs),'YTickLabel',labs,'FontSize',12)
hold on
plot(sontimes(cs(c,1))+[0 0]-ssize/4,[0 1+length(eSets)],'w--')
hold off

subplot(3,2,c+4)
imagesc(sps{cs(c,1)},1:length(eSets),...
    od(cs(c,1)).Hs(:,:,taus==tau)-od(cs(c,2)).Hs(:,:,taus==tau))
%caxis(1*[-1 1])
axis xy
set(gca,'YTick',1:length(labs),'YTickLabel',labs,'FontSize',12)
hold on
plot(sontimes(cs(c,1))+[0 0]-ssize/4,[0 1+length(eSets)],'w--')
hold off
end;
