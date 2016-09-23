function [minPhi,minPart,compMI,entropy,atomPhi,mPlist,Cov_Xs] = getPhiStructNP(data,eSets,PartStruct,tau,maxPart)
%shape of data should be channels/tsteps/trials

if(nargin<5)
    %to speed processing, only look at partitions up to maxPart
    maxPart = size(eSets,2);
end;

eSetNs = sum(eSets,2);
if(size(eSets,1)==1)
    eSets = eSets(eSets>0);
end;

minPart = zeros(size(eSets));
minPhi = zeros(length(eSetNs),1);
compMI = zeros(length(eSetNs),1);
atomPhi = zeros(length(eSetNs),1);
entropy = zeros(length(eSetNs),1);
mPlist = zeros(length(eSetNs),1);

%%
Nelec = size(data,1);
Trs = size(data,3);

Cov_Xs = zeros(Nelec,Nelec,Trs);
Cov_XYs = zeros(Nelec,Nelec,Trs);
Cov_YXs = zeros(Nelec,Nelec,Trs);
Cov_Ys = zeros(Nelec,Nelec,Trs);
for tr = 1:Trs
    [Cov_Xs(:,:,tr), Cov_XYs(:,:,tr), Cov_YXs(:,:,tr), Cov_Ys(:,:,tr)] = ...
        Cov_comp_shrink(data(:,:,tr),tau);
end;

Cov_Xt = mean(Cov_Xs,3);
Cov_XtXtau = mean(Cov_XYs,3);
%Cov_XtauXt = mean(Cov_YXs,3);
Cov_Xtau = mean(Cov_Ys,3);

for V = 1:length(eSetNs)
    ssize = eSetNs(V);
    usePart = PartStruct(ssize).Partitions;
    usePart = usePart(max(usePart,[],2)<=maxPart,:);
    
    plist = 1:size(usePart,1);
    
    Pval = zeros(size(usePart,1),2);
    MI = zeros(size(usePart,1),1);
    H_cond = zeros(size(usePart,1),1);
    phi_MI = zeros(size(usePart,1),2);
    phi_H = zeros(size(usePart,1),2);
    
    for x = 1:size(usePart,1)
        Z = usePart(x,:);
        
        [Pval(x,:),...
            ~,...
            phi_MI(x,:),...
            ~,...
            MI,...
            ~,...
            phi_H(x,:),...
            H_cond,...
            beta] = phi_compNoFixed(Cov_Xt(eSets(V,:),eSets(V,:)),...
            Cov_XtXtau(eSets(V,:),eSets(V,:)),...
            Cov_Xtau(eSets(V,:),eSets(V,:)),...
            Cov_Xtau(eSets(V,:),eSets(V,:)),Z);
    end;
    Pval(1,:) = [inf inf];
    
    minInd = plist(Pval(:,2)==min(Pval(:,2)));
    minPart(V,eSets(V,:)) = usePart(minInd(1),:);
    mPlist(V) = minInd(1);
    minPhi(V) =                Pval(minInd(1),1);
    compMI(V) =                      MI;
    entropy(V) =              H_cond;
    %phi_MI(V) =            phi_MI(phi_MI(:,2)==min(phi_MI(:,2)),1);
    %phi_H(V) =               phi_H(phi_H(:,2)==min(phi_H(:,2)),1);
    atomPhi(V) = Pval(end,1);
end;