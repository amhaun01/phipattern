%clear
%close all
function [NperGr] = qualePlot6b(phiList,Nper,threshold,mlim)%input should be in binary order!

if(nargin>=3)
    phiList(abs(phiList)<threshold) = NaN;
end;
if(nargin<4)
    mlim = max(phiList);
end;

N = size(Nper,2);

sumNper = sum(Nper,2);
%%
Nindex = (2:2^N)-1;
ssumNper = sortrows([sumNper Nper Nindex' phiList(:)],1);
NperGr = struct([]);
tls = 0;
tlevs = zeros(N,1);
%separate groups by order
for x = 1:N
    tlev = nchoosek(N,x);
    if(x<N)
    %NperGr(x).Nper = sortrows(ssumNper(tls+(1:tlev),2:N+1),-1:-1:-N)>0;
    %[NperGr(x).Nper,NperGr(x).phiList] = binmatsorter(ssumNper(tls+(1:tlev),2:N+1),ssumNper(tls+(1:tlev),end));
    NperGr(x).Nper = ssumNper(tls+(1:tlev),2:N+1)>0;
    NperGr(x).phiList = ssumNper(tls+(1:tlev),end);
    
    if(N==4)
        if(x==2)
            pdord = [1 0 1 0;
                     0 1 1 0;
                     0 0 1 1;
                     0 1 0 1;
                     1 0 0 1;
                     1 1 0 0]>0;
                 npL = [];
            for o = 1:length(NperGr(x).phiList)
                for ob = 1:length(NperGr(x).phiList)
                    if(pdord(o,:)==NperGr(x).Nper(ob,:))
                        npL(o) = NperGr(x).phiList(ob);
                    end;
                end;
            end; 
            NperGr(x).phiList = npL;
            NperGr(x).Nper = pdord;
        end;
        
        if(x==3)
            pdord = [1 1 1 0;
                     0 1 1 1;
                     1 0 1 1;
                     1 1 0 1]>0;
                 npL = [];
            for o = 1:length(NperGr(x).phiList)
                for ob = 1:length(NperGr(x).phiList)
                    if(pdord(o,:)==NperGr(x).Nper(ob,:))
                        npL(o) = NperGr(x).phiList(ob);
                    end;
                end;
            end; 
            NperGr(x).phiList = npL;
            NperGr(x).Nper = pdord;
        end;
    end;
    
    else
    NperGr(x).Nper = ssumNper(tls+(1:tlev),2:N+1)>0;
    NperGr(x).phiList = phiList(end);
    end;
    NperGr(x).tlev = tlev;
    tls = tls+tlev;
    tlevs(x) = tlev;
end;

%assign coordinates
for x = 1:N
    %if(x<N)
    for y = 1:NperGr(x).tlev
        pangle = -(y-1)*2*pi/NperGr(x).tlev + mod(x,2)*pi/NperGr(x).tlev + pi/5;
        NperGr(x).coords(y,1:2) = (N - x)*...
                                [cos(pangle) sin(pangle)];
        NperGr(x).coords(y,3) = NperGr(x).phiList(y);
    end;  
    %elseif(x==N)
    %    NperGr(x).coords = [0 0 NperGr(x).phiList];
    %end;
end;
%%
pcoords = [];
pNper = [];
vlist = [];
for n = 2:length(NperGr)
    pcoords = [pcoords; NperGr(n).coords];
    pNper = [pNper;NperGr(n).Nper];
    vlist = [vlist;NperGr(n).phiList'];
end;

for x = 1:length(pNper)
    ll = 'ABCDEFGHI';
    lls{x} = ll(pNper(x,:)>0);
end;
pcoords = pcoords(:,[1 2 3]);
%%
Ns = sum(pNper,2);
dmat = squareform(pdist(pNper,'Euclidean'));
%pcoords = mdscale(dmat,2);pcoords(:,3) = vlist;%pcoords(:,2) = Ns;

amat = squareform(pdist(pNper,'Cityblock'))==1;

[is,js] = find((amat));
Xs = [pcoords(is,1) pcoords(js,1)];
Ys = [pcoords(is,2) pcoords(js,2)];
Zs = [pcoords(is,3) pcoords(js,3)];

cols = 1-min(vlist./mlim,1);
cols = .8*(cols.^.5);

for x = 1:length(is)
    ecols(x) = max(cols([is(x) js(x)]));
end;

for x = 1:length(is)
    plot3(Xs(x,:),Ys(x,:),Zs(x,:),'-','Color',ecols(x)+[0 0 0])
    if(x==1)
        hold on
    end;
end;
%scatter3(pcoords(:,1),pcoords(:,2),pcoords(:,3),100,Ns,'Filled')
hold off
for x = 1:length(vlist)
text(pcoords(x,1),pcoords(x,2),pcoords(x,3),lls{x},...
    'FontSize',12,'Color',cols(x)*[1 1 1] + 0*(1-cols(x))*[0 0 0],...
    'HorizontalAlignment','center','FontWeight','bold')
end;
%axis(1.2*max(abs(pcoords(:,1)+i*pcoords(:,2)))*[-1 1 -1 1])
axis square
%axis off
set(gcf,'Color','w')
