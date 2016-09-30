function eSets = geteSets(Nelec,porder)

if(nargin==1)
    porder = Nelec;
end;

%eSetsA = zeros(2^Nelec,Nelec);
for x = 1:(2^Nelec)
    xs = x;
    for y = 1:Nelec
        eSetsA(x,y) = xs - 2*floor(xs/2);
        xs = floor(xs/2);
    end;
    if(sum(eSetsA(x,:)>porder))
        break
    end;
end;
eSetsA = sortrows([eSetsA sum(eSetsA,2)],Nelec+1);
eSetsA = eSetsA(:,1:Nelec);

eSets = eSetsA(sum(eSetsA,2)>0,:)>0;