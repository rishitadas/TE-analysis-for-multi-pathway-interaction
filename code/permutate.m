%% Shuffle x keeping yz dynamics

function [X,Y,Z] = permutate(x,y,z)

zz = [y,z];
[l,r,dyn]=unique(zz,'rows');
ldyn = [dyn(r),l];

idx = (1:1:length(x))'; 
X=[]; Y= []; Z= [];
for il=1:length(l)
    partidx = idx(dyn==ldyn(il,1));
    X = [X; [partidx, shuffle(x(dyn==ldyn(il,1)))]];
    Y = [Y; ldyn(il,2)*ones(size(partidx,1),1)]; 
    Z = [Z; ldyn(il,3)*ones(size(partidx,1),1)];
end
[~,I] = sort(X(:,1));
X = X(I,2);
Y = Y(I);
Z = Z(I);
end


% Shuffle vector 
function shuffA = shuffle(A)
    shuffA = A(randperm(size(A,1))); 
end
