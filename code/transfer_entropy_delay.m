%% TE based on symbolic approach
function TY_X = transfer_entropy_delay(piXall,piYall,del,flag_norm) 
% This func calculates TE_Y_to_X = I( X_k : Y_k-1 | X_k-1 ) 
% --> mutual information between X and past Y given past X is known
% del = delay in terms of sampled time steps

% X_(k-1) = X1, Y_(k-1) = Y1, X_k = X, Y_k = Y
piX1 = piXall(del:end-1,1); piY1 = piYall(1:end-del,1);
piX = piXall(del+1:end,1); %Y = Yall(2:end,1);

HX1 = compute_entropy(piX1); % Entropy of X_(k-1)
HX1Y1 = compute_entropy([piX1 piY1]); % Joint entropy of X_(k-1),Y_(k-1)
HXX1 = compute_entropy([piX piX1]); % Joint entropy of X_(k),X_(k-1)
HXX1Y1 = compute_entropy([piX piX1 piY1]); % Joint entropy of X_(k),X_(k-1),Y_(k-1) 

TY_X = (HXX1 - HX1) - (HXX1Y1 - HX1Y1);

if (flag_norm>0)
    TY_X = TY_X/(HXX1 - HX1);
end

end