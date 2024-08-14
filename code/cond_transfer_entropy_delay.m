%% conditional TE of symbolized time series
function TY_X_condZ = cond_transfer_entropy_delay(piXall,piYall,piZall,del,del2,flag_norm) 
% This func calculates TE_Y_to_X | Z = I( X_k : Y_k-delta | X_k-1 , Z_k-delta2 ) 
% --> mutual information between X and past Y given past X is known
% del = delay of y (source) ts in terms of sampled time steps
% del2 = delay of z (conditional) ts in terms of sampled time steps
N = length(piXall);

if (del>=del2)
    % X_(k-1) = X1, Y_(k-delta) = Y1, X_k = X, Z_(k-delta2) = Z1
    piX1 = piXall(del:N-1,1); piY1 = piYall(1:N-del,1);
    piX = piXall(del+1:N,1); piZ1 = piZall(del-del2+1:N-del2,1);
else
    % X_(k-1) = X1, Y_(k-delta) = Y1, X_k = X, Z_(k-delta2) = Z1
    piX1 = piXall(del2:N-1,1); piY1 = piYall(del2-del+1:N-del,1);
    piX = piXall(del2+1:N,1); piZ1 = piZall(1:N-del2,1);
end

HX1Z1 = compute_entropy([piX1 piZ1]); % Joint entropy of X_(k-1),Z_(k-del2)
HX1Y1Z1 = compute_entropy([piX1 piY1 piZ1]); % Joint entropy of X_(k-1),Y_(k-del),Z_(k-del2)
HXX1Z1 = compute_entropy([piX piX1 piZ1]); % Joint entropy of X_(k),X_(k-1),Z_(k-del2)
HXX1Y1Z1 = compute_entropy([piX piX1 piY1 piZ1]); % Joint entropy of X_(k),X_(k-1),Y_(k-del),Z_(k-del2)

TY_X_condZ = (HXX1Z1 - HX1Z1) - (HXX1Y1Z1 - HX1Y1Z1);

if (flag_norm>0)
    TY_X_condZ = TY_X_condZ/(HXX1Z1 - HX1Z1);
end

end