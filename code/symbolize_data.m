% Symbolize time series X using embedding dimension m
function piX = symbolize_data(X,m) % Only for m=2,3,4

%% m = 2 case
if (m==2) 
    piX = diff(X)>0;
elseif (m==3)
    X0 = X(1:end-2); % X_t
    X1 = X(2:end-1); % X_t+1
    X2 = X(3:end); % X_t+2
    
    [~,symX] = sort([X0 X1 X2],2); % get m-history symbols in the rows of symX
    symX = symX-1;
    
    [~,~,piX] = unique(symX,'rows'); % assign numbers to each type of symbol
elseif (m==4)
    X0 = X(1:end-3); % X_t
    X1 = X(2:end-2); % X_t+1
    X2 = X(3:end-1); % X_t+2
    X3 = X(4:end); % X_t+3

    [~,symX] = sort([X0 X1 X2 X3],2); % get m-history symbols in the rows of symX
    symX = symX-1;
    
    [~,~,piX] = unique(symX,'rows'); % assign numbers to each type of symbol
else
    disp('We have not coded for m>4 yet');
end