%% Adjoint calculator
% This function calculates the adjoint given simply the transform between
% joints

function adj = adj_calc(g,inv)
    % convert transform to rotation and translation
    R = g(1:3,1:3);
    p = g(1:3,4);
    % get skew symmetric matrix of translation
    p_hat = skew(p);
    if inv == 0
        % package into adjoint
        adj = [R p_hat*R; zeros(3) R];
    else
        adj = [R.' -R.'*p_hat; zeros(3) R.'];
    end
    
end