function psi = LS_ESPRIT_Algorithm(Y_ll_bar, L)
    
    Y_sub1 = Y_ll_bar(1:end-1,:);
    Y_sub2 = Y_ll_bar(2:end,:);
    Y = [Y_sub1; Y_sub2];
    
    [U,~,~] = svds(Y, L);
    U_1 = U(1:size(Y,1)/2, :);
    U_2 = U(size(Y,1)/2+1:end, :);

%     [U_1, ~, ~] = svds(Y_sub1, L);
%     [U_2, ~, ~] = svds(Y_sub2, L);
    Psi = U_1\U_2;
    
    psi = eig(Psi);    
end