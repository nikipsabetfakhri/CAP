function M2 = ConcatMat(M,n_pop,n_states,n_subjects,type)
% Concatenates populations appropriately for Violin plotting
% Matches GUI behavior in CAP_TB.m

% Creates the data matrix (NaNs pad unequal lengths)
M2 = nan(n_pop*n_states, max(cell2mat(n_subjects)));

for i = 1:n_pop
    switch type
        case 'Raw counts'
            tmp = M{i}.raw.state(:,1:n_states)';
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
        case 'Normalized counts'
            tmp = M{i}.frac.state(:,1:n_states)';
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
        case 'Number'
            tmp = M{i}(:,3:3+n_states-1)';
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
        case 'Duration'
            tmp = DiscardNaN(M{i}(:,3:3+n_states-1))';
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
        case {'Betweenness','kin','kout','Resilience'}
            tmp = M{i}';
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
        case 'FD'
            tmp = M{i};
            for j = 1:n_states
                M2(i+(j-1)*n_pop,1:size(tmp,2)) = tmp(j,:);
            end
    end
end
end

function M2 = DiscardNaN(M)
M2 = [];
for i = 1:size(M,1)
    if ~any(isnan(M(i,:)))
        M2 = [M2; M(i,:)];
    end
end
end