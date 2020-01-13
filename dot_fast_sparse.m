function outVec = dot_fast_sparse(vec1, vec2)

% Inputs vec1 and vec2 are assumed to be p x N arrays, where p is the
% vector length, and N is the number of vectors
%
% Output outVec is then a 1 x N array where each entry is the dot product
% of vec1(:,ii) and vec2(:,ii)

N = size(vec1,2);
if size(vec2,2) ~= N
    error('Inputs don''t have the same number of vectors');
end

p = size(vec1,1);
if size(vec2,1) ~= p
    error('Inputs don''t have the same length of vectors');
end

vec1_sparse = sparse(repmat([1:N],p,1),[1:p*N],vec1);

vec2_sparse = sparse([1:p*N],repmat([1:N],p,1),vec2);

outVec = full(diag(vec1_sparse*vec2_sparse));

end