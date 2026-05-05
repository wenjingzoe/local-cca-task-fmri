function s = block_stimulus(n,block_len)
% Deterministic block stimulus
s = zeros(n,1);
for i = 1:n
    if mod(floor((i-1)/block_len),2)==1
        s(i)=1;
    end
end
end
