function output = ndgrid2(input)

    d = numel(input);
    one_vec = ones(1, d);
    siz = cellfun(@numel, input);
    output = zeros(prod(siz), d);
    for k = 1:d
        x = input{k};
        s = one_vec;
        s(k) = numel(x);
        x = reshape(x, s);
        s = siz;
        s(k) = 1;
        x = repmat(x, s);
        output(:, k) = x(:);
    end
    
end