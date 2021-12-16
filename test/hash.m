function hash = hash(A)

    B = getByteStreamFromArray(A);
    md = java.security.MessageDigest.getInstance('SHA-1');
    md.update(B);
    hash = reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []);
    
end