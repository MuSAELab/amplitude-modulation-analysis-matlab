function w = cosinewin(n)
%COSINEWIN Returns a cosine window with N elements 
n_vct = (0.5 : 1 : n-0.5)';
w = sin(pi* n_vct ./ n );
end