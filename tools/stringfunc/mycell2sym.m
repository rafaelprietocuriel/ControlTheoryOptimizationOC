function symval=mycell2sym(cellvec,matrixtype)
% each cell for matrixtype 'matrix' consists of a row of the matrix
% each cell for matrixtype 'vector' consists of an entry of the vector

symval=sym([]);
l=length(cellvec);
switch matrixtype
    case 'vector'
        stringval=['[' deblank(cellvec{1})];
        for ii=2:l
            stringval=[stringval ';' deblank(cellvec{ii})];
        end
        symval=mystr2sym([stringval ']']);
    case 'matrix'
        stringval=['[' deblank(cellvec{1})];
        for ii=2:l
            stringval=[stringval ';' deblank(cellvec{ii})];
        end
        symval=mystr2sym([stringval ']']);
end
