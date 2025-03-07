function arcarg=constraintcombinationindex2arcarg(ocObj,idx)
arcarg=[];
identifier=inequalitycontrolconstraintidentifier(ocObj);
actccomb=identifier(idx==1);

counter=0;
while counter<arcnum(ocObj)
    ccomb=constraintcombination(ocObj,counter);
    try
        if all(strcmp(actccomb,ccomb))
            arcarg=counter;
            break
        end
    end
    counter=counter+1;
end