function newArrival = subsetArrival(arrivals, id)

    newArrival = struct();
    if length(id)==0
        return 
    end
    
    newArrival.phase = arrivals.phase{id};
    newArrival.time = arrivals.time(id);
    newArrival.nReflections = arrivals.nReflections(id);
    newArrival.nConv = arrivals.nConv(id);