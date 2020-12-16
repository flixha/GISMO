function theyAre = are_reflec_and_convers_oneInteract(thisPhase)

    theyAre = true;

    IFinteract{1} = strfind(thisPhase,'_m');
    IFinteract{2} = strfind(thisPhase,'_t');
    IFinteract{3} = strfind(thisPhase,'_M');
    IFinteract{4} = strfind(thisPhase,'_T');
    
    % if the ray crosses the same interface twice, then the latter crossing
    % has to be as a different wave type (P or S) than the former for it to
    % be a reflection. Otherwise it's not a reflection.
    for j=1:1:length(IFinteract)
        nIFinterAc = length(IFinteract{j});
        if nIFinterAc == 2
            if strcmp(thisPhase(IFinteract{j}(1) - 1),...
                    thisPhase(IFinteract{j}(2) - 1))
                theyAre = false;
            else
                theyAre = true;
            end
        end
    end
    
    if theyAre == true
        % start at the bottom. if there is an interaction with that interface,
        % check whether that interaction is a conversion 
        if contains(thisPhase,'_m')
            fstr = strfind(thisPhase,'_m');
            before = thisPhase(fstr(1) - 1);
            after = thisPhase(fstr(1) + 2);
            if strcmp(before, after)
                theyAre = false;
            else
                theyAre = true;
            end
        elseif contains(thisPhase,'_t')
            fstr = strfind(thisPhase,'_t');
            before = thisPhase(fstr(1) - 1);
            after = thisPhase(fstr(1) + 2);
            if strcmp(before, after)
                theyAre = false;
            else
                theyAre = true;
            end
        elseif contains(thisPhase,'_M')
            fstr = strfind(thisPhase,'_M');
            before = thisPhase(fstr(1) - 1);
            after = thisPhase(fstr(1) + 2);
            if strcmp(before, after)
                theyAre = false;
            else
                theyAre = true;
            end
        end
    end