function shortPhase = shorten_phase_descriptor(cat, thisPhase, arrivalTime)

    plength = length(thisPhase);

    shortened = false;
    
    %shorten phase name if it's the direct P
    if arrivalTime == 0
        shortPhase = 'P';
        shortened = true;
        return
    end
    
    %shorten phase name if it's the direct S
    if cat.table.DistFromSlabTop < -8
        if strcmp(thisPhase,'S_mS_tS') || strcmp(thisPhase,'S_mS_tS_MS')
            shortPhase = 'S';
            shortened = true;
            return
        elseif strcmp(thisPhase,'S_mS_tS_MP')
            shortPhase = 'S_MP';
            shortened = true;
            return
        elseif strcmp(thisPhase,'P_mP_tP_MS')
            shortPhase = 'P_MS';
            shortened = true;
            return
        end
    elseif cat.table.DistFromSlabTop < 0
        if strcmp(thisPhase,'S_tS') || strcmp(thisPhase,'S_tS_MS')
            shortPhase = 'S';
            shortened = true;
            return
        elseif strcmp(thisPhase,'S_tS_MP')
            shortPhase = 'S_MP';
            shortened = true;
            return
        elseif strcmp(thisPhase,'P_tP_MS')
            shortPhase = 'P_MS';
            shortened = true;
            return
        elseif strcmp(thisPhase,'P_mP_tP_MP') || strcmp(thisPhase,'P_mP_tP')
            shortPhase = 'P_mP';
            shortened = true;
            return
        elseif strcmp(thisPhase,'S_mS_tS_MS') || strcmp(thisPhase,'S_mS_tS')
            shortPhase = 'S_mS';
            shortened = true;
            return
        elseif strcmp(thisPhase,'P_mS_tS_MS') || strcmp(thisPhase,'P_mS_tS')
            shortPhase = 'P_mS';
            shortened = true;
            return
        elseif strcmp(thisPhase,'P_mS_tS_MS') || strcmp(thisPhase,'P_mS_tS')
            shortPhase = 'P_mS';
            shortened = true;
            return
        end
    elseif cat.table.DistFromSlabTop >= 0
        if strcmp(thisPhase,'P_tP_mP_tP_MP') || strcmp(thisPhase,'P_tP_mP_tP')
            shortPhase = 'P_mP';
            shortened = true;
            return
        elseif strcmp(thisPhase,'S_tS_mS_tS_MS') || strcmp(thisPhase,'S_tS_mS_tS')
            shortPhase = 'S_mS';
            shortened = true;
            return
        elseif strcmp(thisPhase,'P_tP_mS_tS_MS') || strcmp(thisPhase,'P_tP_mS_tS')
            shortPhase = 'P_mS';
            shortened = true;
            return
        elseif strcmp(thisPhase,'S_tS_mP_tP_MP') || strcmp(thisPhase,'S_tS_mP_tP')
            shortPhase = 'S_mP';
            shortened = true;
            return
        end
    end
    
    
    % shorten the last letters if there's no interaction with the
    % overriding Moho
    if ~shortened
        if plength >= 4
            if strcmp(thisPhase(end-3:end),'P_MP')
                shortPhase = thisPhase(1:end-3);
                shortened = true;
            elseif strcmp(thisPhase(end-3:end),'S_MS')
                shortPhase = thisPhase(1:end-3);
                shortened = true;
            end
        end
    end
     
    
    %return long name if it couldn't shorten
    if ~shortened
        shortPhase = thisPhase;
    end   

        