function w = load_miniseed(request)
   %LOAD_MINISEED loads a waveform from MINISEED files
   % combineWaves isn't currently used!
   
   % Glenn Thompson 2016/05/25 based on load_sac
   % request.combineWaves is ignored
   
   if isstruct(request)
      [thisSource, chanInfo, startTime, endTime, ~] = unpackDataRequest(request);
      filenamelist={};
      
      % Work out which files we need
      for i=1:numel(chanInfo)
            thisfilename = getfilename(thisSource,chanInfo(i),startTime);
            found=false;
            listlength = numel(filenamelist);
            for c=1:listlength
                if strcmp(thisfilename,filenamelist{c})
                    found=true;
                end
            end
            if ~found
                filenamelist{listlength+1} = thisfilename;
                listlength=listlength+1;
            end
      end
      wfiles = [];
      
      % Load waveforms from all these files
      for c=1:numel(filenamelist)
         wtmp = [];
         try
            wtmp = mseedfilename2waveform(thisfilename{1}, startTime, endTime);
         catch ME
             warning(ME.message)
             continue
         end
         wtmp = reshape(wtmp, [1 numel(wtmp)]);
         wfiles = [wfiles wtmp];
      end
      if isempty(wfiles)
          w = waveform();
          return
      end
      w = combine(wfiles);
      
      % Extract based on time
      w = extract(w, 'time', startTime, endTime);
      
      % Extract based on ChannelTag
      %w = matchChannelTag(w);
      
   else
      %request should be a filename
      thisFilename = request;
      if exist(thisFilename, 'file')
        w = mseedfilename2waveform(thisFilename);
      else
          w = waveform();
          warning(sprintf('File %s does not exist',thisFilename));
      end
   end
end


function w = mseedfilename2waveform(thisfilename, snum, enum)
    read_with_obspy = false;
    read_with_matlab = false;
    w = waveform();

    try
        s = ReadMSEEDFast(thisfilename); % written by Martin Mityska
     % s = rdseed(thisfilename); % written by Martin Mityska
        read_with_matlab = true;
    catch ME
        % With blockette length error, (e.g., 'Product of known dimensions, 4096, not divisible 
        % into total number of elements, 277504.' - try to read with obspy instead:
        env_info = pyenv;
        % if (strcmp(ME.identifier, 'MATLAB:getReshapeDims:notDivisible') || ...
        %     strcmp(ME.identifier, 'MATLAB:badsubscript')) && ...
        %         strcmp(env_info.Status, "Loaded")
        if ~read_with_matlab && strcmp(env_info.Status, "Loaded")
            % Read with obspy
            stream = py.obspy.read(thisfilename);
            read_with_obspy = true;
            s = cell(stream.traces);
        elseif (contains(ME.message, "MATLAB:FileIO:InvalidFid"))
            warning(ME.message)
            w = waveform();
            return
        else
            throw(ME)
        end
    end
    

    for c=1:numel(s)
        if read_with_obspy
            trace = s{c};
            try
                dtype = string(trace.data.dtype.name);
                switch dtype
                    case "int32"
                        data = int32(trace.data.tolist);
                    case "float32"
                        data = double(trace.data.tolist);
                    otherwise
                        warning("Datatype from obspy miniseed trace not understood")
                end
            catch
                warning("Cannot convert data in file %s read with obspy to Matlab dtype", ...
                        thisfilename);
                warning(ME.message)
                return
            end

            w(c, 1) = waveform( ...
                ChannelTag(char(trace.stats.network), char(trace.stats.station), ...
                           char(trace.stats.location), char(trace.stats.channel)), ...
                trace.stats.sampling_rate, epoch2datenum(trace.stats.starttime.timestamp), ...
                data);
        else
            w(c,1) = waveform( ...
                ChannelTag(s(c).network, s(c).station, s(c).location, s(c).channel), ...
                s(c).sampleRate, epoch2datenum(s(c).startTime), s(c).data);
        end
     end
end