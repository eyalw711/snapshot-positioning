function epochs = ubx_reader_codephase(file)
fid = fopen(file,'rt');

non_unique_epochs = [];
in_epoch = false;
obs_inx = 1;
while 1 %epochs loop
    line = fgetl(fid);
    if (line == -1)
        break;
    end
    
    if ~in_epoch
        answer = findstr(line,'version 1');
        
        if  ~isempty(answer)
            curr_epoch = struct();
            line = fgetl(fid);
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            curr_epoch.recvTOW = str2double(strip(tok))*1e-3;
            curr_epoch.obs = struct();
            curr_epoch.obs.sv = [];
            curr_epoch.obs.cp = [];
            curr_epoch.obs.dopl = [];
            in_epoch = true;
            continue;
        end
        
    else
        
        answer = findstr(line,'UBX-RXM-MEASX');
        if  ~isempty(answer)
            non_unique_epochs = [non_unique_epochs; curr_epoch];
            in_epoch = false;
            obs_inx = 1;
            continue;
        end
        
        answer = findstr(line,'gnssId');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_gnssId = str2num(strip(tok));
            continue;
        end
        
        answer = findstr(line,'svId');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_svId = str2num(strip(tok));
            continue;
        end
        
        answer = findstr(line,'dopplerHz');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_dopl = str2num(strip(tok));
            continue;
        end
        
        answer = findstr(line,'dopplerMS');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_doplMS = str2num(strip(tok));
            continue;
        end
        
        answer = findstr(line,'codePhase');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_cpMes = str2double(strip(tok));
            
            if latest_gnssId == 0
                curr_epoch.obs.sv(obs_inx) = latest_svId;
                curr_epoch.obs.cp(obs_inx) = latest_cpMes;
                curr_epoch.obs.dopl(obs_inx) = latest_dopl;
                curr_epoch.obs.doplMS(obs_inx) = latest_doplMS;
                obs_inx = obs_inx + 1;
            end
            
            continue;
        end
    end
end

fclose(fid);

epochs = [];
for i = 1:numel(non_unique_epochs)
    
    epoch = non_unique_epochs(i);
    
    new_epoch = struct();
    new_epoch.recvTOW = epoch.recvTOW;
    new_epoch.obs = struct();
    new_epoch.obs.sv = [];
    new_epoch.obs.cp = [];
    new_epoch.obs.doplMS = [];
    
    unique_svs = unique(epoch.obs.sv);
    for svnum = unique_svs
        tmp_epoch = epoch;
        tmp_epoch.obs.dopl(epoch.obs.sv ~= svnum) = 0;
        [~,imax] = max(abs(tmp_epoch.obs.dopl));
        new_epoch.obs.sv = [new_epoch.obs.sv svnum];
        new_epoch.obs.cp = [new_epoch.obs.cp epoch.obs.cp(imax)];
        new_epoch.obs.doplMS = [new_epoch.obs.doplMS epoch.obs.doplMS(imax)];
    end
    
    epochs = [epochs; new_epoch];
end
end

