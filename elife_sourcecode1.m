% This matlab code uses spike times of each cell relative to the stimulus, 
% and estimates relevant parameters for obtaining tuning curve, tuning
% width and tuning strength, of individual DSGCs and DSGC population. 
% The "datarun" structure is provided separately, which contains the spike
% times, sampling rate and trigger times to synchronize the spikes to the
% stimulus. 

% Input arguments: 
%                 datarun: datarun structure that can be added to memory by
%                 loading processed data (.mat) files
%                 savedstats: path to already saved stats (.mat file)
%                                                    
% Output arguments: 
%                 retstat: return analyzed ds stat (needs to be saved as
%                 .mat file)
% 
% Requirements: 
%                 (1) Statistics and Machine Learning toolbox
%                 (2) Curve Fitting toolbox
%                 (3) Circular Statistics toolbox (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-)
%                 (4) v2struct.m function   (https://www.mathworks.com/matlabcentral/fileexchange/31532-pack---unpack-variables-to---from-structures-with-enhanced-functionality?focused=3847342&tab=example)
%                 

% Copyright 2018, Suva Roy
%--------------------------------------------------------------------------

%% Note: save retstat as 'DSStat.mat'

function [retstat] = DSStatAggregator(datarun, savedstats)


    % Read in saved data file if it exists 
    if exist(savedstats,'file')~=0
        m = matfile(savedstats); 
        oldStat = m.DSStat; 
    else
        oldStat = [];
    end


    % Extract response params 
    buffer = ceil(mean(diff(datarun.stimulus.triggers))); 
    bin_size = 0.001; % for binning spikes
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datarun,datarun.cell_ids, datarun.triggers(end)+buffer, bin_size) ;

    % If multiple contrasts 
    F = find(strcmp(fieldnames(datarun.stimulus.combinations),'rgb')); 
    G = find(strcmp(fieldnames(datarun.stimulus.combinations),'back_rgb')); 
    contrast_lev = unique(StimComb(:,F));
    contrast_indices = cell(length(contrast_lev),1);
    for i=1:length(contrast_lev)
        [contrast_indices{i},~,] = find(StimComb(:,F)==contrast_lev(i)); 
    end
    back_lev = unique(StimComb(:,G)); 
    contrast_prcnt = contrast_lev.*100./back_lev;

    retstat.stim_dur = buffer; 
    retstat.dirs = datarun.stimulus.params.direction;
    retstat.nrepeats = datarun.stimulus.repetitions;
    retstat.contrast_lev = contrast_lev;
    retstat.contrast_prcnt = contrast_prcnt;

    
    %% Find ds cells by parametrizing direction selective responses 

    % If you want to exclude any specific trials  
    trial_nums = 1:datarun.stimulus.repetitions; % default 

    % Classify DSGCs
    [ds_struct, ds_id, dsi_redef, rasterall, psth, binranges, hf1] = ds_classification_mb(datarun, oldStat, StimComb, buffer, 'trials', trial_nums);
    [~,ds_ind] = intersect(datarun.cell_ids, ds_id);
    ds_ind = ds_ind'; 

    % Add field names values 
    retstat.ds_struct_all = ds_struct; 
    retstat.dsi_redef = dsi_redef;
    retstat.ds_id = ds_id;
    retstat.ds_ind = ds_ind;
    retstat.non_ds_id = setxor(datarun.cell_ids,ds_id);
    [~,retstat.non_ds_ind] = get_cell_indices(datarun, retstat.non_ds_id); 
    retstat.rasterall = rasterall;
    retstat.psth_all = psth;
    retstat.binranges = binranges;
    retstat.trial_nums = trial_nums;


    %% Extract structures for DS cells and plot direction selectivity

    % If you want to choose psth hillocks for DSI calculation 
    hills = true; 

    if ~hills
        ds_struct_subset = get_ds_response_mb_native(datarun, StimComb, ds_ind, rasterall, trial_nums); % all spikes 
    else 
        ds_struct_subset = get_ds_response_mb_modf(datarun, StimComb, ds_ind, ds_struct, rasterall, psth, binranges, trial_nums); % band separated spikes
    end
    retstat.ds_struct_ds = ds_struct_subset;


    %% Run analysis on MB data using all DS cells: select and save ON, ON-OFF DS cell ids for posterity :)

    % Separate ON-OFF DS cells from ON DS cells 

    % Use DS cell ids for segregating ON and ON-OFF responses
    DSset = ds_id;

    % Check existence of oo_ds and o_ds cell ids 
    if (isfield(oldStat,'oo_ds_Ids') && isfield(oldStat,'o_ds_Ids'))~=0  
        commandwindow;
        inp = input('ON-DS and ON-OFF DS cell ids already exist. \n Recreate them? y/n \n','s');
    else
        inp = input('ON-DS and ON-OFF DS cell ids do not exist. \n Create them afresh? (y/n) \n','s');
    end


    if strcmp(inp,'n')==1
        fprintf(['*---------------------------------------------*\n Loading o_ds and oo_ds cell ids. \n',...
            '*---------------------------------------------*\n']);
        if isfield(oldStat,'oo_ds_Ids') 
            retstat.oo_ds_Ids = oldStat.oo_ds_Ids;
            retstat.oo_ds_Inds = oldStat.oo_ds_Inds; 

            oo_ds_Ids = oldStat.oo_ds_Ids;  
        else
            fprintf('oo_ds cell ids do not exist. Skipping creating them....\n');
        end
        if isfield(oldStat,'o_ds_Ids') 
            retstat.o_ds_Ids = oldStat.o_ds_Ids;
            retstat.o_ds_Inds = oldStat.o_ds_Inds; 

            o_ds_Ids = oldStat.o_ds_Ids;
        else
            fprintf('o_ds cell ids do not exist. Skipping creating them....\n');
        end

    else
        % Visually inspect rasters to separate ON-OFF from ON DS cell
        oo_ds_Ids = [];
        o_ds_Ids = []; 

        for rgc = 1:length(DSset)        

            mean_dsi = mean([retstat.ds_struct_ds.mag{1,1,end}(rgc) retstat.ds_struct_ds.mag{1,1,end-1}(rgc)]); 

            plot_raster_aggregate_mb(datarun, StimComb, DSset(rgc), mean_dsi);
            commandwindow;
            UserInStr = input([num2str(rgc),'/',num2str(length(DSset)),': Choose this as oo-DS cell? (y/n)'],'s') ;
            switch UserInStr 
                case 'y'
                    oo_ds_Ids = [oo_ds_Ids DSset(rgc)]; 
                case 'n' 
                    UserInStr2 = input('Choose this as o-DS cell? (y/n)','s') ;
                    if strcmp(UserInStr2,'y') 
                        o_ds_Ids = [o_ds_Ids DSset(rgc)]; 
                    end
            end
        end

        retstat.oo_ds_Ids = oo_ds_Ids; 
        retstat.o_ds_Ids = o_ds_Ids;

        [~,retstat.oo_ds_Inds] = intersect(datarun.cell_ids, oo_ds_Ids); 
        if ~isempty(o_ds_Ids) 
            [~,retstat.o_ds_Inds] = intersect(datarun.cell_ids, o_ds_Ids); 
        else
            retstat.o_ds_Inds = []; 
        end

        warning('Saving new o_ds and oo_ds cell ids. Will overwrite previous files!');  
    end



    %% Extract ds cells and return raster, psth and structure containing response vector direction and magnitudes for DSGCs

    % Classifies DS cells based on either 
    %    (1) user provided threshold "thresh" for summed vector magnitude, OR,
    %    (2) 'ginput' where user selects DS cell clusters using mouse 
    %
    % Input arguments: 
    %               datarun: datarun structure 
    %               oldStat: saved stats 
    %               StimComb: stimulus combinations 
    %               stim_dur: stimulus duration 
    %               optional arguments: supply trial indices as a vector (ex. 1:6)
    % 
    % Output arguments: 
    %               ds_struct: structure containing direction dependent
    %               responses
    %               ds_id: cell id of DSGCs
    %               dsi_redef: if dsi needs to be redefined 
    %               rasterall: cell array of raster of individual cells 
    %               psth: cell array of psth of individual cells 
    %               hf: figure handle 

    function [ds_struct, ds_id, dsi_redef, rasterall, psth, binranges, hf] = ds_classification_mb(datarun, oldStat, StimComb, stim_dur, varargin)

        p = inputParser; 
        addParameter(p, 'trials', []); 
        parse(p,varargin{:}); 
        params = p.Results;

        rasterall = cell(length(datarun.cell_ids),1);
        psth = cell(length(datarun.cell_ids),1);
        trial_nums = params.trials; % choice for all trials or some trials 
        cell_inds = get_cell_indices(datarun, datarun.cell_ids); % get all cell indices 

        commandwindow;
    %     inp = input('Options for classifying DS cells: \n (1) Two highest contrasts \n (2) All contrasts \nEnter 1, 2  \n','s');  
    %     inp = str2double(inp); 
        inp=1; 

        switch inp 

            case 1

                %---------------------------------------------------------------------
                bin = 0.1; %sec
                for k=1:length(datarun.cell_ids)
                    [rasterall{k},psth{k},binranges] = get_raster_aggregate(datarun, StimComb, datarun.cell_ids(k), bin); 
                end
                ds_struct = get_ds_response_mb_native(datarun, StimComb, cell_inds, rasterall, stim_dur, trial_nums); % all spikes 
                dsindex = ds_struct.mag; % DSI
                %---------------------------------------------------------------------

                % Use high constrast responses to classify DSGCs
                if isequal(lower(fieldnames(datarun.stimulus.combinations)), fieldnames(datarun.stimulus.combinations))
                    rgb = datarun.stimulus.params.rgb;
                elseif isequal(upper(fieldnames(datarun.stimulus.combinations)), fieldnames(datarun.stimulus.combinations))
                    rgb = datarun.stimulus.params.RGB;
                end
                tempmat1 = zeros(size(dsindex{1,1,1})); 
                tempmat2 = zeros(size(dsindex{1,1,1})); 
                for i=1:size(dsindex,1)
                    for j=1:size(dsindex,2)
                        if size(dsindex,3)>1                        
                            if iscell(rgb)
                                for k=1:length(rgb)
                                    rgbmat(k)= unique(rgb{k}); 
                                end
                                rgb = rgbmat; 
                            end
                            [B,I] = sort(rgb,'ascend');
                            id1 = I(end); 
                            id2 = I(end-1); 
                            tempmat1 = tempmat1 + dsindex{i,j,id1}; 
                            tempmat2 = tempmat2 + dsindex{i,j,id2}; 
                        else
                            tempmat1 = tempmat1 + dsindex{i,j,1}; 
                            tempmat2 = tempmat2 + dsindex{i,j,1}; 
                        end
                    end
                end 
                dsindex1 = tempmat1./(size(dsindex,1)*size(dsindex,2));
                dsindex2 = tempmat2./(size(dsindex,1)*size(dsindex,2));
                dsi_redef = [dsindex1 dsindex2];

            case 2

                %---------------------------------------------------------------------
                ds_struct = get_ds_response_mb_native(datarun, StimComb, cell_inds, rasterall, stim_dur, trial_nums);
                dsindex = ds_struct.mag; 
                %--------------------------------------------------------------------- 

                % Use all contrast conditions for DS classification 
                tempmat = zeros(size(dsindex{1,1,1})); 
                for c=1:size(dsindex,3)
                    tempmat = tempmat + dsindex{1,1,c}; 
                end
                dsindex1 = tempmat./size(dsindex,3); 
                dsindex2 = dsindex1; 
                dsi_redef = repmat((dsindex1+dsindex2)./2,1,2);     
        end


        % Generate figures
        hf=figure; clf ; 
        plot( log10(dsi_redef(:,1)), log10(dsi_redef(:,2)), '*'); 
        xlabel('log10(TP 1)')
        ylabel('log10(TP 2)')
        xl = min([log10(dsi_redef(:,1)); log10(dsi_redef(:,2))]); 
        set(gca,'xlim',[xl*1.2 0],'ylim',[xl*1.2 0]);
        axis square;
        hold on

        if isfield(oldStat,'ds_id')==1
            ds_id = oldStat.ds_id; 
            [~,tempindx] = intersect(datarun.cell_ids,ds_id); 
            hp = plot( log10(dsi_redef(tempindx,1)), log10(dsi_redef(tempindx,2)), 'ro'); 
            inp = input('Accept selection? y/n ' ,'s'); 

            if strcmp(inp,'n') || strcmp(inp,'')
                delete(hp); 
                [x, y] = ginput;
                plot(x, y);
                IN = inpolygon( log10(dsi_redef(:,1)),  log10(dsi_redef(:,2)), x, y);
                I = find(IN == 1);
                I(I>length(datarun.cell_ids))=[]; % catch mishaps
                plot( log10(dsi_redef(I,1)),  log10(dsi_redef(I,2)), 'ro');
                ds_id = datarun.cell_ids(I);

                inp = input('DS cell ids already exist. Overwrite them? y/n \n','s');
                if strcmp(inp,'y')~=1
                    fprintf('Not saving current selection \n'); 
                    ds_id = oldStat.ds_id;  % revert ds_id to the previously saved values
                end
            end
        else 
            [x, y] = ginput;
            plot(x, y);
            IN = inpolygon( log10(dsi_redef(:,1)),  log10(dsi_redef(:,2)), x, y);
            I = find(IN == 1);
            I(I>length(datarun.cell_ids))=[]; % catch mishaps
            plot( log10(dsi_redef(I,1)),  log10(dsi_redef(I,2)), 'ro');
            ds_id = datarun.cell_ids(I);
            fprintf('DS cell ids do not exist. Creating them. \n'); 
        end
        set(gca,'fontsize',18);
        fprintf('# of DS cells: %d \n',length(ds_id));


        [~,tempindx] = intersect(datarun.cell_ids,ds_id);  
        dsi_redef = dsi_redef(tempindx,:);
    end


    %% Separate ooDSGCs into subtypes 

    if isfield(retstat,'oo_ds_Ids') || exist('oo_ds_Ids','var') 

        [~,oo_ds_Inds] = intersect(datarun.cell_ids, oo_ds_Ids); 

        % Method : hierarchical clustering based on centroid of clusters (k-medoids Clustering) 
        cind = length(ds_struct.U)-1; % next to highest contrast 

        X = [ds_struct.U{1,1,cind}(oo_ds_Inds),ds_struct.V{1,1,cind}(oo_ds_Inds)]; % coordinate location of vector tip (highest contrast)
        for i=1:length(X)
            X(i,:) = X(i,:)./norm(X(i,:)); % normalize each vector
        end
        [clust_indx,C] = kmedoids(X, 4, 'Distance', 'Euclidean');  % assume 4 cardinal directions
        median_ang = atan2d(C(:,2), C(:,1)); % directions of median vectors
        D = [C zeros(4,1)];

        % Check if there are only three major directions, then collapse the
        % irrelevant direction (set threshold for minimum separation to 45deg)
        median_ang_360 = wrapTo360(median_ang); 
        diff_ = bsxfun(@minus, repmat(median_ang_360,1,1), median_ang_360'); 

        [r,c] = find((diff_<30 & diff_>0) | (diff_>330 & diff_<=360)); 
        if ~isempty(r) || ~isempty(c) 
            tempvar = mod((median_ang_360(r)+median_ang_360(c))/2,360);

            median_ang_360(r) = tempvar; 
            median_ang_360(c) = tempvar; clear tempvar;
            median_ang_360 = unique(median_ang_360); 

            clust_indx(union(find(clust_indx==r), find(clust_indx==c))) = min([r c]);
        end
        median_ang = wrapTo180(median_ang_360); % median directions 

        % Separate the oo-DS cells into subtypes 
        uidx = unique(clust_indx); 
        for i=1:length(uidx)
            if median_ang(i)>=45 && median_ang(i)<135 
                oo_ds_Inferior = oo_ds_Ids(clust_indx==uidx(i));  
            elseif median_ang(i)>=-135 && median_ang(i)<-45
                oo_ds_Superior = oo_ds_Ids(clust_indx==uidx(i));  
            elseif median_ang(i)>=-45 && median_ang(i)<45
                oo_ds_Posterior = oo_ds_Ids(clust_indx==uidx(i));  
            elseif ((median_ang(i)>=135 && median_ang(i)<180) || (median_ang(i)>=-180 && median_ang(i)<-135))
                oo_ds_Anterior = oo_ds_Ids(clust_indx==uidx(i));  
            end
        end

        cl = {'r','g','b','m'};
        figure; 
        for i=1:length(oo_ds_Inds)

            subplot(1,2,1); 
            h = compass(ds_struct.U{1,1,cind}(oo_ds_Inds(i)),ds_struct.V{1,1,cind}(oo_ds_Inds(i))); hold on
            a_x = get(h, 'xdata'); 
            b_y = get(h, 'ydata'); 
            set(h, 'xdata', a_x(1:2), 'ydata', b_y(1:2), 'color', [0 0.4470 0.7410]);

            subplot(1,2,2); 
            h = compass(ds_struct.U{1,1,cind}(oo_ds_Inds(i)),ds_struct.V{1,1,cind}(oo_ds_Inds(i))); hold on;
            a_x = get(h, 'xdata'); 
            b_y = get(h, 'ydata'); 
            set(h, 'xdata', a_x(1:2), 'ydata', b_y(1:2), 'color', cl{clust_indx(i)}); 
        end


        if exist('oo_ds_Inferior','var'); retstat.oo_ds_Inferior = oo_ds_Inferior; end 
        if exist('oo_ds_Superior','var'); retstat.oo_ds_Superior = oo_ds_Superior; end
        if exist('oo_ds_Posterior','var'); retstat.oo_ds_Posterior = oo_ds_Posterior; end
        if exist('oo_ds_Anterior','var'); retstat.oo_ds_Anterior = oo_ds_Anterior; end

    else 

        fprintf('Skipping clustering ON-OFF DS cells ....\n');
    end


    %% Return tuning stats 
    retstat = get_ooDS_phaseResp(retstat, hills);

end









%% Compute stimulus combinations and spikes elicited for those combinations      
%
% Input arguments: 
%               datarun: datarun structure
%               timedur: duration after trigger to count spikes        
%               cellids: analyzed cells     
%               binsize: size within which to bin the data
%
% Output arguments: 
%               NumSpikesCell: average total number of spikes for each stimulus type for each cell
%               StimComb: all the stimulus combinations (col 1: spatial,col 2: temporal, col 3: direction)
%
% Copyright: Sneha Ravi, 2012; Greg Field, 2015; Jon Cafaro, 2015; Suva
% Roy, 2016 

function [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datarun,cellids,timedur, bin_size)

    % Initialize output matrices          
    NumSpikesCell = zeros(length(cellids), length(datarun.stimulus.combinations));
    MaxRate = zeros(length(cellids), length(datarun.stimulus.combinations));
    StimComb = zeros(length(datarun.stimulus.combinations), length(fieldnames(datarun.stimulus.params)));

    % Convert upper to lower characters
    s = fieldnames(datarun.stimulus.combinations);
    if isequal(upper(s),s)
        for sc = 1:length(s)
            s(sc) = lower(s(sc)); 
        end
    end

    % Allocate StimComb matrix to memory 
    for i=1:length(s)
        StimComb(:,i) = cell2mat({datarun.stimulus.combinations(:).(s{i})});   
    end

    %Calculation of total spike number for each trial for each cell, then average calculated and placed in NumSpikesCell
    temp_indices = get_cell_indices(datarun, cellids);

    % rgc loops over cells
    for rgc = 1:length(temp_indices)
        NumSpikesAll = zeros(1, length(datarun.stimulus.trials));
        MaxRateAll = zeros(1, length(datarun.stimulus.trials));
        % loop over epochs
        for epoch = 1:length(datarun.stimulus.triggers) 
            % get spike indices between triggers
            if epoch~=length(datarun.stimulus.triggers) % if its not the last trigger
                Spikes_idx = datarun.spikes{rgc} >= datarun.stimulus.triggers(epoch) & datarun.spikes{rgc} < datarun.stimulus.triggers(epoch+1);
                epoch_duration = datarun.stimulus.triggers(epoch+1) - datarun.stimulus.triggers(epoch);
            else
                Spikes_idx = datarun.spikes{rgc} >= datarun.stimulus.triggers(epoch) & datarun.spikes{rgc} < timedur ;
                epoch_duration = timedur - datarun.stimulus.triggers(epoch) ;
            end
            % get spike times
            Spikes_temp = datarun.spikes{rgc}(Spikes_idx);
            % subtract off time of first trigger so that start is zero 
            Spikes_temp = Spikes_temp - datarun.stimulus.triggers(epoch);
            % compute spike number
            NumSpikesAll(epoch) = sum(Spikes_idx);
            MaxRateAll(epoch) = max(histc(Spikes_temp, 0:bin_size:epoch_duration));
        end
        size(NumSpikesAll) ;
        size(MaxRateAll) ;

        NumSpikesCell(rgc,:) = grpstats(NumSpikesAll, datarun.stimulus.trial_list, 'sum'); 
        MaxRate(rgc,:) = grpstats(NumSpikesAll, datarun.stimulus.trial_list, 'max'); 
    end 
end


%% Extract DSGCs and return raster, psth and response vector for each stimulus type 
% 
% Input arguments: 
%               datarun: datarun structure
%               oldStat: saved stats if you desire to use it 
%               StimComb: stimulus combinations 
%               stim_dur: stimulus durations 
%               optional arguments: trials 
%
% Output arguments: 
%               ds_struct: structure containing direction selective
%               responses
%               ds_id: ids of DSGCs
%               dsi_redef: redefined dsi if saved stat is not used 
%               rasterall: cell array of raster of each RGC
%               psth: cell array of psth of each RGC
%               binranges: bins for plotting 
%               hf: figure handle 
%
% Copyright: Jon Cafaro, 2015; Suva Roy, 2016 


function [ds_struct, ds_id, dsi_redef, rasterall, psth, binranges, hf] = ds_classification_mb(datarun, oldStat, StimComb, stim_dur, varargin)

    
% Classifies DS cells based on either 
%    (1) user provided threshold "thresh" for summed vector magnitude, OR,
%    (2) 'ginput' where user selects DS cell clusters using mouse 


    p = inputParser; 
    addParameter(p, 'trials', []); 
    parse(p,varargin{:}); 
    params = p.Results;
       
    rasterall = cell(length(datarun.cell_ids),1);
    psth = cell(length(datarun.cell_ids),1);
    trial_nums = params.trials; % choice for all trials or some trials 
    cell_inds = datarun.cell_inds; % get all cell indices 
    
    inp=1; 
    
    switch inp 
        
        case 1
            
            %---------------------------------------------------------------------
            bin = 0.1; %sec
            for k=1:length(datarun.cell_ids)
                [rasterall{k},psth{k},binranges] = get_raster_aggregate(datarun, StimComb, datarun.cell_ids(k), bin); 
            end
            ds_struct = get_ds_response_mb_native(datarun, StimComb, cell_inds, rasterall, stim_dur, trial_nums); % all spikes 
            dsindex = ds_struct.mag; % DSI
            
            % Use high contrast responses for DSGC classification
            if isequal(lower(fieldnames(datarun.stimulus.combinations)), fieldnames(datarun.stimulus.combinations))
                rgb = datarun.stimulus.params.rgb;
            elseif isequal(upper(fieldnames(datarun.stimulus.combinations)), fieldnames(datarun.stimulus.combinations))
                rgb = datarun.stimulus.params.RGB;
            end
            tempmat1 = zeros(size(dsindex{1,1,1})); 
            tempmat2 = zeros(size(dsindex{1,1,1})); 
            for i=1:size(dsindex,1)
                for j=1:size(dsindex,2)
                    if size(dsindex,3)>1                        
                        if iscell(rgb)
                            for k=1:length(rgb)
                                rgbmat(k)= unique(rgb{k}); 
                            end
                            rgb = rgbmat; 
                        end
                        [B,I] = sort(rgb,'ascend');
                        id1 = I(end); % highest contrast
                        id2 = I(end-1); % 2nd highest contrast
                        tempmat1 = tempmat1 + dsindex{i,j,id1}; 
                        tempmat2 = tempmat2 + dsindex{i,j,id2}; 
                    else
                        tempmat1 = tempmat1 + dsindex{i,j,1}; 
                        tempmat2 = tempmat2 + dsindex{i,j,1}; 
                    end
                end
            end 
            dsindex1 = tempmat1./(size(dsindex,1)*size(dsindex,2));
            dsindex2 = tempmat2./(size(dsindex,1)*size(dsindex,2));
            dsi_redef = [dsindex1 dsindex2];
            
        case 2
            
            ds_struct = get_ds_response_mb_v2(datarun, StimComb, cell_inds, rasterall, stim_dur, trial_nums);
            dsindex = ds_struct.mag;
            
            % Use all contrast conditions for DS classification 
            tempmat = zeros(size(dsindex{1,1,1})); 
            for c=1:size(dsindex,3)
                tempmat = tempmat + dsindex{1,1,c}; 
            end
            dsindex1 = tempmat./size(dsindex,3); 
            dsindex2 = dsindex1; 
            dsi_redef = repmat((dsindex1+dsindex2)./2,1,2);     
            
    end
            
    
    % Generate figures
    hf=figure; clf ; 
    plot( log10(dsi_redef(:,1)), log10(dsi_redef(:,2)), '*'); 
    xlabel('log10(TP 1)')
    ylabel('log10(TP 2)')
    xl = min([log10(dsi_redef(:,1)); log10(dsi_redef(:,2))]); 
    set(gca,'xlim',[xl*1.2 0],'ylim',[xl*1.2 0]);
    axis square;
    hold on

    if isfield(oldStat,'ds_id')==1
        ds_id = oldStat.ds_id; 
        [~,tempindx] = intersect(datarun.cell_ids,ds_id); 
        hp = plot( log10(dsi_redef(tempindx,1)), log10(dsi_redef(tempindx,2)), 'ro'); 
        inp = input('Accept selection? y/n ' ,'s'); 

        if strcmp(inp,'n') || strcmp(inp,'')
            delete(hp); 
            [x, y] = ginput;
            plot(x, y);
            IN = inpolygon( log10(dsi_redef(:,1)),  log10(dsi_redef(:,2)), x, y);
            I = find(IN == 1);
            I(I>length(datarun.cell_ids))=[]; % catch mishaps
            plot( log10(dsi_redef(I,1)),  log10(dsi_redef(I,2)), 'ro');
            ds_id = datarun.cell_ids(I);

            inp = input('DS cell ids already exist. Overwrite them? y/n \n','s');
            if strcmp(inp,'y')~=1
                fprintf('Not saving current selection \n'); 
                ds_id = oldStat.ds_id;  % revert ds_id to the previously saved values
            end
        end
    else 
        [x, y] = ginput;
        plot(x, y);
        IN = inpolygon( log10(dsi_redef(:,1)),  log10(dsi_redef(:,2)), x, y);
        I = find(IN == 1);
        I(I>length(datarun.cell_ids))=[]; % catch mishaps
        plot( log10(dsi_redef(I,1)),  log10(dsi_redef(I,2)), 'ro');
        ds_id = datarun.cell_ids(I);
        fprintf('DS cell ids do not exist. Creating them. \n'); 
    end
    set(gca,'fontsize',18);
    fprintf('# of DS cells: %d \n',length(ds_id));

    [~,tempindx] = intersect(datarun.cell_ids,ds_id);  
    dsi_redef = dsi_redef(tempindx,:);    
end


%% Get raster, psth and bins for individual cell
% 
% Output arguments: 
%               spikecell: raster cell array for different contrasts and directions 
%               psth: psth cell array for different contrasts and directions 
%               binranges: time bins 
% 
% Input arguments: 
%               datarun: datarun structure 
%               StimComb: stimulus combinations 
%               cell_id: cell id 
%               bin: width of time bin
% 
%
% Copyright: Sneha Ravi, 2012; Greg Field, 2015; Jon Cafaro, 2015; Suva
% Roy, 2016 


function [spikecell, psth, binranges] = get_raster_aggregate(datarun, StimComb, cell_id, bin)


    s = fieldnames(datarun.stimulus.combinations);
    if isequal(lower(s),s)
        F = find(strcmp(fieldnames(datarun.stimulus.combinations),'rgb'));
    elseif isequal(upper(s),s)
        F = find(strcmp(fieldnames(datarun.stimulus.combinations),'RGB'));
    end
  
    contrast_lev = unique(StimComb(:,F));
    contrast_indices = cell(length(contrast_lev),1);
    for i=1:length(contrast_lev)
        [contrast_indices{i},~,] = find(StimComb(:,F)==contrast_lev(i)); 
    end

    D = find(strcmp(fieldnames(datarun.stimulus.combinations),'direction'));
    dirs = unique(StimComb(:,D)); 
     
    
    temp_index = find(datarun.cell_ids==cell_id); % cell ID
    stim_dur = mean(diff(datarun.stimulus.triggers));

    binranges = 0:bin:stim_dur; % choose appropriate bin size
    ncounts = zeros(datarun.stimulus.repetitions,length(binranges));

    numCNTR = length(contrast_lev); 
    numDR = length(dirs);
    [spikecell,psth] = deal(cell(numCNTR,numDR));


    % figure(20); hold on; 
    % cnt = 1;
    for cn = 1:numCNTR
        nspikes = zeros(1,numDR); 
        dirmat(cn,:) = StimComb(contrast_indices{cn},D)'; 
        for dr = 1:numDR
            idx_first_rep = contrast_indices{cn}(dr);
            [~, idx_trial_list] = find(datarun.stimulus.trial_list==idx_first_rep);
            temp_triggers = datarun.stimulus.triggers(idx_trial_list);                
            temp_epochs = get_raster(datarun.spikes{temp_index}, temp_triggers,'stop',stim_dur,'plot',false);     
            spikecell{cn,dr} = temp_epochs; 

            % PSTH 
            for te=1:length(temp_epochs)
                ncounts(te,:) = histc(temp_epochs{te},binranges); 
            end 
            ncounts = sum(ncounts,1); % sum over trials
            if sum(ncounts)~=0
                psth{cn,dr} = ncounts;  % normalized PSTH 
            else
                psth{cn,dr} = zeros(1,length(binranges)); 
            end
            
            nspikes(dr) = sum(ncounts); 
        end
        
        % Normalize by sum of spikes collected over all directions 
        nspikes_alldir = sum(nspikes); 
        if nspikes_alldir~=0
            for dr=1:numDR
                psth{cn,dr} = psth{cn,dr}./nspikes_alldir;
            end
        end
    end

    % Order psth and rasterall by directions of movement 
    [tempcell1,tempcell2] = deal(cell(numCNTR,numDR));
    for cn = 1:numCNTR
        [~,I] = sort(dirmat(cn,:));
        for dr = 1:numDR
            tempcell1{cn,dr} = spikecell{cn,I(dr)};
            tempcell2{cn,dr} = psth{cn,I(dr)}; 
        end
    end

    spikecell = tempcell1;
    psth = tempcell2;

end



%% Compute ds_struct redefined if needed 
%
% Calculate average total number of spikes for each stimulus type      
%
% Input arguments: 
%               datarun: datarun structure
%               StimComb: stimulus combinations 
%               cell_inds: cell indices
%               rasterall: cell array of rasters of all RGCs   
%               stim_dur: stim duration 
%               trial_nums: trial set to use 
%
% Output arguments: 
%               ds_struct: structure containing direction selective
%               responses
%
% Copyright: Suva Roy, 2016


function [ds_struct_redef] = get_ds_response_mb_native(datarun, StimComb, cell_inds, rasterall, stim_dur, trial_nums)


    BW = find(strcmp(fieldnames(datarun.stimulus.combinations),'bar_width'));
    barwidth = unique(StimComb(:,BW));
    BS = find(strcmp(fieldnames(datarun.stimulus.combinations),'delta')); 
    barspeed = unique(StimComb(:,BS));  
    D = find(strcmp(fieldnames(datarun.stimulus.combinations),'direction')); 
    dirs = unique(StimComb(:,D));  
    F = find(strcmp(fieldnames(datarun.stimulus.combinations),'rgb')); 
    contrast_lev = unique(StimComb(:,F));

    [mag,dsindex,angle,spkcnt,rho,theta,U,V] = deal(cell(length(barwidth),length(barspeed),length(contrast_lev)));       
    
    %----------------------------------------------------------------------
    % Here's a temporary implementation of modified DSI - multiply DSI with
    % a "regularity index" which defines trial-to-trial regularity in
    % spiking
    warning('This is a temporary implementation of modified DSI');
    bin = 0.001; % 1 ms bins
    trial_bins=0:bin:stim_dur;
    C = nchoosek(trial_nums,2); % all trial pairs

    
    for i = 1:length(barwidth)
        for j = 1:length(barspeed)
            for nc = 1:length(cell_inds)
                for k = 1:length(contrast_lev)
                    clear r; 
                    for dr = 1:length(dirs)                        
                         tempvar = cellfun(@(x) length(x(:,1)), rasterall{cell_inds(nc)}{k,dr}, 'UniformOutput' ,0); 
                         spikecnt_by_trial = cell2mat(tempvar); 
                         r(dr) = sum(spikecnt_by_trial(trial_nums)); % total number of spikes along each direction 
                         
                    end
                    r = r'; 
                    
                    norm_ = sum(r);  
                    r_norm = r./norm_; 
                    r_norm(isnan(r_norm)) = 0; 
                    t = dirs.*pi./180;
                    
                    [X Y] = pol2cart(t, r_norm); 
                    
                    u = sum(X);
                    v = sum(Y);
                    [an,ma] = cart2pol(u, v); 
                    
                    rho{i,j,k}(nc,:) = r_norm;
                    theta{i,j,k}(nc,:) = t;
                    U{i,j,k}(nc,1) = u;
                    V{i,j,k}(nc,1) = v;
                    angle{i,j,k}(nc,1) = an;
                    mag{i,j,k}(nc,1) = ma;
                    
                end
            end
        end
    end
           
    ds_struct_redef = v2struct(mag, dsindex, angle, rho, spkcnt, theta, U, V);
end



%% Compute ds_struct using psth bands to minimize bias in summed vector response by motion independent spikes 
%
% Calculate average total number of spikes for each stimulus type      
%
% Input arguments: 
%               datarun: datarun structure
%               StimComb: stimulus combinations 
%               cell_inds: cell indices
%               ds_struct: needs previously calculated 'ds_struct' 
%               rasterall: cell array of rasters of all RGCs   
%               psth: cell array of psth of all RGCs
%               stim_dur: stim duration 
%               trial_nums: trial set to use 
%
% Output arguments: 
%               ds_struct: structure containing direction selective
%               responses
%
% Copyright: Suva Roy, 2016


function [ds_struct_redef] = get_ds_response_mb_modf(datarun, StimComb, cell_inds, ds_struct, rasterall, psth, binranges, trial_nums)

    

    BW = find(strcmp(fieldnames(datarun.stimulus.combinations),'bar_width'));
    barwidth = unique(StimComb(:,BW));
    BS = find(strcmp(fieldnames(datarun.stimulus.combinations),'delta')); 
    barspeed = unique(StimComb(:,BS));  
    D = find(strcmp(fieldnames(datarun.stimulus.combinations),'direction')); 
    dirs = unique(StimComb(:,D));  
    F = find(strcmp(fieldnames(datarun.stimulus.combinations),'rgb')); 
    contrast_lev = unique(StimComb(:,F));
    

    [mag,dsindex,angle,spkcnt,rho,theta,U,V,on_spikes,off_spikes,all_spikes] = deal(cell(length(barwidth),length(barspeed),length(contrast_lev)));       
    
    
    % Set up bands for limiting spike accumulation region in a trial
    band_size = 1.5; % sec (is bar width and speed specific, modify accordingly) 
    g= gausswin(8); % Gaussian filter  (parameter value: 8 or 9 works best under the given stimulus constraints) 
    g = g./sum(g); 
    
    for i = 1:length(barwidth)
        for j = 1:length(barspeed)
            for nc = 1:length(cell_inds)
                
                %----------------------------------------------------------
                % Determine inclusion bands in PSTH for spike count 
                
                % STEP 1 
                % Obtain preferred direction (PD) of the cell at high contrasts 
                atContrastIndx = length(contrast_lev)-1;
                PD = ds_struct.angle{1,1,atContrastIndx}(cell_inds(nc)) * 180/pi;
                PD = wrapTo360(PD);
                [~,k_indx] = min(abs(dirs-PD)); 
                closestDir = dirs(k_indx); % movement direction closest to PD (degree)
                
                % Determine inclusion band using PSTH @ highest contrast for the preferred direction 
                if sum((psth{cell_inds(nc)}{atContrastIndx,k_indx}))~=0
                    psth_smooth = conv(psth{cell_inds(nc)}{atContrastIndx,k_indx},g,'same'); % choose PSTH corresponding to highest contrast 
                    [pks, loc] = findpeaks(psth_smooth,binranges,'npeaks',2,'SortStr','descend'); % peak separation >1sec
                    [loc, Ind] = sort(loc); 
                    pks = pks(Ind); 
                    switch length(loc)
                        case 0
                            band{k_indx}(1,:) = [(binranges(end)/2)-(band_size/2) (binranges(end)/2)+(band_size/2)]; 
                            band{k_indx}(2,:) = band{k_indx}(1,:); 
                        case 1
                            band{k_indx}(1,:) = [loc(1)-band_size/2 loc(1)+band_size/2];
                            band{k_indx}(2,:) = band{k_indx}(1,:);  % use the same band 
                        case 2
                            band{k_indx}(1,:) = [loc(1)-band_size/2 loc(1)+band_size/2];
                            band{k_indx}(2,:) = [loc(2)-band_size/2 loc(2)+band_size/2];  
                    end
                    % remove outliers
                    band{k_indx}(band{k_indx}<0) = 0; 
                    band{k_indx}(band{k_indx}>binranges(end)) = binranges(end); 
                    % adjust band cross-over
                    if band{k_indx}(1,2)>band{k_indx}(2,1)
                        band{k_indx}(1,2) = (band{k_indx}(1,2)+band{k_indx}(2,1))./2;
                        band{k_indx}(2,1) = band{k_indx}(1,2); 
                    end
                else 
                    band{k_indx} = repmat([(binranges(end)/2)-(band_size/2) (binranges(end)/2)+(band_size/2)], 2, 1); 
                end
                
                
                % STEP 2
                % Assign bands to PSTHs along other directions preserving the band size   
                band_widths =  [diff(band{k_indx}(1,:)); diff(band{k_indx}(2,:))] ; 
                not_k_indx = setdiff(1:length(dirs), k_indx); % index of all other directions
                
                for dr = not_k_indx % For all other directions: Determine peaks and assign bands of same width
                    if sum((psth{cell_inds(nc)}{atContrastIndx,dr}))~=0
                        psth_smooth = conv(psth{cell_inds(nc)}{atContrastIndx,dr},g,'same'); % choose PSTH corresponding to highest contrast 
                        [pks, loc] = findpeaks(psth_smooth,binranges,'npeaks',2,'SortStr','descend'); % peak separation >1sec
                        [loc, Ind] = sort(loc); 
                        pks = pks(Ind); 
                        switch length(loc)
                            case 0
                                band{dr} = band{k_indx};
                            case 1
                                band{dr}(1,:) = [loc(1)-band_widths(1) loc(1)+band_widths(2)]; % supply band widths from master psth 
                                band{dr}(2,:) = band{dr}(1,:);  % use the same band 
                            case 2
                                if diff(loc)<mean(band_widths)
                                    band{dr}(1,:) = [mean(loc)-band_widths(1) mean(loc)];
                                    band{dr}(2,:) = [mean(loc) mean(loc)+band_widths(2)];
                                else
                                    band{dr}(1,:) = [loc(1)-band_widths(1)/2 loc(1)+band_widths(1)/2];
                                    band{dr}(2,:) = [loc(2)-band_widths(2)/2 loc(2)+band_widths(2)/2];  
                                end
                        end
                        % remove outliers
                        band{dr}(band{dr}<0) = 0; 
                        band{dr}(band{dr}>binranges(end)) = binranges(end); 
                        % adjust band cross-over
                        if band{dr}(1,2)>band{dr}(2,1)
                            band{dr}(1,2) = (band{dr}(1,2)+band{dr}(2,1))./2;
                            band{dr}(2,1) = band{dr}(1,2); 
                        end
                    else 
                        band{dr} = repmat([(binranges(end)/2)-(band_size/2) (binranges(end)/2)+(band_size/2)], 2, 1); 
                    end
                end
                %----------------------------------------------------------
                
                % Get spike counts within the bands 
                for k = 1:length(contrast_lev)
                    for dr = 1:length(dirs)      
                        spk_times = cell2mat(rasterall{cell_inds(nc)}{k,dr}(trial_nums));
                        
                        on_spikes{i,j,k}(nc,dr) = sum(spk_times>band{dr}(1,1) & spk_times<band{dr}(1,2));
                        off_spikes{i,j,k}(nc,dr) = sum(spk_times>=band{dr}(2,1) & spk_times<band{dr}(2,2));
                        all_spikes{i,j,k}(nc,dr) = on_spikes{i,j,k}(nc,dr) + off_spikes{i,j,k}(nc,dr); 
                        
                    end
                   
                    r = all_spikes{i,j,k}(nc,:)'; 
                    
                    norm_ = sum(r);  
                    r_norm = r./norm_; 
                    r_norm(isnan(r_norm)) = 0; 
                    t = dirs.*pi./180;
                    
                    [X Y] = pol2cart(t, r_norm); 
                    
                    u = sum(X);
                    v = sum(Y);
                    [an,ma] = cart2pol(u, v); 
                    
                    rho{i,j,k}(nc,:) = r_norm;
                    theta{i,j,k}(nc,:) = t;
                    U{i,j,k}(nc,1) = u;
                    V{i,j,k}(nc,1) = v;
                    angle{i,j,k}(nc,1) = an;
                    mag{i,j,k}(nc,1) = ma;
                    
                end
            end
        end
    end
           
    ds_struct_redef = v2struct(mag, dsindex, angle, rho, spkcnt, theta, U, V, on_spikes, off_spikes, all_spikes);
end



%% Return figure handle for raster plots (identifying ooDSGCs from oDSGCs)
%
% Input arguments: 
%               datarun: datarun structure
%               StimComb: stimulus combinations 
%               cell_id: cell ids
%               mean_dsi: dsi averaged over contrasts or speeds 
%
% Output arguments: 
%               handles: figure handle 
%
% Copyright: Suva Roy, 2016

function [handles] = plot_raster_aggregate_mb(datarun, StimComb, cell_id, mean_dsi)


    s = fieldnames(datarun.stimulus.combinations);
    if isequal(lower(s),s)
        F = find(strcmp(fieldnames(datarun.stimulus.combinations),'rgb'));
    elseif isequal(upper(s),s)
        F = find(strcmp(fieldnames(datarun.stimulus.combinations),'RGB'));
    end

    contrast_lev = unique(StimComb(:,F));
    contrast_indices = cell(length(contrast_lev),1);
    for i=1:length(contrast_lev)
        [contrast_indices{i},~,] = find(StimComb(:,F)==contrast_lev(i)); 
    end


    temp_index = find(datarun.cell_ids==cell_id); % cell ID
    stim_dur = mean(diff(datarun.stimulus.triggers));

    if exist('hf100','var')
        clear hf100;
    end

    hf100=figure(100); clf(gcf);
    set(hf100,'position',[73 65 1300 740]);
    suptitle(['Cell id: ',num2str(cell_id)]);

    R = length(contrast_lev); 
    C = length(contrast_indices{1});
    dirmat = zeros(R,C);
    handles = cell(R,C);
    spikecell = cell(R,C);

    for cl = 1:R
        dirmat(cl,:) = StimComb(contrast_indices{cl},3)'; % Note: 3rd col of StimComb has "directions" parameter  
        for dr = 1:C
            idx_first_rep = contrast_indices{cl}(dr);
            [~, idx_trial_list] = find(datarun.stimulus.trial_list==idx_first_rep);
            temp_triggers = datarun.stimulus.triggers(idx_trial_list);                
            temp_epochs = get_raster(datarun.spikes{temp_index}, temp_triggers,'stop',stim_dur,'plot',false);     
            spikecell{cl,dr} = temp_epochs; 
        end
    end

    % Order raster by directions of movement 
    cnt = 1; 
    for cl = 1:R
        [~,I] = sort(dirmat(cl,:));
        for dr = 1:C
            handles{cl,dr} = subplot(R,C,cnt);  
            
            raster = spikecell{cl,I(dr)}; 
            start_time = 0; 
            end_time = stim_dur; 
            for j = 1:length(raster)
                SpikeTime = raster{j};
                SpikeTime = SpikeTime';
                X = [SpikeTime; SpikeTime];
                Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
                line(handles{cl,dr}, X ,Y, 'color', 'b', 'linewidth', 1);
                hold on;
            end
            set(handles{cl,dr},'xlim',[start_time  end_time],'ylim',[0 length(raster)]);
            
            cnt = cnt+1;
        end
    end
    suptitle(['Cell id: ',num2str(cell_id),';  DSI: ',num2str(mean_dsi)]);
end



%% Return direction tuning stats 
% 
% Input arguments: 
%               retstat: uses predetermined structure for appending new
%               fields
%               hills: preset option for using ON and OFF phases only 
%
% Output arguments: 
%               retstat: returns fields containing tuning stats 
%
% Copyright: Suva Roy, 2016


function [retstat] = get_ooDS_phaseResp(retstat, hills)

    %---- Set up parameters used in calculation ----------------% 

    ds_id = retstat.ds_id;
    ds_ind = retstat.ds_ind;
    non_ds_id = retstat.non_ds_id;     
    non_ds_ind = retstat.non_ds_ind; 

    if isfield(retstat, 'o_ds_Ids') && isfield(retstat, 'oo_ds_Ids') 
        o_ds_Ids = retstat.o_ds_Ids;
        o_ds_Inds = retstat.o_ds_Inds;
        oo_ds_Ids = retstat.oo_ds_Ids;
        oo_ds_Inds = retstat.oo_ds_Inds;
        actual_ds_Ids = [o_ds_Ids oo_ds_Ids]; % These cells will be used to compute stats
        actual_ds_Inds = [o_ds_Inds oo_ds_Inds]; 

    else
        actual_ds_Ids = ds_id; % These cells will be used to compute stats
        actual_ds_Inds = ds_ind; 
    end

    rasterall = retstat.rasterall;
    psth = retstat.psth_all;
    binranges = retstat.binranges;
    contrast_lev = retstat.contrast_lev;
    ds_struct_ds = retstat.ds_struct_ds;
    ds_struct = retstat.ds_struct_all;

    stim_dur = retstat.stim_dur;
    dirs = retstat.dirs; 
    trial_nums = retstat.trial_nums; 

    R = length(retstat.contrast_lev); % contrasts
    C = length(dirs); % directions

    inp = 'y';
    if strcmp(inp, 'y')

        %---- Calculation for all DS and non-DS cells --------------------% 

        % Accumulate total number of spikes for all contrast conditions and
        % directions for non-DS cells (summed over nrepeats)
        for nc = 1:length(non_ds_ind)
            for cont = 1:R  % contrast 
                for dr = 1:C % directions
                    total_spikes_nonDS{nc}(cont,dr) = length(cell2mat(rasterall{non_ds_ind(nc)}{cont,dr})); 
                end
            end
        end

        % Compute tuning width and tuning strength of all DS cells: Apply
        % inclusion bands to PSTH for this calculation 
        dirs_cent = [dirs(2:end) 360]-180; % center around 0 deg (Will be used for centering tuning curves around 0 degree) 
        for nc = 1:length(actual_ds_Inds)

            if hills 
                %----------------------------------------------------------
                % Determine inclusion bands in PSTH for spike count 

                % Set up bands for limiting spike accumulation region in a trial
                band_size = 2.0; % sec 
                g= gausswin(8); % Gaussian filter  (parameter value: 8 or 9 works best) 
                g = g./sum(g);


                % STEP 1 
                % Obtain preferred direction (PD) of the cell at high contrast 
                atContrastIndx = length(contrast_lev)-1;
                PD = ds_struct.angle{1,1,atContrastIndx}(actual_ds_Inds(nc)) * 180/pi;
                PD = wrapTo360(PD);
                [~,k_indx] = min(abs(dirs-PD)); 
                closestDir = dirs(k_indx); % movement direction closest to PD (degree)

                % Determine inclusion band using PSTH @ highest contrast
                % for the preferred direction 
                if sum((psth{actual_ds_Inds(nc)}{atContrastIndx,k_indx}))~=0
                    psth_smooth = conv(psth{actual_ds_Inds(nc)}{atContrastIndx,k_indx},g,'same'); % choose PSTH corresponding to highest contrast 

                    psth_smooth_rec(k_indx,:) = psth_smooth; 

                    [pks, loc] = findpeaks(psth_smooth,binranges,'npeaks',2,'SortStr','descend'); % peak separation >1sec
                    [loc, Ind] = sort(loc); 
                    pks = pks(Ind); 
                    switch length(loc)
                        case 0 
                            band{k_indx}(1,:) = [(binranges(end)/2)-(band_size/2) (binranges(end)/2)+(band_size/2)]; 
                            band{k_indx}(2,:) = band{k_indx}(1,:); 
                        case 1
                            band{k_indx}(1,:) = [loc(1)-band_size/2 loc(1)+band_size/2];
                            band{k_indx}(2,:) = band{k_indx}(1,:);  % use the same band 
                        case 2
                            band{k_indx}(1,:) = [loc(1)-band_size/2 loc(1)+band_size/2];
                            band{k_indx}(2,:) = [loc(2)-band_size/2 loc(2)+band_size/2];  
                    end
                    % remove outliers
                    band{k_indx}(band{k_indx}<0) = 0; 
                    band{k_indx}(band{k_indx}>binranges(end)) = binranges(end); 
                    % adjust band cross-over
                    if band{k_indx}(1,2)>band{k_indx}(2,1)
                        band{k_indx}(1,2) = (band{k_indx}(1,2)+band{k_indx}(2,1))./2;
                        band{k_indx}(2,1) = band{k_indx}(1,2); 
                    end
                else 
                    band{k_indx} = repmat([(binranges(end)/2)-(band_size/2) (binranges(end)/2)+(band_size/2)], 2, 1); 
                end


                % STEP 2
                % Now assign bands to PSTHs along other directions preserving the band size   
                band_widths =  [diff(band{k_indx}(1,:)); diff(band{k_indx}(2,:))] ; 
                not_k_indx = setdiff(1:length(dirs), k_indx); % index of all other directions

                for dr = not_k_indx % For all other directions: Determine location of peaks and assign bands of same width
                    if sum((psth{actual_ds_Inds(nc)}{atContrastIndx,dr}))~=0
                        psth_smooth = conv(psth{actual_ds_Inds(nc)}{atContrastIndx,dr},g,'same'); % choose PSTH corresponding to highest contrast 
                        psth_smooth_rec(dr,:) = psth_smooth; 
                        [pks, loc] = findpeaks(psth_smooth,binranges,'npeaks',2,'SortStr','descend'); % peak separation >1sec
                        [loc, Ind] = sort(loc); 
                        pks = pks(Ind); 
                        switch length(loc)
                            case 0
                                band{dr} = band{k_indx};
                            case 1
                                band{dr}(1,:) = [loc(1)-band_widths(1) loc(1)+band_widths(2)]; % supply band widths from master psth 
                                band{dr}(2,:) = band{dr}(1,:);  % use the same band 
                            case 2
                                if diff(loc)<mean(band_widths)
                                    band{dr}(1,:) = [mean(loc)-band_widths(1) mean(loc)];
                                    band{dr}(2,:) = [mean(loc) mean(loc)+band_widths(2)];
                                else
                                    band{dr}(1,:) = [loc(1)-band_widths(1)/2 loc(1)+band_widths(1)/2];
                                    band{dr}(2,:) = [loc(2)-band_widths(2)/2 loc(2)+band_widths(2)/2];  
                                end
                        end
                        % remove outliers
                        band{dr}(band{dr}<0) = 0; 
                        band{dr}(band{dr}>binranges(end)) = binranges(end); 
                        % adjust band cross-over
                        if band{dr}(1,2)>band{dr}(2,1)
                            band{dr}(1,2) = (band{dr}(1,2)+band{dr}(2,1))./2;
                            band{dr}(2,1) = band{dr}(1,2); 
                        end
                    else 
                        band{dr} = repmat([(binranges(end)/2)-(band_size/2) (binranges(end)/2)+(band_size/2)], 2, 1); 
                    end
                end

            else 

                
                % Use raster cut times (includes the entire trial length)
                for j=1:C 
                    [pks,loc,~] = findpeaks(psth{actual_ds_Inds(nc)}{end,j},binranges,'WidthReference','HalfHeight','SortStr','descend');
                    if isempty(loc)
                        cuttime(j) = stim_dur./2; 
                    elseif ~isempty(loc) && length(loc)==1 
                        cuttime(j) = stim_dur./2; 
                    elseif ~isempty(pks) && length(pks)>1
                        cuttime(j) = ((loc(1)+loc(2))./2); 
                    end
                    band{j}(1,:) = [0 cuttime(j)];
                    band{j}(2,:) = [cuttime(j) binranges(end)];
                end

            end
            

            for k = 1:length(contrast_lev)

                for dr = 1:length(dirs)      
                    % Get spike count 
                    spk_times = cell2mat(rasterall{actual_ds_Inds(nc)}{k,dr}(trial_nums));
                    on_spikes{nc}(k,dr) = sum(spk_times>band{dr}(1,1) & spk_times<band{dr}(1,2));
                    off_spikes{nc}(k,dr) = sum(spk_times>=band{dr}(2,1) & spk_times<band{dr}(2,2));

                    all_spikes{nc}(k,dr) = on_spikes{nc}(k,dr)+off_spikes{nc}(k,dr); 
                    residual_spikes{nc}(k,dr) = length(spk_times) - all_spikes{nc}(k,dr); 
                end                            

                % Get angular separation between ON and OFF preferred
                % directions
                on_X{nc}(k) = sum(on_spikes{nc}(k,:).*cosd(dirs)); 
                on_Y{nc}(k) = sum(on_spikes{nc}(k,:).*sind(dirs));
                off_X{nc}(k) = sum(off_spikes{nc}(k,:).*cosd(dirs)); 
                off_Y{nc}(k) = sum(off_spikes{nc}(k,:).*sind(dirs));
                [th_ON{nc}(k), ~] = cart2pol(on_X{nc}(k), on_Y{nc}(k));
                [th_OFF{nc}(k), ~] = cart2pol(off_X{nc}(k), off_Y{nc}(k));
                th_ON{nc}(k) = wrapTo360(th_ON{nc}(k) * 180/pi);                 
                th_OFF{nc}(k) = wrapTo360(th_OFF{nc}(k) * 180/pi); 
                deltaangle_on_off{nc}(k) = abs(th_ON{nc}(k) - th_OFF{nc}(k)); 

                % Get circular std of ON and OFF response 
                bin_angle = diff(dirs(1:2));
                [angdev, circstdon] = circ_std(dirs * (pi/180), on_spikes{nc}(k,:), bin_angle*pi/180, 2);
                on_circstd{nc}(k) = circstdon*180/pi;
                [angdev, circstdoff] = circ_std(dirs * (pi/180), off_spikes{nc}(k,:), bin_angle*pi/180, 2);
                off_circstd{nc}(k) = circstdoff*180/pi;

                % Get circular std of Total response (degree)
                nd = wrapTo360(dirs(k_indx)-180); % null direction
                nullindx = find(dirs==nd); 
                [angdev, circstd] = circ_std(dirs * (pi/180), all_spikes{nc}(k,:), bin_angle*pi/180, 2);
                allDS_circstd{nc}(k) = circstd*180/pi; 

                % Get PD and ND at highest contrast (will be used outside this
                % loop)
                rr(k,:) = all_spikes{nc}(k,:); % will be used later for determining the response along PD and NULL directions
                if k==length(contrast_lev)
                    PD_angRad = ds_struct.angle{1,1,k}(actual_ds_Inds(nc));
                    PD_ang360 = wrapTo360(PD_angRad*180/pi);  % preferred direction in degree
                    max_resp_angs = dirs([find(dirs<PD_ang360, 1, 'last') find(dirs>PD_ang360,1,'first')]); % two nearest PDs
                    if length(max_resp_angs)==1 % catches error during 360 deg cross-over
                        max_resp_angs = [max_resp_angs 0];
                    end
                    [~,max_resp_indices,~] = intersect(dirs, max_resp_angs) ;
                    ND_ang360 = wrapTo360(PD_ang360+180); % null direction in degree
                    min_resp_angs = wrapTo360(max_resp_angs+180); % two nearest NDs
                    [~,min_resp_indices,~] = intersect(dirs, min_resp_angs) ;                        
                end

                % Get PD and ND at highest contrast separately for ON and OFF phases of response(will be used outside this
                % loop)
                if k==length(contrast_lev)
                    % ON
                    PD_ang360_ON = th_ON{nc}(k); 
                    max_resp_angs_ON = dirs([find(dirs<PD_ang360_ON, 1, 'last') find(dirs>PD_ang360_ON,1,'first')]); % two nearest PDs
                    if length(max_resp_angs_ON)==1 % catches error during 360 deg cross-over
                        max_resp_angs_ON = [max_resp_angs_ON 0];
                    end
                    [~,max_resp_indices_ON,~] = intersect(dirs, max_resp_angs_ON) ;
                    ND_ang360_ON = wrapTo360(PD_ang360_ON+180); % null direction in degree
                    min_resp_angs_ON = wrapTo360(max_resp_angs_ON+180); % two nearest NDs
                    [~,min_resp_indices_ON,~] = intersect(dirs, min_resp_angs_ON) ;

                    % OFF 
                    PD_ang360_OFF = th_OFF{nc}(k);
                    max_resp_angs_OFF = dirs([find(dirs<PD_ang360_OFF, 1, 'last') find(dirs>PD_ang360_OFF,1,'first')]); % two nearest PDs
                    if length(max_resp_angs_OFF)==1 % catches error during 360 deg cross-over
                        max_resp_angs_OFF = [max_resp_angs_OFF 0];
                    end
                    [~,max_resp_indices_OFF,~] = intersect(dirs, max_resp_angs_OFF) ;
                    ND_ang360_OFF = wrapTo360(PD_ang360_OFF+180); % null direction in degree
                    min_resp_angs_OFF = wrapTo360(max_resp_angs_OFF+180); % two nearest NDs
                    [~,min_resp_indices_OFF,~] = intersect(dirs, min_resp_angs_OFF) ;
                end

            end

            % Get tuning curves and tuning width by centering peak around PD
            % obtained from directional response @ highest contrast. The
            % location of this peak should be invariant to contrast 
            [~,ind_temp] = max(cosd(max_resp_angs - PD_ang360)); % index of the max_resp_angs that is closest to the PD (@ highest contrast)
            dir_near_PD = max_resp_angs(ind_temp); 
            [~,ind_max] = intersect(dirs, dir_near_PD); % index of "dirs" to use as proxy for PD

            for k=1:length(contrast_lev)
                % Get tuning curves and tuning width for Total Response
                clear tempvar;
                tempvar = all_spikes{nc}(k,:); 
                ctr = round(length(dirs)/2); 
                if ind_max>ctr
                    tempvar = [tempvar(ind_max-ctr+1:end)  tempvar(1:(ind_max-ctr))]; 
                elseif ind_max<=ctr
                    tempvar = [tempvar(end-(ctr-ind_max)+1:end) tempvar(1:end-(ctr-ind_max))];
                end

                try 
                    [~, ~, widths, ~] = findpeaks(tempvar-min(tempvar),dirs_cent,'WidthReference','halfheight','SortStr','descend');
                    allDS_tuningwid{nc}(k) = widths(1);
                catch ME   
                    allDS_tuningwid{nc}(k) = 0;
                end
                allDS_tuningCurves{nc}(k,:) = tempvar;

                % Get tuning curves and tuning width for ON Response
                clear tempvar;
                tempvar = on_spikes{nc}(k,:);  
                ctr = round(length(dirs)/2); 
                if ind_max>ctr
                    tempvar = [tempvar(ind_max-ctr+1:end)  tempvar(1:(ind_max-ctr))]; 
                elseif ind_max<=ctr
                    tempvar = [tempvar(end-(ctr-ind_max)+1:end) tempvar(1:end-(ctr-ind_max))];
                end

                try 
                    [~, ~, widths, ~] = findpeaks(tempvar-min(tempvar),dirs_cent,'WidthReference','halfheight','SortStr','descend');
                    on_phase_tuningwid{nc}(k) = widths(1);
                catch ME   
                    on_phase_tuningwid{nc}(k) = 0;
                end
                on_phase_tuningCurves{nc}(k,:) = tempvar;
                
                % Get tuning curves and tuning width for OFF Response 
                clear tempvar;
                tempvar = off_spikes{nc}(k,:);  
                ctr = round(length(dirs)/2); 
                if ind_max>ctr
                    tempvar = [tempvar(ind_max-ctr+1:end)  tempvar(1:(ind_max-ctr))]; 
                elseif ind_max<=ctr
                    tempvar = [tempvar(end-(ctr-ind_max)+1:end) tempvar(1:end-(ctr-ind_max))];
                end

                try 
                    [~, ~, widths, ~] = findpeaks(tempvar-min(tempvar),dirs_cent,'WidthReference','halfheight','SortStr','descend');
                    off_phase_tuningwid{nc}(k) = widths(1);
                catch ME   
                    off_phase_tuningwid{nc}(k) = 0;
                end
                off_phase_tuningCurves{nc}(k,:) = tempvar;
            end

            % Get tuning strength (note: PD is determined from highest contrast)
            for k=1:length(contrast_lev)
                pd_spkcnt{nc}(k,1) = sum(rr(k,max_resp_indices).*cosd(max_resp_angs-PD_ang360))./sum(cosd(max_resp_angs-PD_ang360));  % approximate number of spikes in the PD
                nd_spkcnt{nc}(k,1) = sum(rr(k,min_resp_indices).*cosd(min_resp_angs-ND_ang360))./sum(cosd(min_resp_angs-ND_ang360));  % approximate number of spikes in the ND
                tuning_strength{nc}(k,1) = (pd_spkcnt{nc}(k,1) - nd_spkcnt{nc}(k,1))./((pd_spkcnt{nc}(k,1) + nd_spkcnt{nc}(k,1)));
                if isinf(tuning_strength{nc}(k,1)) || isnan(tuning_strength{nc}(k,1)) 
                    tuning_strength{nc}(k,1) = 0; 
                end
            end

            % Get tuning strength for ON and OFF phases of response (note: PD is determined from highest contrast separately for ON and OFF)
            rrON = on_spikes{nc}; 
            for k=1:length(contrast_lev)
                pd_spkcnt_ON{nc}(k,1) = sum(rrON(k,max_resp_indices_ON).*cosd(max_resp_angs_ON-PD_ang360_ON))./sum(cosd(max_resp_angs_ON-PD_ang360_ON));  % approximate number of spikes in the PD
                nd_spkcnt_ON{nc}(k,1) = sum(rrON(k,min_resp_indices_ON).*cosd(min_resp_angs_ON-ND_ang360_ON))./sum(cosd(min_resp_angs_ON-ND_ang360_ON));  % approximate number of spikes in the ND
                tuning_strength_ON{nc}(k,1) = (pd_spkcnt_ON{nc}(k,1) - nd_spkcnt_ON{nc}(k,1))./((pd_spkcnt_ON{nc}(k,1) + nd_spkcnt_ON{nc}(k,1)));
                if isinf(tuning_strength_ON{nc}(k,1)) || isnan(tuning_strength_ON{nc}(k,1)) 
                    tuning_strength_ON{nc}(k,1) = 0; 
                end
            end
            rrOFF = off_spikes{nc};
            for k=1:length(contrast_lev)
                pd_spkcnt_OFF{nc}(k,1) = sum(rrOFF(k,max_resp_indices_OFF).*cosd(max_resp_angs_OFF-PD_ang360_OFF))./sum(cosd(max_resp_angs_OFF-PD_ang360_OFF));  % approximate number of spikes in the PD
                nd_spkcnt_OFF{nc}(k,1) = sum(rrOFF(k,min_resp_indices_OFF).*cosd(min_resp_angs_OFF-ND_ang360_OFF))./sum(cosd(min_resp_angs_OFF-ND_ang360_OFF));  % approximate number of spikes in the ND
                tuning_strength_OFF{nc}(k,1) = (pd_spkcnt_OFF{nc}(k,1) - nd_spkcnt_OFF{nc}(k,1))./((pd_spkcnt_OFF{nc}(k,1) + nd_spkcnt_OFF{nc}(k,1)));
                if isinf(tuning_strength_OFF{nc}(k,1)) || isnan(tuning_strength_OFF{nc}(k,1)) 
                    tuning_strength_OFF{nc}(k,1) = 0; 
                end
            end
        end



        % Return stats: All o-DS and oo-DS cells
        retstat.stim_dur = stim_dur; 
        retstat.dirs = dirs; 
        retstat.dirs_cent = dirs_cent;

        retstat.total_spikes_pd_DS = pd_spkcnt;
        retstat.total_spikes_nd_DS = nd_spkcnt;
        retstat.tuning_strength = tuning_strength;
        retstat.non_ds_id = non_ds_id; 
        retstat.non_ds_ind = non_ds_ind;
        retstat.total_spikes_nonDS = total_spikes_nonDS;  
        retstat.allDS_tuningCurves = allDS_tuningCurves;
        retstat.allDS_tuningwid = allDS_tuningwid;
        retstat.on_off_angdisp_ds = allDS_circstd;

        retstat.on_phase_spkcnt = on_spikes; 
        retstat.off_phase_spkcnt = off_spikes; 
        retstat.onoff_phase_spkcnt = all_spikes;
        retstat.residual_spkcnt = residual_spikes; 
        retstat.on_X = on_X;
        retstat.on_Y = on_Y;
        retstat.off_X = off_X;
        retstat.off_Y = off_Y;    
        retstat.pd_ang_ON = th_ON; % 
        retstat.pd_ang_OFF = th_OFF;
        retstat.on_off_deltaangle = deltaangle_on_off;    
        retstat.on_phase_angdisp_ds = on_circstd; 
        retstat.off_phase_angdisp_ds = off_circstd;   
        retstat.on_phase_tuningCurves = on_phase_tuningCurves; 
        retstat.off_phase_tuningCurves = off_phase_tuningCurves; 
        retstat.on_phase_tuningwid = on_phase_tuningwid;
        retstat.off_phase_tuningwid = off_phase_tuningwid;

        retstat.on_phase_tuningstrength = tuning_strength_ON;
        retstat.off_phase_tuningstrength = tuning_strength_OFF;


        %------- Separate calculation for SUPERIOR, INFERIOR, POSTERIOR, ANTERIOR cells conditional on their existence ------------% 
        if (isfield(retstat,'oo_ds_Inferior') || isfield(retstat,'oo_ds_Superior') ... 
                || isfield(retstat,'oo_ds_Posterior') || isfield(retstat,'oo_ds_Anterior'))


            % Set up field names of DS subtypes 
            names = {'oo_ds_Inferior','oo_ds_Superior','oo_ds_Posterior','oo_ds_Anterior'}; 
            cnt =1;
            for u=1:length(names)
                if isfield(retstat, names{u})==1 
                    strname{cnt} = names{u}(7:end); 
                    namesAc{cnt} = names{u}(7:9);
                    cnt=cnt+1; 
                end
            end


            % Extract indices within DS cell population
            for u=1:length(strname)
                if strcmp(strname{u},'Inferior') 
                    [~,indInf] = intersect(actual_ds_Ids, retstat.oo_ds_Inferior); 
                elseif strcmp(strname{u},'Superior') 
                    [~,indSup] = intersect(actual_ds_Ids, retstat.oo_ds_Superior); 
                elseif strcmp(strname{u},'Posterior') 
                    [~,indPos] = intersect(actual_ds_Ids, retstat.oo_ds_Posterior);   
                elseif strcmp(strname{u},'Anterior') 
                    [~,indAnt] = intersect(actual_ds_Ids, retstat.oo_ds_Anterior);
                end
            end


            % Extract stats for Inferior, Superior, Anterior, Posterior
            % subtypes 

            for u=1:length(namesAc)

                if strcmp(namesAc{u},'Inf') 
                    if ~isempty(indInf) 
                        bool_cond = true; 
                        indType = indInf; 
                    end
                elseif strcmp(namesAc{u},'Sup')
                    if ~isempty(indSup)
                        bool_cond = true; 
                        indType = indSup; 
                    end
                elseif strcmp(namesAc{u},'Pos') 
                    if ~isempty(indPos)
                        bool_cond = true; 
                        indType = indPos; 
                    end
                elseif strcmp(namesAc{u},'Ant') 
                    if ~isempty(indAnt)
                        bool_cond = true; 
                        indType = indAnt; 
                    end
                end

                if bool_cond
                    retstat.(strname{u}).stim_dur = stim_dur; 
                    retstat.(strname{u}).dirs = dirs; 
                    retstat.(strname{u}).dirs_cent = dirs_cent;

                    retstat.(strname{u}).on_phase_spkcnt = on_spikes{indType}; 
                    retstat.(strname{u}).off_phase_spkcnt = on_spikes{indType}; 
                    retstat.(strname{u}).onoff_phase_spkcnt = all_spikes{indType}; 
                    retstat.(strname{u}).on_X = on_X{indType};
                    retstat.(strname{u}).on_Y = on_Y{indType};
                    retstat.(strname{u}).off_X = off_X{indType};
                    retstat.(strname{u}).off_Y = off_Y{indType};  
                    retstat.(strname{u}).pd_ang_ON = th_ON{indType};  
                    retstat.(strname{u}).pd_ang_OFF = th_OFF{indType};  
                    retstat.(strname{u}).on_off_deltaangle = deltaangle_on_off{indType};
                    retstat.(strname{u}).on_phase_angdisp_ds = on_circstd{indType};
                    retstat.(strname{u}).off_phase_angdisp_ds = off_circstd{indType};
                    retstat.(strname{u}).on_phase_tuningCurves = on_phase_tuningCurves{indType};
                    retstat.(strname{u}).off_phase_tuningCurves = off_phase_tuningCurves{indType};
                    retstat.(strname{u}).on_phase_tuningwid = on_phase_tuningwid{indType};
                    retstat.(strname{u}).off_phase_tuningwid = off_phase_tuningwid{indType};  
                end

            end
        end

        else
        warning('Skipping creating ON, OFF stats !');
    end


end



%% Return spike raster 
%
% Input arguments: 
%               spike_times: spike times in sec 
%               trial_begin_times: list of times at which each new trial begins 
%               optional arguments: 
%                                   start           0
%                                   stop           []
%                                   first_trial     1
%                                   tic_color       'k'
%                                   tic_thickness   1
%                                   line_space      0.2
%                                   axis_range      []
%                                   first_tic       1
%
% Output arguments: 
%               return_raster: cell array storing the spike times of each trial
%
%
% Copyright: Sneha Ravi, 2013


function return_raster = get_raster(spike_times, trial_begin_times, varargin)

    % SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParameter('start', 0, @isnumeric);
    p.addParameter('stop', [], @isnumeric);
    p.addParameter('tic_color', [0 0 0]);
    p.addParameter('line_space', 0.2, @isnumeric)
    p.addParameter('tic_thickness', 1);
    p.addParameter('first_trial', 1, @isnumeric)
    p.addParameter('labels', true);
    p.addParameter('axis_range', [0 20 0 50]);
    p.addParameter('first_tic', 1, @isnumeric) ;

    % parse inputs
    p.parse(varargin{:});
    params = p.Results;

    if isempty(params.stop)
        params.stop = mean(diff(trial_begin_times));
    end

    % setup the figure
    plot_axes = gca; 
    xlabel('seconds');
    ylabel('trials');


    % calculate the number of trials to plot
    num_trials = length(trial_begin_times)-params.first_trial+1;

    % initialize returned variable and counter
    return_raster = cell(num_trials, 1);
    trial_plot_line = params.first_tic; % counter

    for i = params.first_trial:num_trials
        % shift spikes
        ith_shifted_times = spike_times - trial_begin_times(i);  
        % get spikes in trial frame
        ith_trial_times = ith_shifted_times(ith_shifted_times >= params.start & ith_shifted_times <= params.stop);
        % plot the spike times
        if params.plot
            if ~isempty(ith_trial_times)
                line([ith_trial_times, ith_trial_times], [trial_plot_line-1+params.line_space, trial_plot_line-params.line_space],'color', params.tic_color, 'LineWidth', params.tic_thickness)
                axis(params.axis_range);
                set(plot_axes,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1);
            end
        end
        % keep track of trial plotted
        trial_plot_line = trial_plot_line + 1; 
        % store spike times
        return_raster{i} = ith_trial_times;
    end
end












