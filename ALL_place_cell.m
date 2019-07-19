function ALL_place_cell
%runs place cell function on every cell

count = 0;

%LOCATE SUBJECT FOLDERS
%
file_list_subjects = dir('neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:}, 1);

%ITERATE THROUGH SUBJECTS
%
for subject = 1:length_subjects
    current_rat = file_names_subjects{:}(subject,1).name
    file_list_session_type = dir(['neurodata/' num2str(current_rat) '/']);
    file_list_session_type(1:2) = [];
    file_names_session_types = {file_list_session_type([file_list_session_type(:).isdir])};
    length_session_types = size(file_names_session_types{:}, 1);
    
    %ITERATE THROUGH SESSION TYPE
    %
    for session_type = 1%:length_session_types
        current_session_type = file_names_session_types{:}(session_type,1).name
        file_list_sessions = dir(['neurodata/' num2str(current_rat) '/' num2str(current_session_type) '/*.mat']);
        %any2csv(file_list_sessions, '|', 1)
        length_sessions = size(file_list_sessions,1);
        
        %ITERATE THROUGH SESSIONS
        %
        for session = 1:length_sessions
            day_file = file_list_sessions(session,1).name
            
            %load session
            load(['/Users/ampm/Documents/MATLAB/PCia/neurodata/' num2str(current_rat) '/' num2str(current_session_type) '/' num2str(day_file)]);

            if isempty(clusters)
                continue
            end
           
            %INSERT FUNCTION HERE
            clusters(:,1)'
            clusts = unique(eptrials(~isnan(eptrials(:,4)), 4));
            clusts = clusts(clusts>1)'
            for cf = clusts 
                
                %{
                count = count +1
                if count < 10
                    continue
                end
                %}
                
                
                figure; hold on
                
                vis_count = 0;
                for visit = sort(context_ids)'
                    vis_count = vis_count+1;
                    %current visit
                    visit_idx = eptrials(:,6)==visit;
                    
                    %plot heatmap
                    subplot(1,4, visit);
                    [~,~,~,rm] = rate_mtx(eptrials(visit_idx,:), cf, 50); 
                    rm = inpaint_nans(smooth2a(rm,1));
                    imagesc(rm); 
                    
                    if vis_count == 1
                        rm_min = min(rm(:));
                        rm_max = max(rm(:));
                    else     
                        rm_min = min([rm_min min(rm(:))]);
                        rm_max = max([rm_max max(rm(:))]);
                    end
                        
                        
                        
                    axis square; axis off
                    title(num2str(mode(eptrials(visit_idx,6))))
                      
                    %place cell info
                    [H_hi, cc_hi] = place_cell(eptrials(visit_idx,:), cf, 50, 0, 1);
                    [H_lo, cc_lo] = place_cell(eptrials(visit_idx,:), cf, 50, 0, 0);
                    
                    %plot place cell coords
                    if H_hi == 1

                        for pc = 1:length(cc_hi)
                        
                            %remove padding from coords
                            cc_hi{pc} = cc_hi{pc} - ones(size(cc_hi{pc}));
                            %correct for bins
                            cc_hi{pc}(1,:) = cc_hi{pc}(1,:) + repmat(.5, size(cc_hi{pc}(1,:)));
                            cc_hi{pc}(2,:) = cc_hi{pc}(2,:) + repmat(.5, size(cc_hi{pc}(2,:)));
                            %smooth coords
                            smooth_wndw = 8;
                            cc_hi_sm(1,:) = smooth([cc_hi{pc}(1,:) cc_hi{pc}(1,1)],smooth_wndw);
                            cc_hi_sm(2,:) = smooth([cc_hi{pc}(2,:) cc_hi{pc}(2,1)],smooth_wndw);

                            %plot
                            hold on; plot(cc_hi_sm(1,:), cc_hi_sm(2,:), 'k-', 'linewidth', 2.5)
                        end
                    end
                    
                    if H_lo == 1
                        
                        for pc = 1:length(cc_hi)
                        
                            %remove padding from coords
                            cc_lo{pc} = cc_lo{pc} - ones(size(cc_lo{pc}));
                            %correct for bins
                            cc_lo{pc}(1,:) = cc_lo{pc}(1,:) + repmat(.5, size(cc_lo{pc}(1,:)));
                            cc_lo{pc}(2,:) = cc_lo{pc}(2,:) + repmat(.5, size(cc_lo{pc}(2,:)));
                            %smooth coords
                            smooth_wndw = 8;
                            cc_lo_sm(1,:) = smooth([cc_lo{pc}(1,:) cc_lo{pc}(1,1)],smooth_wndw);
                            cc_lo_sm(2,:) = smooth([cc_lo{pc}(2,:) cc_lo{pc}(2,1)],smooth_wndw);

                            %plot
                            hold on; plot(cc_lo_sm(1,:), cc_lo_sm(2,:), 'k-', 'linewidth', 2.5)
                        end
                    end
                    
                    if vis_count == 4
                        for viz = 1:4
                            subplot(1,4, viz)
                            colorbar
                            caxis([rm_min*1.1      rm_max*.9])
                        end
                    end
                    clearvars cc_hi_sm cc_lo_sm
                    
                end
                
                %print(strcat('/Users/ampm/Desktop/pcia_figs/all_heatmaps/place/', num2str(current_rat),'_' ,num2str(current_session_type),'_' ,day_file(1:end-4), '_', num2str(cf*100), '.pdf'), '-dpdf')
                print(strcat('/Users/ampm/Desktop/pcia_figs/all_heatmaps/obj_rate/', num2str(current_rat),'_' ,num2str(current_session_type),'_' ,day_file(1:end-4), '_', num2str(cf*100), '.pdf'), '-dpdf')
                close all
                
                
                
                %{
                [H] = place_cell(eptrials, cf, 50, 2, 1); 
                H_hi = H; display(H_hi); 
                [H] = place_cell(eptrials, cf, 50, 2, 0); 
                H_lo = H; 
                display(H_lo); 
                %}
            end
        
            
        end
    end
end


end