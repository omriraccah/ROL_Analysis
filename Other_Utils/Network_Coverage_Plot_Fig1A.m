%{
   For Figure 1A: This script will plot subject-specific brain mesh, network parcellations,
   and electrodes on top of those parcellations (black over networks and grey over the rest of the brain)
%}

% Example:  Network_Coverage_Plot_Fig1A('S12_38_LK',parcOut, {'l','lm','li'}) 
function Coverage_Plot_Fig2A(subject_str,elec, views)

if ~exist(['/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Results_Figures',filesep,'fig_1A'],'dir')
            mkdir(pwd,'fig_1A');
end

 parcOut=elec2Parc(subject_str,'Y7');

color_matrix = zeros(length(parcOut),3);
showOnly = {};
counter = 0;

% Create color matrix
for p = 1:length(parcOut)
%     
%     % get electrode number
%     electrode = parcOut{p,1};
%     electrode = electrode(1,1);
%     electrode = str2double(electrode(6:end));
        
%     % loop through math-active electrodes 
    for e = 1:length(elec.mathA_LPC_selective_elecs)
        
        
        if electrode == elec.mathA_LPC_selective_elecs(e)
            
            % make red
            color_matrix(p,:) = [0.674509803921569,0.145098039215686,0.184313725490196];
            counter = counter+1;

            showOnly{1,counter} = parcOut{p,2};
            
        end
        
    end
    
    % loop through attention electrodes 
    for e = 1:length(elec.mathA_LPC_elecs_Attention)
        
        
        if electrode == elec.mathA_LPC_elecs_Attention(e)
            
            % make pink
            color_matrix(p,:) = [0.909803921568627,0.541176470588235,0.482352941176471];
            counter = counter+1;

            showOnly{1,counter} = parcOut{p,2};
            
        end
        
    end
    
   
    
    % loop through autobio electrodes
    for e = 1:length(elec.autoA_PMC_elecs)
        
        
        if electrode == elec.autoA_PMC_elecs(e)
            
            % make yellow
            color_matrix(p,:) = [0.972549019607843,0.803921568627451,0.345098039215686];
            counter = counter+1;

            showOnly{1,counter} = parcOut{p,2};
            
        end
        
    end
    
    % loop through math-deactive electrodes 
    for e = 1:length(elec.mathD_PMC_elecs)
        
        
        if electrode == elec.mathD_PMC_elecs(e)
            
            % make blue
            color_matrix(p,:) = [0.117647058823529,0.662745098039216,0.788235294117647];
            counter = counter+1;

            showOnly{1,counter} = parcOut{p,2};
            
        end
        
    end
%     
    %loop through math-deactive electrodes 
    for e = 1:length(elec.Visual_elecs)
        
        
        if electrode == elec.Visual_elecs(e)
            
           % make green
            color_matrix(p,:) = [0.0392156862745098,0.654901960784314,0.403921568627451];
            counter = counter+1;

            showOnly{1,counter} = parcOut{p,2};
            
        end
        
    end
    
%      % loop through math-deactive electrodes 
%     for e = 1:length(elec.autoA_mathD_PMC_elecs)
%         
%         
%         if electrode == elec.autoA_mathD_PMC_elecs(e)
%             
%             % make black
%             color_matrix(p,:) = [0 0 0];
%             counter = counter+1;
% 
%             showOnly{1,counter} = parcOut{p,2};
%             
%         end
%         
%     end
  
    

end
    
    
% look through views
for i = 1:length(views)

   %[averts, label, col]=read_annotation(fullfile(getFsurfSubDir(),'fsaverage','label','lh.Yeo2011_7Networks_N1000.annot'));
    %Network you want to plot
    id = [4 8];
    %Make verteces gray
    %parc_col = .7.*255.*ones(size(col.table(:,1:3)));
    %Color parcel of interest
   parc_col(id,:)=col.table(id,1:3);
    cfg=[];
    cfg.view= views{i}; % l / lm / li
    cfg.showLabels='n';
   %cfg.overlayParcellation='Y7';
    cfg.pullOut = 1.5;
    %cfg.parcellationColors = parc_col;
    cfg.title= [];
    cfg.onlyShow = showOnly;
    cfg.elecSize=3;
    cfg.elecShape='sphere';
    %elecNames = parcOut(:,2); % will need to change for each subject
    cfg.elecColorScale =[0 1];
    cfg.elecColors= color_matrix(:,1:3);   
    %cfg.elecNames=elecNames;
    cfgOut=plotPialSurf(subject_str,cfg);

    print('-opengl','-r300','-dpng',strcat([pwd,filesep,'fig_1A',filesep,'Fig_1A_',subject_str(1:7),'_final',char(views{i})]));
end






