%{
   For Figure 2C: This script will plot subject-specific brain mesh, network parcellations,
   and math active and deactive electrored that pass the ROL post-hoc)
%}

% Example:  Math_Activation_Plot_Fig2C(subject_str,EOI_All_keep_in,{'l','lm','li'})
function Math_Activation_Plot_Fig2A_new(subject_str,sig_PMC_elecs,sig_LPC_elecs,views,PMC_Elecs, SPL_Elecs,Elec_List)

if ~exist(['/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Results_Figures',filesep,'fig_2C'],'dir')
            mkdir(pwd,'fig_2C');
end

parcOut=elec2Parc(subject_str,'Y7');

color_matrix = zeros(length(parcOut),3);
elecNames = cell(length(parcOut),1);

p_to_plot = [];

% Create color matrix
for p = 1:length(parcOut)
    
    index = find(ismember(Elec_List(:,1),parcOut{p,1}));
    elec = Elec_List{index,2};
    elec = str2num(elec(6:end));
    
    if ismember(elec,SPL_Elecs) || ismember(elec,PMC_Elecs)

        color_matrix(p,:) = [0    0    0];
      %  p_to_plot = [p_to_plot, p];

        
        if ismember(elec,sig_PMC_elecs) 
            
            color_matrix(p,:) = [0,0.513725490196078,0.611764705882353];
            p_to_plot = [p_to_plot, p];
            
        elseif ismember(elec,sig_LPC_elecs) 
            
            color_matrix(p,:) = [0.674509803921569,0.145098039215686,0.184313725490196];
            p_to_plot = [p_to_plot, p];

        end
        
        
    else
        
        color_matrix(p,:) =  [.5    .5    .5];
        
    end
    
    elecNames{p,1} = parcOut{p,2};

end

    
% look through views
for i = 1:length(views)

 
    [averts, label, col]=read_annotation(fullfile(getFsurfSubDir(),'fsaverage','label','lh.Yeo2011_17Networks_N1000.annot'));
    %Network you want to plot
   % id = [4 8];
    %Make verteces gray
   % parc_col = .7.*255.*ones(size(col.table(:,1:3)));
    %Color parcel of interest
    col.table(4,1:3) = [0 170 14]; % change green to light green
    col.table(8,1:3) = [230 62 78]; % change green to light green
  %  parc_col(id,:)=col.table(id,1:3);
    cfg=[];
    cfg.view= views{i}; % l / lm / li
%     
%     brainView.light=[0 0 0];
%     brainView.hem='l';
%     brainView.eyes=[25 0]
%     cfg.view=brainView
   
    cfg.showLabels='n';
    %cfg.overlayParcellation='Y17';
    cfg.pullOut = 0;
  %  cfg.parcellationColors = parc_col;
    cfg.elecSize=2.4;
    cfg.title= [];
    cfg.onlyShow = parcOut(p_to_plot,1);
    cfg.elecShape='sphere';
    %elecNames = parcOut(:,2); % will need to change for each subject
    cfg.elecColorScale =[0 1];
    cfg.elecColors= color_matrix(:,1:3);   
    %cfg.elecNames=elecNames;
    delete(findall(gcf,'Type','light'))
    cfgOut=plotPialSurf(subject_str,cfg);

    print('-opengl','-r300','-dpng',strcat([pwd,filesep,'fig_2C',filesep,'Fig_Math_Resposes_2C_',subject_str,'_view:',char(views{i})]));
end