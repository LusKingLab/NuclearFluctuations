%% Verify movie quality.

clear all; clc; close all;

% Fill out and verify the following until the dashed line:
rootpath = './'; % Current folder where data is stored and where you
                        % want to write analysis data.
verifypath = [rootpath,'verification_epe1_clr4_swi6_sm1']; % Name of folder in which
                                                  % to write verification
                                                  % data.
mkdir(verifypath); % Make verification folder.
run setup_header_epe1_clr4_swi6_sm1.m; % Designate setup_header file with list of 
                              % data to process, colors, etc.
allfiles = dir(fullfile('./','data')); % .mat files to be analyzed 
                                              % should be in the following
                                              % folder sequence (example): 
                                              % ./nuclei/data/wt/wt.mat
allnamestmp = {allfiles(:).name};
allnames = {};
    
for i = 1:length(allnamestmp)
    if startsWith(allnamestmp(i),'.') == 0   
        allnames(i) = allnamestmp(i); % Ignore cached file that does not 
                                      % contain pertinent data and cannot 
                                      % be read. - JJ 11/19/18
    end
end

allnames = allnames(~cellfun('isempty',allnames)); % Delete empty cells.
% Select only the folder names of data you want to compare.
allnames = {allnames{1},... 
            allnames{2},...
            allnames{3},...
            allnames{4},...
            allnames{12},...
            allnames{13},...
            allnames{14},...
            allnames{15},...
            allnames{16},...
            allnames{17}};
                               
% -------------------------------------------------------------------------

% Setup parameters; do not change...                             
[points,faces,edges,neighbors] = TriSphere(3); % TRISPHERE: Returns the 
                                               % triangulated model of a 
                                               % sphere using the 
                                               % icosaedron subdivision 
                                               % method.
neighbors(1:12,6) = (1:12)';
zrange = find(abs(points(:,3))<0.5);
gnmovienames = {};
goodnuclei = [];

% Set plot colors for final cumulative RMSF plot using colors from setup
% header. 
allcolors = [230,159,0;...   % clr4
             230,159,0;...   % clr4 MBC
             0,114,178;...   % epe1
             0,114,178;...   % epe1 MBC
             0,158,115;...   % swi6
             213,94,0;...    % swi6-sm1
             213,94,0;...    % swi6-sm1 MBC
             0,158,115;...   % swi6 MBC  
             0,0,0;...       % wt
             0,0,0]/...      % wt MBC
             255;

% Open figure windows.
f1 = figure(3001);
f2 = figure(3002);
f3 = figure(3003);
f4 = figure(3004);
f5 = figure(3005);
f6 = figure(3006);
f7 = figure(3007);
f8 = figure(3008);

% Display each data set.
for itype = 1:length(allnames)
    clf(f1);
    clf(f2);
    clf(f3);
    clf(f4);
    clf(f5);
    clf(f7);
    clf(f8);

    typermsf = [];
    
    % Find .mat files for each data set.
    moviefiles = dir(fullfile(rootpath,'data',allnames{itype},'*.mat'));

    movienamestmp = {moviefiles.name};
    movienames = {};
    
    for i = 1:length(movienamestmp)
        if startsWith(movienamestmp(i),'.') == 0   
            movienames(i) = movienamestmp(i); % Ignore cached file that 
                                              % does not contain pertinent 
                                              % data and cannot be read. 
                                              % - JJ 11/19/18
        end
    end
    
    movienames = movienames(~cellfun('isempty',movienames));
    
    legendids = movienames;
    moviecolor = parula(length(movienames));
    for imovie = 1:length(movienames)
        load(fullfile(rootpath,'data',allnames{itype},movienames{imovie}));
        display(['processing ',movienames{imovie}]);
        goodnucleitmp=zeros(1,nm.num_nuc);
        %%
        rmsf=zeros(length(zrange),nm.num_nuc);
        dr2s=zeros(length(zrange),nm.num_nuc);
        drs=zeros(nm.num_nuc,length(zrange),nm.endframe);
        existflags=zeros(nm.num_nuc,nm.endframe);
        xs=zeros(nm.num_nuc,nm.endframe);
        ys=zeros(nm.num_nuc,nm.endframe);
        zs=zeros(nm.num_nuc,nm.endframe);
        ozs=zeros(nm.num_nuc,nm.endframe);
        dcs=zeros(nm.num_nuc,nm.endframe);
        for inuc=1:nm.num_nuc
            r_s=zeros(length(zrange),nm.endframe);
            dr_s=zeros(length(zrange),nm.endframe);
            for iframe=1:nm.endframe
                nuc=nm.nuclei{iframe,inuc};
                allr=nuc.r_new;
                neighbor_r=allr(neighbors);
                dr2=sum((allr*ones(1,6)-neighbor_r).^2,2)/6;

                r=nuc.r_new(zrange);
                r_s(:,iframe)=r;
                dr_s(:,iframe)=dr2(zrange);
                existflags(inuc,iframe)=nuc.exitflag;
                xs(inuc,iframe)=nuc.origin_new(1);
                ys(inuc,iframe)=nuc.origin_new(2);
                zs(inuc,iframe)=nuc.origin_new(3);
            end
%             drs(inuc,:,:)=r_s-mean(r_s,2)*ones(1,size(r_s,2));
            xs(inuc,:)=xs(inuc,:)-mean(xs(inuc,:));
            ys(inuc,:)=ys(inuc,:)-mean(ys(inuc,:));
            ozs(inuc,:)=zs(inuc,:);
            zs(inuc,:)=zs(inuc,:)-mean(zs(inuc,:));
            dcs(inuc,:)=sqrt(xs(inuc,:).^2+ys(inuc,:).^2+zs(inuc,:).^2)*...
                p2um;
            rmsf(:,inuc)=std(r_s,1,2)*p2um;
            dr2s(:,inuc)=max(dr_s,[],2)';
            
            if max(ozs(inuc,:))<=8 && min(ozs(inuc,:))>=3 ...
                    && max(rmsf(:,inuc))<0.3 && mean(rmsf(:,inuc))<0.1 ...
                     && max(dcs(inuc,:))<0.6 ...
                     && max(dr2s(:,inuc))<0.5
                goodnucleitmp(inuc)=1;
            else
                goodnucleitmp(inuc)=0;
            end
        end
        gnmovienames=[gnmovienames;{allnames{itype},' ',...
            movienames{imovie}}];
        goodnuclei=[goodnuclei;{goodnucleitmp}];
        
        % Do good nuclei exist?
        figure(f1);
        plot(existflags(find(goodnucleitmp),:)','color',...
            moviecolor(imovie,:),'LineWidth',3);hold on;
        plot(existflags(find(~goodnucleitmp),:)','color',...
            moviecolor(imovie,:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('existflag');axis([0 101 -2 4])
        set(gca,'Fontsize',30,'LineWidth',3);

        % Z position of nuclei over time-lapse movie.
        figure(f2);
        plot(ozs(find(goodnucleitmp),:)','color',...
            moviecolor((imovie),:));hold on;
        plot(ozs(find(~goodnucleitmp),:)','color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('z-position (slice)');axis([0 101 0 11])
        set(gca,'Fontsize',30,'LineWidth',3);
        
        % X,Y drift.
        figure(f3);
        plot(dcs(find(goodnucleitmp),:)','color',...
            moviecolor((imovie),:));hold on;
        plot(dcs(find(~goodnucleitmp),:)','color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('drift (\mum)');axis([0 101 0 2])
        set(gca,'Fontsize',30,'LineWidth',3);
        
        % RMSF by surface angle.
        figure(f4)
        plot(rmsf(:,find(goodnucleitmp)),'color',...
            moviecolor((imovie),:));hold on;
        plot(rmsf(:,find(~goodnucleitmp)),'color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('angles');ylabel('rmsf (\mum)');...
            axis([0 length(zrange) 0 0.4])
        set(gca,'Fontsize',30,'LineWidth',3);
        
        % Cumulative probability RMSF by movie.
        figure(f5)
        bins=0:0.0025:0.15;
        rmsf1=rmsf(:,find(goodnucleitmp));
        [counts]=hist(rmsf1(:),bins);
        cumcounts=cumsum(counts)./sum(counts);
        plot(bins,cumcounts,'linewidth',2,'color',...
            moviecolor((imovie),:));hold on;
        xlabel('rmsf(\mum)');ylabel('cumulative probability');...
            axis([0 0.1 0 1])
        set(gca,'Fontsize',30,'LineWidth',3);
        legend(legendids,'Interpreter','none','FontSize',6,'Location',...
            'southeast'); % Set Interpreter property 
            % to 'none'; the default for text fields is LaTeX, e.g., '_' in
            % string becomes an subscript in legend text. - JJ, 11/15/2017
        
        figure(f7)
        plot(dr2s(:,find(goodnucleitmp)),'color',...
            moviecolor((imovie),:));hold on;
        plot(dr2s(:,find(~goodnucleitmp)),'color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('sum dr square over 6 (pixel^2)');...
            axis([0 length(zrange) 0 max(max(dr2s))])
        set(gca,'Fontsize',30,'LineWidth',3);
        
        % Plot all together in single figure window.
        figure(f8);set(f8,'Position',[0 0 1500 1000])
        subplot(2,3,1)
        plot(existflags(find(goodnucleitmp),:)','color',...
            moviecolor((imovie),:));hold on;
        plot(existflags(find(~goodnucleitmp),:)','color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('exist flag');axis([0 101 -2 4]);...
            title('existflag');
        set(gca,'Fontsize',30,'LineWidth',3);
        
        subplot(2,3,2)
        plot(ozs(find(goodnucleitmp),:)','color',...
            moviecolor((imovie),:));hold on;
        plot(ozs(find(~goodnucleitmp),:)','color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('z-position (slice)');...
            axis([0 101 0 11]);title('z-position');
        set(gca,'Fontsize',30,'LineWidth',3);
        
        subplot(2,3,3)
        plot(dcs(find(goodnucleitmp),:)','color',...
            moviecolor((imovie),:));hold on;
        plot(dcs(find(~goodnucleitmp),:)','color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('drift (\mum)');axis([0 101 0 2]);...
            title('drift');
        set(gca,'Fontsize',30,'LineWidth',3);
        
        subplot(2,3,4)
        plot(rmsf(:,find(goodnucleitmp)),'color',...
            moviecolor((imovie),:));hold on;
        plot(rmsf(:,find(~goodnucleitmp)),'color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('angles');ylabel('rmsf (\mum)');...
            axis([0 length(zrange) 0 0.4]);title('rmsf');
        set(gca,'Fontsize',30,'LineWidth',3);
        
        subplot(2,3,5)
        bins=0:0.0025:0.15;
        rmsf1=rmsf(:,find(goodnucleitmp));
        [counts]=hist(rmsf1(:),bins);
        cumcounts=cumsum(counts)./sum(counts);
        plot(bins,cumcounts,'linewidth',2,'color',...
            moviecolor((imovie),:));hold on;
        xlabel('rmsf(\mum)');ylabel('cumulative probability');...
            axis([0 0.1 0 1]);title('cum rmsf');
        set(gca,'Fontsize',30,'LineWidth',3);
        legend(legendids,'Interpreter','none','FontSize',6,'Location',...
            'southeast'); % Set Interpreter property 
            % to 'none'; the default for text fields is LaTeX, e.g., '_' in
            % string becomes an subscript in legend text. - JJ, 11/15/2017
        
        subplot(2,3,6)
        plot(dr2s(:,find(goodnucleitmp)),'color',...
            moviecolor((imovie),:));hold on;
        plot(dr2s(:,find(~goodnucleitmp)),'color',...
            moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('sum dr square over 6 (pixel^2)');...
            axis([0 length(zrange) 0 1]);title('outlier');
        set(gca,'Fontsize',30,'LineWidth',3);

        typermsf = [typermsf;rmsf1(:)];
    end
    
    % Save figures in designated folders.
    name = allnames{itype};
    mkdir([verifypath,'/',name]);
    print(f1,[verifypath,'/',name,'/_existflag'],'-dpng');
    print(f2,[verifypath,'/',name,'/_zpos'],'-dpng');
    print(f3,[verifypath,'/',name,'/_drift'],'-dpng');
    print(f4,[verifypath,'/',name,'/_rmsf'],'-dpng');
    print(f5,[verifypath,'/',name,'/_cumsumrmsf'],'-dpng');
    print(f7,[verifypath,'/',name,'/_sumdr2'],'-dpng');
    print(f8,[verifypath,'/',name],'-dpng');
         
    figure(f6)
    bins = 0:0.0025:0.15;
    [counts] = hist(typermsf(:),bins);
    cumcounts = cumsum(counts)./sum(counts);
    if ~isempty(strfind(name,'MBC'))
        plot(bins,cumcounts,'linestyle','-.','color',...
            allcolors(itype,:),'LineWidth',3)
        hold on;
    else
        plot(bins,cumcounts,'color',allcolors(itype,:),'LineWidth',3)
        hold on;
    end
    xlabel('RMSF (\mum)');ylabel('Cumulative Probability');...
        axis([0 0.1 0 1])
    set(gca,'Fontsize',30,'LineWidth',3);
    legend(allnames,'Interpreter','none','FontSize',10,'Location',...
        'southeast'); % Set Interpreter property 
        % to 'none'; the default for text fields is LaTeX, e.g., '_' in
        % string becomes an subscript in legend text. - JJ, 11/15/2017

end

FigureFormat(f6);
print(f6,[verifypath,'/s','_allcumsumrmsf'],'-dpng');
save(fullfile(rootpath,'goodnuclei.mat'),'gnmovienames','goodnuclei');
