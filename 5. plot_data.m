%% Plot data.

% clear all; close all;clc;
run setup_header_epe1_clr4_swi6_sm1.m;
rootpath='./';
load([rootpath,'/strainall_epe1_clr4_swi6_sm1.mat']);
resultpath=[rootpath,'result_epe1_clr4_swi6_sm1/'];
mkdir(resultpath)
num_names=length(nameAll);

% attention !!!!!!! size filter
sizeth=1000;%1.18;

%% Number of nuclei collected.

figure('Position',[0 0 1200 800]);

meandata = zeros(num_names,1);
sedata = zeros(num_names,1);
stddata = zeros(num_names,1);

% All nuclei (good + bad).
for itype = 1:length(nameAll2)
    meandata(itype) = length(strainall(itype).nuclei);
    bar(itype,meandata(itype),'FaceColor',colorAll(itype,:),'LineWidth',...
        2)
    alpha(.3)
    hold on
end

% Just the good nuclei...
meandata1 = meandata;
for itype=1:length(nameAll2)
    if ~isempty(strainall(itype).nuclei)
        meandata(itype)=sum([strainall(itype).nuclei.good]);
        bar(itype,meandata(itype),'FaceColor',colorAll(itype,:),...
            'LineWidth',2)
        hold on
    end
end

meandata2=meandata;

axis([0 num_names+1 0 500])
set(gca,'XTick',1:num_names,'XTicklabel',nameAll2(1:num_names),...
    'XTickLabelRotation',45);
ylabel('Nuclei analyzed');
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/number'],'-dpng');
print(gcf,[resultpath,'/number'],'-dsvg');
savefig(gcf,[resultpath,'/number']);

%% RMSF PDFs; RMSF CDFs plots for all strains/conditions.
% figure('Position',[0 0 1200 800]);

% Make folders for plots.
mkdir([resultpath,'/rmsf']);
mkdir([resultpath,'/cumrmsf']);
hold off;
% Set data range.
zrange=find(abs(points(:,3))<0.5);
xval=0.00:0.0025:0.1;
rmsfcounts=zeros(length(strainall),length(xval));
rmsfcumcounts=zeros(length(strainall),length(xval));

% Compile RMSF values.
for itype=1:length(strainall)
    str=strainall(itype);
    rmsftmp=[];
    for inuc=1:length(str.nuclei)
        if str.nuclei(inuc).good==1 & str.nuclei(inuc).size<=sizeth
            rmsftmp=[rmsftmp,str.nuclei(inuc).rmsf(zrange)];
        end
    end
    rmsfctmp=hist(rmsftmp,xval);
    rmsfcounts(itype,:)=rmsfctmp/sum(rmsfctmp);
    rmsfcumcounts(itype,:)=cumsum(rmsfcounts(itype,:));
end
    
%% Plot RMSF of each data set w/ & w/out MBC.

analyzename = 'rmsf_indiv_MBC';
mkdir([resultpath,'/',analyzename]);
for j=1:num_names/2
    rmsfchose=[j j+num_names/2];
    plot(xval'*ones(size(rmsfchose)),rmsfcounts(rmsfchose(1),:)',...
        'LineWidth',3,'Color',colorAll(rmsfchose(1),:))
    hold on
    plot(xval'*ones(size(rmsfchose)),rmsfcounts(rmsfchose(2),:)',...
        'LineWidth',3,'Color',colorAll(rmsfchose(1),:),'LineStyle',...
        '--')
    line([0.025 0.025],[0 0.3],'Color',[.3 .3 .3],'LineWidth',3,...
        'LineStyle',':')
    hold off
%     legend(nameAll2(rmsfchose))
    xlabel('RMSF (\mum)')
    ylabel('Probability')
    axis([0 0.08 0 0.25])
    set(gca,'Fontsize',30,'LineWidth',3)
    set(gca,'box','off')
    set(gca,'YTick',[0.05 0.15 0.25])
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(j))],...
        '-dpng')
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(j))],...
        '-dsvg')
    savefig([resultpath,'/',analyzename,'/',catname(nameAll(j))])
end
    
%% Plot cumulative RMSF of each data set w/ & w/out MBC.

analyzename = 'cumrmsf_indiv_MBC';
mkdir([resultpath,'/',analyzename]);
for j=1:num_names/2
    rmsfchose=[j j+num_names/2];
    plot(xval'*ones(size(rmsfchose)),rmsfcumcounts(rmsfchose(1),:)',...
        'LineWidth',3,'Color',colorAll(rmsfchose(1),:))
    hold on
    plot(xval'*ones(size(rmsfchose)),rmsfcumcounts(rmsfchose(2),:)',...
        'LineWidth',3,'Color',colorAll(rmsfchose(1),:),'LineStyle',...
        '--')
    line([0.025 0.025],[0 1],'Color',[.3 .3 .3],'LineWidth',3,...
        'LineStyle',':')
    hold off
%     legend(nameAll2(rmsfchose),'Location','southeast')
    xlabel('RMSF (\mum)')
    ylabel('Cumulative probability')
    axis([0 0.08 0 1])
    set(gca,'Fontsize',30,'LineWidth',3)
    set(gca,'box','off')
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(j))],...
        '-dpng')
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(j))],...
        '-dsvg')
    savefig([resultpath,'/',analyzename,'/',catname(nameAll(j))])
end
    
%% Plot RMSF for all datasets w/out MBC together.

figure('Position',[0 0 1000 500]);  
rmsfchose=1:num_names/2; % use this line with MBC data
% rmsfchose=1:num_names;
for ir=1:length(rmsfchose)
    plot(xval',rmsfcounts(rmsfchose(ir),:)','color',...
        colorAll(rmsfchose(ir),:),'LineWidth',5);
    hold on
end
line([0.025 0.025],[0 0.3],'Color',[.3 .3 .3],'LineWidth',3,...
        'LineStyle',':')
hold off
legend(nameAll2(rmsfchose),'Location','northeast','box','off')
xlabel('Root mean squared fluctuation (\mum)')
ylabel('Probability')
axis([0 0.08 0 0.16])
set(gca,'Fontsize',30,'LineWidth',2)
set(gca,'box','off')
print(gcf,[resultpath,'/rmsf/',catname(nameAll(rmsfchose))],'-dpng');
print(gcf,[resultpath,'/rmsf/',catname(nameAll(rmsfchose))],'-dsvg');
savefig([resultpath,'/rmsf/',catname(nameAll(rmsfchose))]);
    
%% Plot cumulative RMSF for all datasets w/out MBC together.

figure('Position',[0 0 800 500]);
rmsfchose=1:num_names/2; % use this line with MBC data
% rmsfchose=1:num_names;
for ir=1:length(rmsfchose)
    plot(xval',rmsfcumcounts(rmsfchose(ir),:)','color',...
        colorAll(rmsfchose(ir),:),'LineWidth',5)
    hold on
end
line([0.025 0.025],[0 1],'Color',[.3 .3 .3],'LineWidth',3,...
    'LineStyle',':')
hold off;
legend(nameAll2(rmsfchose),'Location','southeast','box','off')
xlabel('Root mean squared fluctuation (\mum)')
ylabel('Cumulative probability')
axis([0 0.08 0 1]);
set(gca,'Fontsize',30,'LineWidth',2)
set(gca,'box','off')
print(gcf,[resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))],...
    '-dpng');
print(gcf,[resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))],...
    '-dsvg');
savefig([resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))]);
    
%% Plot MBC-treated RMSFs together.

figure('Position',[0 0 800 500]);
rmsfchose = num_names/2+1:num_names;
% rmsfchose = [7,8,10,11,12];
for ir=1:length(rmsfchose)
    plot(xval',rmsfcounts(rmsfchose(ir),:)','color',...
        colorAll(rmsfchose(ir),:),'LineWidth',3,'LineStyle','-.')
    hold on
end
line([0.025 0.025],[0 1],'Color',[.3 .3 .3],'LineWidth',3,...
    'LineStyle',':')
hold off    
% legend(nameAll2(rmsfchose),'Location','northeast')
xlabel('RMSF (\mum)')
ylabel('Probability')
axis([0 0.06 0 .25]);
set(gca,'Fontsize',25,'LineWidth',2)
set(gca,'box','off')
print(gcf,[resultpath,'/rmsf/',catname(nameAll(rmsfchose))],'-dpng');
print(gcf,[resultpath,'/rmsf/',catname(nameAll(rmsfchose))],'-dsvg');
savefig([resultpath,'/rmsf/',catname(nameAll(rmsfchose))]);
  
%% Plot cumulative RMSFs MBC together   
    
figure('Position',[0 0 800 500]);
rmsfchose=num_names/2+1:num_names;
% rmsfchose = [7,8,10,11,12];
for ir=1:length(rmsfchose)
    plot(xval', rmsfcumcounts(rmsfchose(ir),:)','color',...
        colorAll(rmsfchose(ir),:),'LineStyle','-.','LineWidth',3)
    hold on
end
line([0.025 0.025],[0 1],'Color',[.3 .3 .3],'LineWidth',2,...
    'LineStyle',':')
hold off    
% legend(nameAll2(rmsfchose),'Location','southeast')
xlabel('RMSF (\mum)')
ylabel('Cumulative probability')
axis([0 0.1 0 1]);
set(gca,'Fontsize',25,'LineWidth',2)
set(gca,'box','off')
print(gcf,[resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))],...
    '-dpng');
print(gcf,[resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))],...
    '-dsvg');
savefig([resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))]);

%% Plot RMSFs of select groups together.

figure('Position',[0 0 1000 500]);
rmsfchose = [1,5,4]; % Strain indices (from nameAll) to plot.
for ir=1:length(rmsfchose)
    plot(xval',rmsfcounts(rmsfchose(ir),:)','color',...
        colorAll(rmsfchose(ir),:),'LineWidth',6)
    hold on
end
line([0.025 0.025],[0 1],'Color',[.3 .3 .3],'LineWidth',3,...
    'LineStyle',':')
hold off    
% legend(nameAll2(rmsfchose),'Location','northeastoutside')
xlabel('Root mean squared fluctuation (\mum)')
ylabel('Probability')
axis([0.02 0.08 0 0.16]);
set(gca,'Fontsize',40,'LineWidth',4)
set(gca,'box','off')
print(gcf,[resultpath,'/rmsf/',catname(nameAll(rmsfchose))],'-dpng');
print(gcf,[resultpath,'/rmsf/',catname(nameAll(rmsfchose))],'-dsvg');
savefig([resultpath,'/rmsf/',catname(nameAll(rmsfchose))]);

%% Plot cumulative RMSF of select groups together.

figure('Position',[0 0 800 500]);
rmsfchose = [1,5,4]; % Strain indices (from nameAll) to plot.
for ir=1:length(rmsfchose)
    plot(xval',rmsfcumcounts(rmsfchose(ir),:)','color',...
        colorAll(rmsfchose(ir),:),'LineWidth',6)
    hold on
end
line([0.025 0.025],[0 1],'Color',[.3 .3 .3],'LineWidth',3,...
    'LineStyle',':')
hold off;
legend(nameAll2(rmsfchose),'Location','southeast','box','off')
xlabel('Root mean squared fluctuation (\mum)')
ylabel('Cumulative probability')
axis([0.02 0.08 0 1]);
set(gca,'Fontsize',40,'LineWidth',4)
set(gca,'box','off')
print(gcf,[resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))],...
    '-dpng');
print(gcf,[resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))],...
    '-dsvg');
savefig([resultpath,'/cumrmsf/',catname(nameAll(rmsfchose))]);

%% Plot RMSFs of MBC and non-MBC data together.

analyzename = 'rmsf_combined';
mkdir([resultpath,'/',analyzename]);

figure('Position',[0 0 1000 500]);

% Plot first the lower limit of flutuation resolution.
line([0.025 0.025],[0 1],'Color',[.5 .5 .5],'LineWidth',4,...
    'LineStyle',':')
hold on

% Plot MBC data as dashed lines.
rmsfchose2 = [6,8];
for ir = 1:length(rmsfchose2)
    plot(xval',rmsfcounts(rmsfchose2(ir),:)','Color',...
        colorAll(rmsfchose2(ir),:),'LineWidth',5,'LineStyle','-.')
end

% Plot non-MBC data as solid lines.
rmsfchose = [1,3];
for ir=1:length(rmsfchose)
    plot(xval',rmsfcounts(rmsfchose(ir),:)','Color',...
        colorAll(rmsfchose(ir),:),'LineWidth',5)
end
hold off 

% legend(nameAll2(rmsfchose),'Location','northeastoutside')
xlabel('Root mean squared fluctuation (\mum)')
ylabel('Probability')
axis([0.02 0.08 0 0.26]);
set(gca,'Fontsize',40,'LineWidth',4)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(rmsfchose))],...
    '-dpng');
print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(rmsfchose))],...
    '-dsvg');
savefig([resultpath,'/',analyzename,'/',catname(nameAll(rmsfchose))]);

%% Plot cumulative RMSFs of MBC and non-MBC data together.

analyzename = 'cumrmsf_combined';
mkdir([resultpath,'/',analyzename]);

figure('Position',[0 0 800 500]);

% Plot MBC data first as dashed lines.
rmsfchose2 = [6,8];
for ir = 1:length(rmsfchose2)
    plot(xval', rmsfcumcounts(rmsfchose2(ir),:)','color',...
        colorAll(rmsfchose2(ir),:),'LineStyle','-.','LineWidth',5)
    hold on
end

% Plot non-MBC data on top as solid lines.
rmsfchose = [1,3];
for ir=1:length(rmsfchose)
    plot(xval',rmsfcumcounts(rmsfchose(ir),:)','color',...
        colorAll(rmsfchose(ir),:),'LineWidth',5)
end

line([0.025 0.025],[0 1],'Color',[.5 .5 .5],'LineWidth',4,...
    'LineStyle',':')
hold off 

% legend(nameAll2(rmsfchose),'Location','southeast','box','off')
legend('','','','MBC-treated','Location','southeast','box','off')
xlabel('Root mean squared fluctuation (\mum)')
ylabel('Cumulative probability')
axis([0.02 0.08 0 1]);
set(gca,'Fontsize',40,'LineWidth',4)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(rmsfchose))],...
    '-dpng');
print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(rmsfchose))],...
    '-dsvg');
savefig([resultpath,'/',analyzename,'/',catname(nameAll(rmsfchose))]);
    
%% RMSF vs. nuclear angle plot.

analyzename='rmsf_angle';
mkdir([resultpath,'/',analyzename]);
hold off
zrange=find(abs(points(:,3))<0.5);
nbins=40;
dx=180/nbins;
xval=dx:dx:180;
rmsfangle=cell(length(xval),num_names);
%     meanrmsfangle=zeros(length(angles),num_names);
%     stdrmsfangle=zeros(length(angles),num_names);
for itype=1:length(strainall)
    str=strainall(itype);
    for inuc=1:length(str.nuclei)
        if str.nuclei(inuc).good==1
            celltheta=str.nuclei(inuc).orientation;
            nuctheta=acos(points(zrange,1)*cos(celltheta) + ...
                points(zrange,2)*sin(celltheta))/pi*180;
            rmsf=str.nuclei(inuc).rmsf(zrange);
            for ith=1:length(nuctheta)
                intangle=floor(nuctheta(ith)/dx)+1;
                if intangle==nbins+1
                    intangle=1;
                end
                rmsfangle{intangle,itype}=[rmsfangle{intangle,itype},...
                    rmsf(ith)];
            end
        end
    end
end
meanrmsfangle=cellfun(@mean,rmsfangle);
numrmsfangle=cellfun(@length,rmsfangle);
stdrmsfangle=cellfun(@std,rmsfangle);
sermsfangle=stdrmsfangle./sqrt(numrmsfangle);
    
%% Plot RMSF vs. angle, for each data set w/ & w/out MBC, w/ small error
% bars (s.e.m.).

for j=1:num_names/2
    rmsfchose=[j j+num_names/2];
    % Plot RMSF w/out MBC.
    plot(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose(1)),...
        'color',colorAll(rmsfchose(1),:),'LineWidth',3)
    hold on
    % Error bars.
    errorbar(xval'*ones(size(rmsfchose(1))),meanrmsfangle(:,rmsfchose(1)),...
        sermsfangle(:,rmsfchose(1)),'color',colorAll(rmsfchose(1),:),...
        'LineWidth',3)
    % Plot RMSF w/ MBC.
    plot(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose(2)),...
        'color',colorAll(rmsfchose(1),:),'LineWidth',1)
    errorbar(xval'*ones(size(rmsfchose(2))),meanrmsfangle(:,rmsfchose(2)),...
        sermsfangle(:,rmsfchose(2)),'color',colorAll(rmsfchose(1),:),...
        'LineWidth',1.5)
    line([0 180],[0.025 0.025],'Color',[.3 .3 .3],'LineWidth',2,...
        'LineStyle',':')
%     legend(nameAll2(rmsfchose))
    xlabel('Nuclear radial angle (degrees)')
    ylabel('RMSF (\mum)')
    axis([min(xval) max(xval) 0.024 0.05]);
    set(gca,'XTick',[0 45 90 135 180])
    set(gca,'Fontsize',30,'LineWidth',3)
    set(gca,'box','off')
    print(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))],'-depsc');
    savefig(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))]);
    hold off;
end

%% Plot RMSF vs. angle, for each data set w/out MBC, w/ band of standard 
% deviation.

analyzename = 'rmsf_angle_std';
mkdir([resultpath,'/',analyzename]);
figure('Position',[1 1 590 450])
for j=1:num_names/2
    rmsfchose=[j j+num_names/2];
    % Band of +/- standard deviation above/below mean values.
    f1 = fill([xval'*ones(size(rmsfchose(1)));...
        flipud(xval'*ones(size(rmsfchose(1))))],...
        [meanrmsfangle(:,rmsfchose(1)) - stdrmsfangle(:,rmsfchose(1));...
        flipud(meanrmsfangle(:,rmsfchose(1)) + ...
        stdrmsfangle(:,rmsfchose(1)))],colorAll(rmsfchose(1),:),...
        'LineStyle','none');
%     alpha(f1,.5)
    hold on
    % Plot mean RMSF w/out MBC.
    plot(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose(1)),...
        'color',colorAll(rmsfchose(1),:),'LineWidth',3)
    line([0 180],[0.025 0.025],'Color',[.5 .5 .5],'LineWidth',4,...
        'LineStyle',':')
    xlabel({'Nuclear radial angle','(degrees)'})
    ylabel({'Root mean squared','fluctuation (\mum)'})
    axis([min(xval) max(xval) 0.024 0.06]);
    set(gca,'XTick',[0 45 90 135 180])
    set(gca,'Fontsize',34,'LineWidth',3)
    set(gca,'box','off');
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(j))],...
        '-dpng');
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(j))],...
        '-dsvg');
    savefig(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(j))]);
    hold off;
end

%% Plot RMSF vs. angle, for each data set w/ & w/out MBC, w/ band of 
% standard deviation.

analyzename2 = 'rmsf_angle_std_MBC';
mkdir([resultpath,'/',analyzename2]);
for j=1:num_names/2
    rmsfchose=[j j+num_names/2];
    % Band of +/- standard deviation above/below mean values.
    f1 = fill([xval'*ones(size(rmsfchose(1)));...
        flipud(xval'*ones(size(rmsfchose(1))))],...
        [meanrmsfangle(:,rmsfchose(1)) - stdrmsfangle(:,rmsfchose(1));...
        flipud(meanrmsfangle(:,rmsfchose(1)) + ...
        stdrmsfangle(:,rmsfchose(1)))],colorAll(rmsfchose(1),:),...
        'LineStyle','none');
    alpha(f1,.5)
    hold on
    % Plot mean RMSF w/out MBC.
    plot(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose(1)),...
        'color',colorAll(rmsfchose(1),:),'LineWidth',3)
    % Band of +/- standard deviation of MBC data above/below mean values.
    f2 = fill([xval'*ones(size(rmsfchose(2)));...
        flipud(xval'*ones(size(rmsfchose(2))))],...
        [meanrmsfangle(:,rmsfchose(2)) - stdrmsfangle(:,rmsfchose(2));...
        flipud(meanrmsfangle(:,rmsfchose(2)) + ...
        stdrmsfangle(:,rmsfchose(2)))],colorAll(rmsfchose(1),:),...
        'LineStyle','none');
    alpha(f2,.25)
    % Plot mean RMSF w/ MBC.
    plot(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose(2)),...
        'color',colorAll(rmsfchose(1),:),'LineWidth',1.5)
    line([0 180],[0.025 0.025],'Color',[.3 .3 .3],'LineWidth',2,...
        'LineStyle',':')
%     legend(nameAll2(rmsfchose))
    xlabel('Nuclear radial angle (degrees)')
    ylabel('RMSF (\mum)')
    axis([min(xval) max(xval) 0.024 0.06]);
    set(gca,'XTick',[0 45 90 135 180])
    set(gca,'Fontsize',34,'LineWidth',3)
    set(gca,'box','off');
    print(gcf,[resultpath,'/',analyzename2,'/',...
        catname(nameAll(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'/',analyzename2,'/',...
        catname(nameAll(rmsfchose))],'-dsvg');
    savefig(gcf,[resultpath,'/',analyzename2,'/',...
        catname(nameAll(rmsfchose))]);
    hold off;
end

%% Plot RMSF vs. angle as polar plot.

analyzename3 = 'rmsf_angle_polar_plot';
mkdir([resultpath,'/',analyzename3]);
% Convert xval (degrees) to radians.
xval_rad = deg2rad(xval);

% For plotting lower limit of fluctuation resolution.
% rad = linspace(0,pi,100);
% r = 0.025;

for j=1:num_names/2
    rmsfchose=[j j+num_names/2];
    % Plot mean RMSF vs angle.
    polarplot(xval_rad'*ones(size(rmsfchose)),...
        meanrmsfangle(:,rmsfchose(1)),'Color',colorAll(rmsfchose(1),:),...
        'LineWidth',5)
    hold on
    % Plot standard deviation below mean as a dashed line.
    polarplot(xval_rad'*ones(size(rmsfchose(1))),...
        meanrmsfangle(:,rmsfchose(1)) - stdrmsfangle(:,rmsfchose(1)),...
        'Color',colorAll(rmsfchose(1),:),'LineWidth',1)
    % Plot standard deviation above mean as a dashed line.
    polarplot(xval_rad'*ones(size(rmsfchose(1))),...
        meanrmsfangle(:,rmsfchose(1)) + stdrmsfangle(:,rmsfchose(1)),...
        'Color',colorAll(rmsfchose(1),:),'LineWidth',1)
    % Plot lower limit of fluctuation resolution.
%     polarplot(rad,r + zeros(size(rad)),'Color',[.3 .3 .3],'LineWidth',2,...
%         'LineStyle',':')
    thetalim([0 180])
    rlim([0 0.07])
    set(gca,'RTickLabelRotation',45)
    set(gca,'Fontsize',30,'LineWidth',3)
    set(gca,'box','off')  
    print(gcf,[resultpath,'/',analyzename3,'/',...
        catname(nameAll(j))],'-dpng');
    print(gcf,[resultpath,'/',analyzename3,'/',...
        catname(nameAll(j))],'-dsvg');
    savefig(gcf,[resultpath,'/',analyzename3,'/',...
        catname(nameAll(j))]);
    hold off;
end

%% Plot RMSF vs. angle all together w/out MBC.

% figure('Position',[0 0 1200 800]);  
rmsfchose=1:num_names/2; % Use this line with MBC data
% rmsfchose=1:num_names;
for i = 1:length(rmsfchose)
    plot(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,rmsfchose(i)),...
        'color',colorAll(rmsfchose(i),:),'LineWidth',3)
    hold on
    errorbar(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,...
        rmsfchose(i)),sermsfangle(:,rmsfchose(i)),'color',...
        colorAll(rmsfchose(i),:),'LineWidth',3,'HandleVisibility','off');
end
% legend(nameAll2(rmsfchose),'Location','northeastoutside');
line([0 180],[0.025 0.025],'Color',[.3 .3 .3],'LineWidth',3,...
    'LineStyle',':')
xlabel('Nuclear radial angle (degrees)')
ylabel('RMSF (\mum)')
axis([min(xval) max(xval) 0.024 0.05]);
set(gca,'XTick',[0 45 90 135 180])
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose))],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose))],'-dsvg');
savefig(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose))]);

%% RMSF vs. angle w/ MBC, altogether.
 
% figure('Position',[0 0 1200 800]);  
rmsfchose=num_names/2+1:num_names;
for i = 1:length(rmsfchose)
    plot(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,rmsfchose(i)),...
        'color',colorAll(rmsfchose(i),:),'LineWidth',1)
    hold on
    errorbar(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,...
        rmsfchose(i)),sermsfangle(:,rmsfchose(i)),'color',...
        colorAll(rmsfchose(i),:),'LineWidth',1,'HandleVisibility','off');
end
%     legend(nameAll2(rmsfchose),'Location','northeastoutside');
    line([0 180],[0.025 0.025],'Color',[.3 .3 .3],'LineWidth',2,...
        'LineStyle',':')
    xlabel('Nuclear radial angle (degrees)')
    ylabel('RMSF (\mum)')
    axis([min(xval) max(xval) 0.02 0.055]);
    set(gca,'Fontsize',30,'LineWidth',3)
    print(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))],'-depsc');
    savefig(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))]);
    
%% Plot RMSF vs. angle of select groups together of non-MBC data.

analyzename = 'rmsf_angle';
mkdir([resultpath,'/',analyzename]);
% figure('Position',[0 0 1200 800]);  
rmsfchose = [1,2,3,4]; % Strain indices (from nameAll) to plot.
for i = 1:length(rmsfchose)
    plot(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,rmsfchose(i)),...
        'color',colorAll(rmsfchose(i),:),'LineWidth',3)
    hold on
    errorbar(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,...
        rmsfchose(i)),sermsfangle(:,rmsfchose(i)),'color',...
        colorAll(rmsfchose(i),:),'LineWidth',3,'HandleVisibility','off');
end
% legend(nameAll2(rmsfchose),'Location','northeastoutside');
line([0 180],[0.025 0.025],'Color',[.3 .3 .3],'LineWidth',3,...
    'LineStyle',':')
xlabel('Nuclear radial angle (degrees)')
ylabel('RMSF (\mum)')
axis([min(xval) max(xval) 0.024 0.05]);
set(gca,'XTick',[0 45 90 135 180])
set(gca,'YTick',[0.03 0.04 0.05])
set(gca,'Fontsize',25,'LineWidth',2)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose))],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose))],'-dsvg');
savefig(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose))]);

%% Plot RMSF vs. angle of select groups together of MBC data.

analyzename = 'rmsf_angle_MBC_zoom';
mkdir([resultpath,'/',analyzename]);
% figure('Position',[0 0 1200 800]);  
rmsfchose2 = [7,8,10,11,12];
for i = 1:length(rmsfchose2)
    plot(xval'*ones(size(rmsfchose2(i))),meanrmsfangle(:,rmsfchose2(i)),...
        'color',colorAll(rmsfchose2(i),:),'LineWidth',1)
    hold on
    errorbar(xval'*ones(size(rmsfchose2(i))),meanrmsfangle(:,...
        rmsfchose2(i)),sermsfangle(:,rmsfchose2(i)),'color',...
        colorAll(rmsfchose2(i),:),'LineWidth',1,'HandleVisibility','off');
end
% legend(nameAll2(rmsfchose),'Location','northeastoutside');
line([0 180],[0.025 0.025],'Color',[.3 .3 .3],'LineWidth',3,...
    'LineStyle',':')
xlabel('Nuclear radial angle (degrees)')
ylabel('RMSF (\mum)')
axis([min(xval) max(xval) 0.024 0.038]);
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose2))],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose2))],'-dsvg');
savefig(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose2))]);

%% Plot RMSF vs. angle of select groups w/ & w/out MBC combined.

analyzename = 'rmsf_angle_combined';
mkdir([resultpath,'/',analyzename]);

% Plot non-MBC data first.
rmsfchose = [1,3]; % Strain indices (from nameAll) to plot.
for i = 1:length(rmsfchose)
    plot(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,rmsfchose(i)),...
        'color',colorAll(rmsfchose(i),:),'LineWidth',5)
    hold on
    errorbar(xval'*ones(size(rmsfchose(i))),meanrmsfangle(:,...
        rmsfchose(i)),sermsfangle(:,rmsfchose(i)),'color',...
        colorAll(rmsfchose(i),:),'LineWidth',4,'HandleVisibility','off');
end

% Plot MBC data.
rmsfchose2 = [6,8];
for i = 1:length(rmsfchose2)
    plot(xval'*ones(size(rmsfchose2(i))),meanrmsfangle(:,rmsfchose2(i)),...
        'color',colorAll(rmsfchose2(i),:),'LineWidth',1.5)
    hold on
    errorbar(xval'*ones(size(rmsfchose2(i))),meanrmsfangle(:,...
        rmsfchose2(i)),sermsfangle(:,rmsfchose2(i)),'color',...
        colorAll(rmsfchose2(i),:),'LineWidth',1.5,'HandleVisibility','off');
end
% legend(nameAll2(rmsfchose),'Location','northeastoutside');
line([0 180],[0.025 0.025],'Color',[.5 .5 .5],'LineWidth',4,...
    'LineStyle',':')
xlabel({'Nuclear radial angle','(degrees)'})
ylabel({'Average root mean','squared fluctuation (\mum)'})
axis([min(xval) max(xval) 0.024 0.05]);
set(gca,'XTick',[0 45 90 135 180])
set(gca,'YTick',[0.03 0.04 0.05])
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose2))],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose2))],'-dsvg');
savefig(gcf,[resultpath,'/',analyzename,'/',...
    catname(nameAll(rmsfchose2))]);

%% Fluctuations subtracted: RMSF vs. angle plot.

analyzename='fluc_sub_rmsf_angle';
mkdir([resultpath,'/',analyzename]);
hold off;
zrange=find(abs(points(:,3))<0.5);
inzrange=abs(points(:,3))<0.5;
nbins=40;
dx=180/nbins;
xval=dx:dx:180;
rmsfangle=cell(length(xval),num_names);
%     meanrmsfangle=zeros(length(angles),num_names);
%     stdrmsfangle=zeros(length(angles),num_names);
for itype=1:length(strainall)
    str=strainall(itype);
    for inuc=1:length(str.nuclei)
        if str.nuclei(inuc).good==1
            celltheta=str.nuclei(inuc).orientation;
            goodindtmp=inzrange & str.nuclei(inuc).nonfluczone;
            nuctheta=acos(points(goodindtmp,1)*cos(celltheta) + ...
                points(goodindtmp,2)*sin(celltheta))/pi*180;
            rmsf=str.nuclei(inuc).rmsf(goodindtmp);
            for ith=1:length(nuctheta)
                intangle=floor(nuctheta(ith)/dx)+1;
                if intangle==nbins+1
                    intangle=1;
                end
                rmsfangle{intangle,itype}=[rmsfangle{intangle,itype},...
                    rmsf(ith)];
            end
        end
    end
end
meanrmsfangle=cellfun(@mean,rmsfangle);
numrmsfangle=cellfun(@length,rmsfangle);
stdrmsfangle=cellfun(@std,rmsfangle);
sermsfangle=stdrmsfangle./sqrt(numrmsfangle);

%% Plot fluctuation-subtracted RMSF vs. angle, individually w/ & w/out MBC.
  
for j=1:num_names/2
    rmsfchose=[j j+num_names/2];
    % Plot RMSF w/out MBC.
    plot(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose(1)),...
        'color',colorAll(rmsfchose(1),:),'LineWidth',3)
    hold on
    % Error bars.
    errorbar(xval'*ones(size(rmsfchose(1))),meanrmsfangle(:,rmsfchose(1)),...
        sermsfangle(:,rmsfchose(1)),'color',colorAll(rmsfchose(1),:),...
        'LineWidth',3)
    % Plot RMSF w/ MBC.
    plot(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose(2)),...
        'color',colorAll(rmsfchose(1),:),'LineWidth',1)
    errorbar(xval'*ones(size(rmsfchose(2))),meanrmsfangle(:,rmsfchose(2)),...
        sermsfangle(:,rmsfchose(2)),'color',colorAll(rmsfchose(1),:),...
        'LineWidth',1)
    line([0 180],[0.025 0.025],'Color',[.3 .3 .3],'LineWidth',2,...
        'LineStyle',':')
    legend(nameAll2(rmsfchose))
    xlabel('Nuclear radial angle (degrees)')
    ylabel('RMSF (\mum)')
    axis([min(xval) max(xval) 0.02 0.055]);
    set(gca,'Fontsize',30,'LineWidth',3)
    print(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))],'-depsc');
    savefig(gcf,[resultpath,'/',analyzename,'/',...
        catname(nameAll(rmsfchose))]);
    hold off;
end
   
%%
    rmsfchose=1:num_names/2;
    plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));
    errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
    legend(nameAll2(rmsfchose));
    xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0.02 0.06]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
    
    rmsfchose=num_names/2+1:num_names;
    plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));
    errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
    legend(nameAll2(rmsfchose));
    xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0.02 0.06]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
    end

%% Number of fluctuations collected.

figure('Position',[0 0 1200 800]);
fparam = struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,...
    'durationub',200,'anglelb',30,'angleub',80);

meandata1 = zeros(num_names,1);
meandata2 = zeros(num_names,1);
sedata = zeros(num_names,1);
stddata = zeros(num_names,1);

% All fluctuations.
for itype = 1:length(nameAll2)
    datatmp = (strainall(itype).flucs);
    gooddata = [strainall(itype).flucs.good] & ...
        [strainall(itype).flucs.inplane];
    meandata1(itype) = sum(gooddata);
    bar(itype,meandata1(itype),'FaceColor',colorAll(itype,:),...
        'LineWidth',2)
    alpha(.3)
    hold on
end

% MT-induced fluctuations.
    for itype = 1:length(nameAll2)
        if ~isempty(strainall(itype).nuclei)
            datatmp = (strainall(itype).flucs);
            gooddata = [strainall(itype).flucs.good] & ...
                [strainall(itype).flucs.inplane];
            [s_ind, ds_ind] = SelectMtFluc(datatmp(gooddata),fparam);
            meandata2(itype) = sum(s_ind);
            bar(itype,meandata2(itype),'facecolor',colorAll(itype,:),...
                'LineWidth',2)
            hold on
        end
    end
    
axis([0 num_names+1 0 3000])
set(gca,'XTick',1:num_names,'XTicklabel',nameAll2(1:num_names),...
    'XTickLabelRotation',45)
ylabel('Fluctuations detected')
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/flucnumber'],'-dpng')
print(gcf,[resultpath,'/flucnumber'],'-dsvg')
savefig(gcf,[resultpath,'/flucnumber'])

%% Mean fluctuation size (standard deviation above surface mean) surface 
% kymograph.

figure('Position',[0 0 1200 800]);
analyzename = 'mean_fluc_map';
mkdir([resultpath,'/',analyzename]);
hold off;
%     stdimgs=zeros(19,);
for itype = 1:length(strainall)
    str = strainall(itype);
    stdtmp = 0;
    numgood = 0;
    for inuc = 1:length(str.nuclei)
        if str.nuclei(inuc).good == 1
            stdtmp = stdtmp + str.nuclei(inuc).stdimg;
            numgood = numgood + 1;
        end
    end
    stdtmp = stdtmp/numgood;
    imagesc(stdtmp,[0.025 0.05])
%     imagesc(stdtmp,[0.03 0.05])
    h = colorbar;
%     CUDmap = [0,114,178   % blue
%               7,119,183
%               14,125,187
%               21,130,192
%               29,136,196
%               36,141,201
%               43,147,205
%               50,152,210
%               57,158,214
%               64,163,219
%               72,169,224
%               79,174,228
%               86,180,233
%               80,178,224
%               74,177,216
%               67,175,208
%               61,174,199
%               55,172,191
%               49,170,182
%               43,169,174
%               37,167,165
%               31,166,157
%               24,164,149
%               18,163,140
%               12,161,132
%               6,159,123
%               0,158,115
%               22,164,110
%               44,171,106
%               65,177,101
%               87,183,97
%               109,190,93
%               131,196,88
%               152,202,84
%               174,209,79
%               196,215,75
%               218,222,70
%               240,228,66
%               239,222,60
%               238,216,55
%               237,211,49
%               237,205,44
%               236,199,38
%               235,193,33
%               234,188,27
%               233,182,22
%               232,176,16
%               232,170,11
%               231,165,5
%               230,159,0
%               228,156,12
%               226,154,24
%               224,151,36
%               222,148,48
%               221,145,59
%               219,143,71
%               217,140,83
%               215,137,95
%               213,134,107
%               211,132,119
%               209,129,131
%               208,126,143
%               206,124,155
%               204,121,167]/255;
%     colormap(CUDmap)
    
    colormap(colormap_CUD)
    ylabel(h, 'Mean fluctuation height (\mum)')
    axis image;
    title([nameAll2(itype)]);
    xlabel({'Longitudinal angle','(degrees)'})
    ylabel({'Latitudinal angle','(radians)'})
    set(gca,'XTick',(0:16:64));
    xticklabels({'0','90','180','270','360'})
    set(gca,'YTick',([2,10,18]));
    yticklabels({'\pi/4','0','-\pi/4'})
    set(gca,'Fontsize',60,'LineWidth',3)
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(itype))],...
        '-dpng');
    print(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(itype))],...
        '-dsvg');
    savefig(gcf,[resultpath,'/',analyzename,'/',catname(nameAll(itype))]);
end

%% Fluctuation size (standard deviation) surface kymograph (individual).
% Plots surface kymograph (heat map) of individual nuclei.


analyzename = 'fluc_map_indiv_select';
mkdir([resultpath,'/',analyzename]);
hold off

% for itype = 1:length(strainall)

%     for inuc = 1:length(str.nuclei)
%         if str.nuclei(inuc).good == 1

%%

itype = 9; % strain type from strainall
inuc = 122; % nucleus number from strainall.nuclei
    str = strainall(itype);
    stdtmp = 0;
    figure('Position',[0 0 1200 800]);
    stdtmp = str.nuclei(inuc).stdimg;
    imagesc(stdtmp,[0.02 0.085])
    h = colorbar;    
    colormap(colormap_bluered)
    ylabel(h, 'Fluctuation height (\mum)')
    axis image;
    title([nameAll2(itype)]);
    xlabel({'Longitudinal angle','(degrees)'})
    ylabel({'Latitudinal angle','(degrees)'})
    set(gca,'XTick',(0:16:64));
    xticklabels({'0','90','180','270','360'})
    set(gca,'YTick',([2,10,18]));
    yticklabels({'45','0','-45'})
    set(gca,'Fontsize',60,'LineWidth',3)           
    print(gcf,[resultpath,'/',analyzename,'/',...
        char(strcat(nameAll(itype),'_',num2str(inuc)))],...
        '-dpng');
    print(gcf,[resultpath,'/',analyzename,'/',...
        char(strcat(nameAll(itype),'_',num2str(inuc)))],...
        '-dsvg');
    savefig(gcf,[resultpath,'/',analyzename,'/',...
        char(strcat(nameAll(itype),'_',num2str(inuc)))]);


%% Nuclear radius size.

figure('Position',[0 0 1200 800]);
analyzename = 'radius';
yval = 0.5:0.05:1.8;
mkdir([resultpath,'/',analyzename]);
% datalength = num_names/2; % use with MBC data
datalength = num_names;
histdata = zeros(datalength,length(yval));
cumhistdata = zeros(datalength,length(yval));
meandata = zeros(datalength,1);
sedata = zeros(datalength,1);
stddata = zeros(datalength,1);
radii = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = [strainall(itype).nuclei.size];
        goodtmp = find([strainall(itype).nuclei.good] & ...
            [strainall(itype).nuclei.size]<sizeth);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata(itype) = mean(datatmp);
        stddata(itype) = std(datatmp);
        sedata(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'FaceColor',colorAll(itype,:),...
            'LineWidth',3)
        hold on
        radii(itype).radius = datatmp';
    end
end
e = errorbar(1:datalength,meandata,stddata,'k');
e.LineWidth = 3;
e.LineStyle = 'none';
axis([0 datalength+1 min(yval) 1.6]);
set(gca,'XTick',1:datalength,'XTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45);
ylabel('Nuclear radius (\mum)');
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
hold off;
    
%%  Cell Width.

figure('Position',[0 0 1200 800]);

analyzename='cell_width';
yval=2:0.1:5;
mkdir([resultpath,'/',analyzename]);
% datalength = num_names/2; % Use with MBC data
datalength = num_names;
histdata = zeros(datalength,length(yval));
cumhistdata = zeros(datalength,length(yval));
meandata = zeros(datalength,1);
sedata = zeros(datalength,1);
stddata = zeros(datalength,1);
cellwidth = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = [strainall(itype).nuclei.cellwidth];
        goodtmp = find([strainall(itype).nuclei.good]);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata(itype) = mean(datatmp);
        stddata(itype) = std(datatmp);
        sedata(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'FaceColor',colorAll(itype,:),...
            'LineWidth',3)
        hold on
        cellwidth(itype).cellwidth = datatmp';
    end
end
e = errorbar(1:datalength,meandata,stddata,'k');
e.LineWidth = 3;
e.LineStyle = 'none';
axis([0 datalength+1 min(yval) max(yval)])
set(gca,'XTick',1:datalength,'XTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45);
ylabel('Cell width (\mum)')
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
    
%%  Cell length.

figure('Position',[0 0 1200 800]);

analyzename = 'cell_length';
yval = 6:0.2:12;
mkdir([resultpath,'/',analyzename]);
% datalength = num_names/2; % Use with MBC data
datalength = num_names;
histdata = zeros(datalength,length(yval));
cumhistdata = zeros(datalength,length(yval));
meandata = zeros(datalength,1);
sedata = zeros(datalength,1);
stddata = zeros(datalength,1);
celllength = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = [strainall(itype).nuclei.celllength];
        goodtmp = find([strainall(itype).nuclei.good]);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata(itype) = mean(datatmp);
        stddata(itype) = std(datatmp);
        sedata(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'facecolor',colorAll(itype,:),...
            'LineWidth',3)
        hold on
        celllength(itype).celllength = datatmp';
    end
end
e = errorbar(1:datalength,meandata,stddata,'k');
e.LineWidth = 3;
e.LineStyle = 'none';
axis([0 datalength+1 min(yval) max(yval)+2]);
set(gca,'xTick',1:datalength,'xTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45)
ylabel('Cell length (\mum)')
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);

%%  Cell volume.

figure('Position',[0 0 1200 800]);

analyzename = 'cell_volume';
yval=30:2:150;
mkdir([resultpath,'/',analyzename]);
% datalength = num_names/2; % Use with MBC data
datalength = num_names;
histdata = zeros(datalength,length(yval));
cumhistdata = zeros(datalength,length(yval));
meandata = zeros(datalength,1);
sedata = zeros(datalength,1);
stddata = zeros(datalength,1);
cellvolume = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = [strainall(itype).nuclei.cellvolume];
        goodtmp = find([strainall(itype).nuclei.good]);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata(itype) = mean(datatmp);
        stddata(itype) = std(datatmp);
        sedata(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'facecolor',colorAll(itype,:),...
            'LineWidth',3)
        hold on
        cellvolume(itype).cellvolume = datatmp';
    end
end

e = errorbar(1:datalength,meandata,stddata,'k');
e.LineWidth = 3;
e.LineStyle = 'none';
axis([0 datalength+1 min(yval) max(yval)+50]);
set(gca,'XTick',1:datalength,'XTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45);
ylabel('Cell volume (\mum^3)');
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);

%%  Nuclear volume.

figure('Position',[0 0 1200 800]);

analyzename = 'nucleus_volume';
yval=0:0.05:12;
mkdir([resultpath,'/',analyzename]);
% datalength = num_names/2; % Use with MBC data
datalength = num_names;
histdata = zeros(datalength,length(yval));
cumhistdata = zeros(datalength,length(yval));
meandata = zeros(datalength,1);
sedata = zeros(datalength,1);
stddata = zeros(datalength,1);
nucleusvolume = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = [strainall(itype).nuclei.nucleusvolume];
        goodtmp = find([strainall(itype).nuclei.good]);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata(itype) = mean(datatmp);
        stddata(itype) = std(datatmp);
        sedata(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'facecolor',colorAll(itype,:),...
            'LineWidth',3)
        hold on
        nucleusvolume(itype).nucleusvolume = datatmp';
    end
end
e = errorbar(1:datalength,meandata,stddata,'k');
e.LineWidth = 3;
e.LineStyle = 'none';
axis([0 datalength+1 min(yval) max(yval) ]);
set(gca,'XTick',1:datalength,'XTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45);
ylabel('Nuclear volume (\mum^3)');
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
    
%% Nuclear to cell volume ratio.

figure('Position',[0 0 1200 800]);

analyzename = 'nc_ratio';
yval = 0:0.005:0.15;
mkdir([resultpath,'/',analyzename]);
% datalength = num_names/2; % Use with MBC data.
datalength = num_names;
histdata = zeros(datalength,length(yval));
cumhistdata = zeros(datalength,length(yval));
meandata = zeros(datalength,1);
sedata = zeros(datalength,1);
stddata = zeros(datalength,1);
ncratio = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = [strainall(itype).nuclei.nucleusvolume]./...
            [strainall(itype).nuclei.cellvolume];
        goodtmp = find([strainall(itype).nuclei.good]);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata(itype) = mean(datatmp);
        stddata(itype) = std(datatmp);
        sedata(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'facecolor',colorAll(itype,:),...
            'LineWidth',3)
        hold on
        ncratio(itype).ncratio = datatmp';
    end
end
e = errorbar(1:datalength,meandata,stddata,'k');
e.LineWidth = 3;
e.LineStyle = 'none';
axis([0 datalength+1 min(yval) max(yval)])
set(gca,'XTick',1:datalength,'XTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45)
ylabel('N/C ratio')
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
    
%% Fluctuation location and max height angle.

fparam=struct('maxheightlb',0,'meanheightlb',0,'durationlb',0,...
    'durationub',inf,'anglelb',0,'angleub',90);

analyzename='max_height_angle';
mkdir([resultpath,'/',analyzename]);
datalength = num_names;
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        fluctmp = straintmp.flucs;
        fluctmp = fluctmp([fluctmp.inplane]==1);
        [ s_ind, ds_ind ] = SelectMtFluc( fluctmp,fparam);
        fluctmp = fluctmp(find(s_ind));
        datatmp1 = acos(cos([fluctmp.avglat]).*cos([fluctmp.avglong]))/...
            pi*180;
        datatmp2 = [fluctmp.maxheight];
        goodtmp = find([fluctmp.good]);
        datatmp1 = datatmp1(goodtmp);
        datatmp2 = datatmp2(goodtmp);
        xval = 0:9:180;
        histdata = hist(datatmp1,xval);
        histdata = histdata/sum(histdata);
        fig = figure('Position',[0 0 1200 800]);
        left_color = [0 0 0];
        right_color = colorAll(itype,:);
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        yyaxis right
        bh = bar(xval,histdata,'FaceColor',colorAll(itype,:),'LineWidth',2)
        alpha(.6)
        xlim([0 180])
        ylim([0 0.15])
        ylabel('Fluctuation probability','Color',colorAll(itype,:))
        yyaxis left
        scatter(datatmp1,datatmp2,30,'MarkerFaceColor','k',...
            'MarkerEdgeColor','k')
        ylim([0.05 0.3])
        xlabel('Nuclear angle (degrees)')
        ylabel('Maximum fluctuation height (\mum)','Color','k')
%         title('Fluctuation distribution');
        legend(nameAll2(itype),'Location','northeast','box','off')
        set(gca,'Fontsize',30,'LineWidth',3)
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-dpng');
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-depsc');
    end
end

%% Fluctuation location and mean height angle w/ MT-fluctuations.

fparam=struct('maxheightlb',0,'meanheightlb',0,'durationlb',20,...
    'durationub',200,'anglelb',0,'angleub',90);
fparam2=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,...
    'durationub',200,'anglelb',30,'angleub',80);

analyzename = 'mean_height_angle';
% titlestr = {'Fluctuation distribution blue: non MT fluc, red: MT fluc'};
mkdir([resultpath,'/',analyzename]);
datalength=num_names/2;
for itype=1:datalength
    straintmp=strainall(itype);
    if ~isempty(straintmp.nuclei)
       
        % Plot MT fluctuations first.
        fluctmp=straintmp.flucs;
        [ s_ind, ds_ind ] = SelectMtFluc( fluctmp,fparam2);
        fluctmp=fluctmp(find(s_ind));
        datatmp1=acos(cos([fluctmp.avglat]).*cos([fluctmp.avglong]))/pi*180;
        datatmp2=[fluctmp.meanheight];
        goodtmp=find([fluctmp.good]);
        datatmp1=datatmp1(goodtmp);
        datatmp2=datatmp2(goodtmp);
        xval1=0:num_names/2+1:180;
        histdata=hist(datatmp1,xval1);
        histdata=histdata/sum(histdata);
        
        % Plot non-MT fluctuations.
        fluctmp=straintmp.flucs;
        fluctmp=fluctmp([fluctmp.inplane]==1);
        [ s_ind, ds_ind ] = SelectMtFluc( fluctmp,fparam);
        fluctmp=fluctmp(find(s_ind));
        datatmp3=acos(cos([fluctmp.avglat]).*cos([fluctmp.avglong]))/pi*180;
        datatmp4=[fluctmp.meanheight];
        goodtmp=find([fluctmp.good]);
        datatmp3=datatmp3(goodtmp);
        datatmp4=datatmp4(goodtmp);
        xval2=0:9:180;
        histdata2=hist(datatmp3,xval2);
        histdata2=histdata2/sum(histdata2);
        
        fig = figure('Position',[0 0 1200 800]);
        left_color = [0 0 0];
        right_color = colorAll(itype,:);
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        yyaxis right
        bh1 = bar(xval2,histdata2,'FaceColor',colorAll(itype,:),...
            'FaceAlpha',0.2,'LineWidth',2)
        bh2 = bar(xval1,histdata,'FaceColor',colorAll(itype,:),...
            'LineWidth',2)
        ylim([0 0.15])
        ylabel('Fluctuation probability')

        yyaxis left
        scatter(datatmp1,datatmp2,100,'MarkerFaceColor',colorAll(itype,:),...
            'MarkerEdgeColor',colorAll(itype,:))
        scatter(datatmp3,datatmp4,30,'MarkerFaceColor','k',...
            'MarkerEdgeColor','k')
        ylim([0.04 0.3])
        xlim([0 180])
        xlabel('Nuclear angle (degrees)')
        ylabel('Mean fluctuation height (\mum)')
%         title(titlestr);
        legend(nameAll2(itype),'Location','northeast','box','off')
        set(gca,'Fontsize',30,'LineWidth',3)
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-dpng');
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-depsc');
        savefig(gcf,[resultpath,'/',analyzename,'/',straintmp.name]);
    end
end

%% Global fluctuation location.

figure('Position',[0 0 1200 800]);
fparam = struct('maxheightlb',0,'meanheightlb',0,'durationlb',0,...
    'durationub',inf,'anglelb',0,'angleub',90);
analyzename = 'fluclocation';
mkdir([resultpath,'/',analyzename]);
datalist = [1,2,3,4,5];
% Non-MBC data:
for itype = 1:length(datalist)
    straintmp = strainall(datalist(itype));
    if ~isempty(straintmp.flucs)
        fluctmp = straintmp.flucs;
        %             fluctmp=fluctmp([fluctmp.inplane]==1);
        [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
        datatmp1 = [fluctmp.avglong]/pi*180;
        datatmp2 = [fluctmp.avglat]/pi*180;
        scatter(datatmp1(s_ind),datatmp2(s_ind),46,'MarkerFaceColor',...
            colorAll(itype,:),'MarkerEdgeColor',...
            colorAll(datalist(itype),:))
        hold on
        scatter(datatmp1(ds_ind),datatmp2(ds_ind),46,'MarkerFaceColor',...
            colorAll(datalist(itype),:))
        line([0 360],[45 45],'Color',[.4 .4 .4],'LineWidth',4)
        line([0 360],[-45 -45],'Color',[.4 .4 .4],'LineWidth',4)
        hold off
        xlabel('Longitude (degrees)')
        ylabel('Latitude (degrees)')
        axis([0 180 -100 100])
        axis equal
%         title('fluctuations location');
        legend(nameAll2(itype),'Location','northeast')
        set(gca,'Fontsize',40,'LineWidth',3)
        set(gca,'box','off')
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-dpng');
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-dsvg');
        savefig(gcf,[resultpath,'/',analyzename,'/',straintmp.name]);
    end
end

% MBC data:
datalist = [6,7,8,9,10];
for itype = 1:length(datalist)
    straintmp = strainall(datalist(itype));
    if ~isempty(straintmp.flucs)
        fluctmp = straintmp.flucs;
        %             fluctmp=fluctmp([fluctmp.inplane]==1);
        [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
        datatmp1 = [fluctmp.avglong]/pi*180;
        datatmp2 = [fluctmp.avglat]/pi*180;
        scatter(datatmp1(s_ind),datatmp2(s_ind),'MarkerFaceColor',...
            'none','MarkerEdgeColor',colorAll(datalist(itype),:),...
            'LineWidth',1.5)
        hold on
        scatter(datatmp1(ds_ind),datatmp2(ds_ind),'MarkerFaceColor',...
            'none')
        line([0 360],[45 45],'Color',[.4 .4 .4],'LineWidth',4)
        line([0 360],[-45 -45],'Color',[.4 .4 .4],'LineWidth',4)
        hold off
        xlabel('Longitude (degrees)')
        ylabel('Latitude (degrees)')
        axis([0 180 -100 100])
        axis equal
%         title('fluctuations location');
        legend(nameAll2(datalist(itype)),'Location','northeast')
        set(gca,'Fontsize',40,'LineWidth',3)
        set(gca,'box','off')
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-dpng');
        print(gcf,[resultpath,'/',analyzename,'/',straintmp.name],'-dsvg');
        savefig(gcf,[resultpath,'/',analyzename,'/',straintmp.name]);
    end
end

%% Number of fluctuations per cell per minute.

figure('Position',[0 0 1200 800]);
fparam1=struct('maxheightlb',0,'meanheightlb',0,'durationlb',20,...
    'durationub',200,'anglelb',0,'angleub',90);
fparam2=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,...
    'durationub',200,'anglelb',30,'angleub',80);

analyzename='fluc_per_min';
% titlestr='All fluctuations (dark and light) & MT-induced fluctuations (dark)';
yval=(0:1:4)/0.5/250*60;
mkdir([resultpath,'/',analyzename]);
datalength=num_names;
histdata=zeros(datalength,length(yval));
cumhistdata=zeros(datalength,length(yval));
meandata1=zeros(datalength,1);
sedata1=zeros(datalength,1);
stddata1=zeros(datalength,1);
meandata2=zeros(datalength,1);
sedata2=zeros(datalength,1);
stddata2=zeros(datalength,1);
   
% MT fluctuations.
mtflucnum = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = zeros(1,length(straintmp.nuclei));
        for inucs = 1:length(straintmp.nuclei)
            fluctmp = straintmp.nuclei(inucs).flucs;
            if ~isempty(fluctmp)
                [ s_ind, ds_ind ] = ...
                    SelectMtFluc(fluctmp([fluctmp.inplane]==1),fparam2);
            else
                s_ind = [];
            end
            datatmp(inucs) = sum(s_ind)/0.5/250*60;
        end
        goodtmp = find([strainall(itype).nuclei.good]);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata1(itype) = mean(datatmp);
        stddata1(itype) = std(datatmp);
        sedata1(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata1(itype),'facecolor',colorAll(itype,:),...
            'LineWidth',3)
        hold on
        mtflucnum(itype).mtflucnum = datatmp'; % Save matrix of values.
    end
end

% Remove zeros.
for i = 1:length(mtflucnum)
    mtflucnum(i).mtflucnum = ...
        mtflucnum(i).mtflucnum(mtflucnum(i).mtflucnum > 0);
end

e1 = errorbar(1:datalength,meandata1,stddata1,'k')
e1.LineWidth = 3;
e1.LineStyle = 'none';

% non-MT fluctuations.
nonmtflucnum = [];
for itype = 1:datalength
    straintmp = strainall(itype);
    if ~isempty(straintmp.nuclei)
        datatmp = zeros(1,length(straintmp.nuclei));
        for inucs = 1:length(straintmp.nuclei)
            fluctmp = straintmp.nuclei(inucs).flucs;
            if ~isempty(fluctmp)
                [ s_ind, ds_ind ] = ...
                    SelectMtFluc( fluctmp([fluctmp.inplane]==1),fparam1);
            else
                s_ind = [];
            end
            datatmp(inucs) = sum(s_ind)/0.5/250*60;
        end
        goodtmp = find([strainall(itype).nuclei.good]);
        datatmp = datatmp(goodtmp);
        histdata(itype,:) = hist(datatmp,yval);
        histdata(itype,:) = histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:) = cumsum(histdata(itype,:));
        meandata2(itype) = mean(datatmp);
        stddata2(itype) = std(datatmp);
        sedata2(itype) = std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata2(itype),'facecolor',colorAll(itype,:),...
            'FaceAlpha',0.3,'LineWidth',3)
        hold on
        nonmtflucnum(itype).nonmtflucnum = datatmp'; % Save matrix of 
                                                     % values.
    end
end

% Remove zeros.
for i = 1:length(nonmtflucnum)
    nonmtflucnum(i).nonmtflucnum = ...
        nonmtflucnum(i).nonmtflucnum(nonmtflucnum(i).nonmtflucnum > 0);
end

e2 = errorbar(1:datalength,meandata2,stddata2,'k')
e2.LineWidth = 3;
e2.LineStyle = 'none';
% title(titlestr);
axis([0 datalength+1 min(yval) max(yval)+6])
set(gca,'XTick',1:datalength,'XTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45)
ylabel('Fluctuation frequency (1/s)')
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
    
%% MT rise fall (original)

figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,'durationub',200,...
    'anglelb',30,'angleub',80);

analyzename='MTFlucTime';
xname='time (s)';
xval=0:10:150;
mkdir([resultpath,'/',analyzename]);

datalength=num_names;
histdata1=zeros(datalength,length(xval));
cumhistdata1=zeros(datalength,length(xval));
histdata2=zeros(datalength,length(xval));
cumhistdata2=zeros(datalength,length(xval));
    
%     datalength=num_names/2; set to only use non MBC data
meandata1=zeros(datalength,1);
sedata1=zeros(datalength,1);
stddata1=zeros(datalength,1);
meandata2=zeros(datalength,1);
sedata2=zeros(datalength,1);
stddata2=zeros(datalength,1);
    
for itype=1:datalength
    straintmp=strainall(itype);
    if ~isempty(straintmp.nuclei)
        fluctmp=straintmp.flucs;
        fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
        [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
        datatmp=[fluctmp.risetime];
        goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&s_ind);
        datatmp=datatmp(goodtmp);
        histdata1(itype,:)=hist(datatmp,xval);
        histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
        cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
        meandata1(itype)=mean(datatmp);
        stddata1(itype)=std(datatmp);
        sedata1(itype)=std(datatmp)./sqrt(length(datatmp));
        B=barh(datalength+1-itype,meandata1(itype),'facecolor',colorAll(itype,:));hold on;
        set(get(B,'child'),'facea',.3);
    end
end
%     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
herrorbar(meandata1,datalength:-1:1,sedata1,'k');hold on;
    
for itype=1:datalength
    straintmp=strainall(itype);
    if ~isempty(straintmp.nuclei)
        fluctmp=straintmp.flucs;
        fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
        [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
        datatmp=[fluctmp.falltime];
        goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&s_ind);
        datatmp=datatmp(goodtmp);
        histdata2(itype,:)=hist(datatmp,xval);
        histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
        cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
        meandata2(itype)=mean(datatmp);
        stddata2(itype)=std(datatmp);
        sedata2(itype)=std(datatmp)./sqrt(length(datatmp));
        B=barh(datalength+1-itype,meandata2(itype),'facecolor',colorAll(itype,:));hold on;
    end
end
%     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
herrorbar(meandata2,datalength:-1:1,sedata2,'k');hold on;

title('MT induced fluc: rise time (light+dark) and fall time (dark)');
axis([min(xval) max(xval) 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
xlabel(xname);
FigureFormat(gcf)
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-deps');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
hold off;
    
%     datalength=num_names/2;
%     histdata=zeros(datalength,length(xval));
%     cumhistdata=zeros(datalength,length(xval));
%     meandata=zeros(datalength,1);
%     sedata=zeros(datalength,1);
%     stddata=zeros(datalength,1);
%     
%     for itype=1:datalength
%         straintmp=strainall(itype+num_names/2);
%         if ~isempty(straintmp.nuclei)
%             [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
%             datatmp=[straintmp.flucs.risetime];
%             goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&s_ind);
%             datatmp=datatmp(goodtmp);
%             histdata1(itype+num_names/2,:)=hist(datatmp,xval);
%             histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
%             cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
%             meandata(itype)=mean(datatmp);
%             stddata(itype)=std(datatmp);
%             sedata(itype)=std(datatmp)./sqrt(length(datatmp));
%             B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
%             set(get(B,'child'),'facea',.3);
%         end
%     end
%     %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
%     herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
%     
%     for itype=1:datalength
%         straintmp=strainall(itype+num_names/2);
%         if ~isempty(straintmp.nuclei)
%             [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
%             datatmp=[straintmp.flucs.falltime];
%             goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&s_ind);
%             datatmp=datatmp(goodtmp);
%             histdata2(itype+num_names/2,:)=hist(datatmp,xval);
%             histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
%             cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
%             meandata(itype)=mean(datatmp);
%             stddata(itype)=std(datatmp);
%             sedata(itype)=std(datatmp)./sqrt(length(datatmp));
%             B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
%         end
%     end
%     %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
%     herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
%     title('MT induced fluc: rise time (light+dark) and fall time (dark)');
%     axis([min(xval) max(xval) 0 datalength+1]);
%     set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
%     xlabel(xname);
%     FigureFormat(gcf)
%     print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-dpng');
%     print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-deps');
%     savefig(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC']);
%     hold off;
    

%% non MT rise fall
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,'durationub',200,...
    'anglelb',30,'angleub',80);
for i=1
    analyzename='NonMTFlucTime';
    xname='time (s)';
    xval=0:10:150;
    mkdir([resultpath,'\',analyzename]);
    
    datalength=num_names;
    histdata1=zeros(datalength,length(xval));
    cumhistdata1=zeros(datalength,length(xval));
    histdata2=zeros(datalength,length(xval));
    cumhistdata2=zeros(datalength,length(xval));
    
    datalength=num_names/2;
    meandata1=zeros(datalength,1);
    sedata1=zeros(datalength,1);
    stddata1=zeros(datalength,1);
    meandata2=zeros(datalength,1);
    sedata2=zeros(datalength,1);
    stddata2=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.risetime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype,:)=hist(datatmp,xval);
            histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
            cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
            meandata1(itype)=mean(datatmp);
            stddata1(itype)=std(datatmp);
            sedata1(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata1(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata1,datalength:-1:1,sedata1,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.falltime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype,:)=hist(datatmp,xval);
            histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
            cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
            meandata2(itype)=mean(datatmp);
            stddata2(itype)=std(datatmp);
            sedata2(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata2(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata2,datalength:-1:1,sedata2,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    datalength=num_names/2;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.risetime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype+num_names/2,:)=hist(datatmp,xval);
            histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
            cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.falltime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype+num_names/2,:)=hist(datatmp,xval);
            histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
            cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC']);
    hold off;
    
end
%% MT rise fall 2

figure('Position',[0 0 1200 800]);
% fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',40,...
%     'durationub',100,'anglelb',30,'angleub',80);
fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,...
    'durationub',200,'anglelb',30,'angleub',80);

analyzename='MTFlucTime2';
yname='Time (s)';
yval=0:10:150;
mkdir([resultpath,'/',analyzename]);

% datalength=num_names;
% histdata1=zeros(datalength,length(xval));
% cumhistdata1=zeros(datalength,length(xval));
% histdata2=zeros(datalength,length(xval));
% cumhistdata2=zeros(datalength,length(xval));

datalength=num_names/2; % Use when MBC data is included.
% datalength=num_names;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
MTrisetime = [];
for itype=1:datalength
    straintmp=strainall(itype);
    if ~isempty(straintmp.nuclei)
        fluctmp=straintmp.flucs;
        fluctmp=fluctmp([fluctmp.inplane]==1);
        [s_ind, ds_ind] = SelectMtFluc(fluctmp,fparam);
        datatmp=[fluctmp.risetime];
        goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&s_ind);
        datatmp=datatmp(goodtmp);
        histdata1(itype,:)=hist(datatmp,yval);
        histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
        cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
        meandata(itype)=mean(datatmp);
        stddata(itype)=std(datatmp);
        sedata(itype)=std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'facecolor',colorAll(itype,:),...
            'LineWidth',3)
        alpha(.3)
        MTrisetime(itype).MTrisetime = datatmp';
        hold on
    end
end
e1 = errorbar(1:datalength,meandata,sedata,'k');
e1.LineWidth = 3;
e1.LineStyle = 'none';

MTfalltime = [];
for itype=1:datalength
    straintmp=strainall(itype);
    if ~isempty(straintmp.nuclei)
        fluctmp=straintmp.flucs;
        fluctmp=fluctmp([fluctmp.inplane]==1);
        [s_ind, ds_ind] = SelectMtFluc(fluctmp,fparam);
        datatmp=[fluctmp.falltime];
        goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&s_ind);
        datatmp=datatmp(goodtmp);
        histdata2(itype,:)=hist(datatmp,yval);
        histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
        cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
        meandata(itype)=mean(datatmp);
        stddata(itype)=std(datatmp);
        sedata(itype)=std(datatmp)./sqrt(length(datatmp));
        bar(itype,meandata(itype),'facecolor',colorAll(itype,:),...
            'LineWidth',3)
        MTfalltime(itype).MTfalltime = datatmp';
        hold on
    end
end
e2 = errorbar(1:datalength,meandata,sedata,'k')
e2.LineWidth = 3;
e2.LineStyle = 'none';

title('Rise time (light), fall time (dark)')
axis([0 datalength+1 min(yval) 80])
set(gca,'XTick',1:datalength,'XTicklabel',nameAll2(1:datalength),...
    'XTickLabelRotation',45)
ylabel(yname)
set(gca,'Fontsize',30,'LineWidth',3)
set(gca,'box','off')
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-depsc');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
    
%     datalength=num_names/2;
%     histdata=zeros(datalength,length(xval));
%     cumhistdata=zeros(datalength,length(xval));
%     meandata=zeros(datalength,1);
%     sedata=zeros(datalength,1);
%     stddata=zeros(datalength,1);
%     
%     for itype=1:datalength
%         straintmp=strainall(itype+num_names/2);
%         if ~isempty(straintmp.nuclei)
%             [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
%             datatmp=[straintmp.flucs.risetime];
%             goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&s_ind);
%             datatmp=datatmp(goodtmp);
%             histdata1(itype+num_names/2,:)=hist(datatmp,xval);
%             histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
%             cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
%             meandata(itype)=mean(datatmp);
%             stddata(itype)=std(datatmp);
%             sedata(itype)=std(datatmp)./sqrt(length(datatmp));
%             B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
% %             set(get(B,'child'),'facea',.3);
%             alpha(.3)
%         end
%     end
%     %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
%     herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
%     
%     for itype=1:datalength
%         straintmp=strainall(itype+num_names/2);
%         if ~isempty(straintmp.nuclei)
%             [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
%             datatmp=[straintmp.flucs.falltime];
%             goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&s_ind);
%             datatmp=datatmp(goodtmp);
%             histdata2(itype+num_names/2,:)=hist(datatmp,xval);
%             histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
%             cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
%             meandata(itype)=mean(datatmp);
%             stddata(itype)=std(datatmp);
%             sedata(itype)=std(datatmp)./sqrt(length(datatmp));
%             B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
%         end
%     end
% %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
% herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
% 
% title('MT-induced Fluctuations: Rise time (light + dark); Fall time (dark)');
% axis([min(xval) max(xval) 0 datalength+1]);
% set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
% xlabel(xname);
% FigureFormat(gcf)
% print(gcf,[resultpath,'/',analyzename,'/',analyzename,'MBC'],'-dpng');
% print(gcf,[resultpath,'/',analyzename,'/',analyzename,'MBC'],'-deps');
% savefig(gcf,[resultpath,'/',analyzename,'/',analyzename,'MBC']);  

%% all rise fall
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0,'durationlb',20,'durationub',200,...
    'anglelb',0,'angleub',90);

analyzename='AllFlucTime';
xname='time (s)';
xval=0:10:150;
mkdir([resultpath,'/',analyzename]);

% datalength=num_names;
% histdata1=zeros(datalength,length(xval));
% cumhistdata1=zeros(datalength,length(xval));
% histdata2=zeros(datalength,length(xval));
% cumhistdata2=zeros(datalength,length(xval));

datalength=num_names/2; % Use w/ MBC data
% datalength=num_names;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);

for itype=1:datalength
    straintmp=strainall(itype);
    if ~isempty(straintmp.nuclei)
        fluctmp=straintmp.flucs;
        fluctmp=fluctmp([fluctmp.inplane]==1);
        [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
        datatmp=[fluctmp.risetime];
        goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&s_ind);
        datatmp=datatmp(goodtmp);
        histdata1(itype,:)=hist(datatmp,xval);
        histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
        cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
        meandata(itype)=mean(datatmp);
        stddata(itype)=std(datatmp);
        sedata(itype)=std(datatmp)./sqrt(length(datatmp));
        B=barh(datalength+1-itype,meandata(itype),'facecolor',...
            colorAll(itype,:))
        alpha(0.3)
        hold on
    end
end

herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;

for itype=1:datalength
    straintmp=strainall(itype);
    if ~isempty(straintmp.nuclei)
        fluctmp=straintmp.flucs;
        fluctmp=fluctmp([fluctmp.inplane]==1);
        [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
        datatmp=[fluctmp.risetime];
        goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&s_ind);
        datatmp=datatmp(goodtmp);
        histdata2(itype,:)=hist(datatmp,xval);
        histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
        cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
        meandata(itype)=mean(datatmp);
        stddata(itype)=std(datatmp);
        sedata(itype)=std(datatmp)./sqrt(length(datatmp));
        B=barh(datalength+1-itype,meandata(itype),'facecolor',...
            colorAll(itype,:))
        hold on
    end
end
herrorbar(meandata,datalength:-1:1,sedata,'k')
hold on
title('MT induced fluc: rise time (light+dark) and fall time (dark)');
axis([min(xval) max(xval) 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
xlabel(xname);
set(gca,'Fontsize',30,'LineWidth',3)
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-dpng');
print(gcf,[resultpath,'/',analyzename,'/',analyzename],'-deps');
savefig(gcf,[resultpath,'/',analyzename,'/',analyzename]);
hold off;

% datalength=num_names/2;
% histdata=zeros(datalength,length(xval));
% cumhistdata=zeros(datalength,length(xval));
% meandata=zeros(datalength,1);
% sedata=zeros(datalength,1);
% stddata=zeros(datalength,1);
% 
% for itype=1:datalength
%     straintmp=strainall(itype+num_names/2);
%     if ~isempty(straintmp.nuclei)
%         [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
%         datatmp=[straintmp.flucs.risetime];
%         goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&s_ind);
%         datatmp=datatmp(goodtmp);
%         histdata1(itype+num_names/2,:)=hist(datatmp,xval);
%         histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
%         cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
%         meandata(itype)=mean(datatmp);
%         stddata(itype)=std(datatmp);
%         sedata(itype)=std(datatmp)./sqrt(length(datatmp));
%         B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
%         set(get(B,'child'),'facea',.3);
%     end
% end
% %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
% herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
% 
% for itype=1:datalength
%     straintmp=strainall(itype+num_names/2);
%     if ~isempty(straintmp.nuclei)
%         [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
%         datatmp=[straintmp.flucs.falltime];
%         goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&s_ind);
%         datatmp=datatmp(goodtmp);
%         histdata2(itype+num_names/2,:)=hist(datatmp,xval);
%         histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
%         cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
%         meandata(itype)=mean(datatmp);
%         stddata(itype)=std(datatmp);
%         sedata(itype)=std(datatmp)./sqrt(length(datatmp));
%         B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
%     end
% end
% herrorbar(meandata,datalength:-1:1,sedata,'k')
% hold on
% title('MT induced fluc: rise time (light+dark) and fall time (dark)');
% axis([min(xval) max(xval) 0 datalength+1]);
% set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
% xlabel(xname);
% FigureFormat(gcf)
% print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-dpng');
% print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-deps');
% savefig(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC']);
% hold off;
    
%%
% fparam=struct('maxheightlb',0.1,'meanheightlb',0,'durationlb',20,'durationub',200,...
%     'anglelb',0,'angleub',90);
