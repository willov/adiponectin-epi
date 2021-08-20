
%% Setup things
% close all

cd('Scripts')

modelName='adr_endo';
fileDir='../Results/PPL/191220';
load('../Results/opt(28.3531).mat')
bestparam=optParam;

toEstimateOn={'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP','highCa_ATP','EPI_ATP','CL_Ca'};
toValidateOn={'CL_ATP'};

[model, expData, estimation, validation, dgf]=Init(modelName, toEstimateOn, toValidateOn);

design=matlab.lang.makeValidName({'Ca+ATP';'Ca+noATP';'noCa+ATP';'noCa+noATP';'highCa+ATP';...
    'EPI+ATP';'CL+ATP';'CL+Ca'});
%   Ca,           ATP,         cAMP,   Epi, CL
stimulus=[...
    0.0015         3           0.1     0    0;   %'Ca+ATP'
    0.0015         0           0.1     0    0;   %'Ca+noATP'
    0              3           0.1     0    0;   %'noCa+ATP'
    0              0           0.1     0    0;   %'noCa+noATP'
    0.015          3           0.1     0    0;   %'highCa+ATP'
    0              3           0       5    0;   %'EPI+ATP'
    0              3           0       0    1;   %'CL+ATP'
    0.0015         3           0       0    1];   %'CL+Ca'

experiments=table(stimulus(:,1:3), stimulus(:,4:end), 'VariableNames',{'Pipette','Agonist'},'RowNames',design);

%% Select what to plot
choice = input('What to Plot?  \n 1. Estimation (Fig 3)  \n 2. Validation (Fig 4)   \n 3. Prediction (Fig 5)   \n 4. Insights (Fig 6)  \n Choice: ');

%% Load parameters
if choice ~=5
    costFun=@(param) CostFunction(param , model, expData, estimation, dgf);
    bestCost=costFun(optParam);
    sprintf('Cost: %.1f, dgf: %i, Chi2: %.1f', bestCost, dgf, chi2inv(0.95, dgf))
    
    if ~exist([fileDir '/allParams.mat'],'file')
        files=dir([fileDir '/**/*.mat']);
        nFiles=size(files,1);
        params=nan(nFiles,length(bestparam));
        for p=1:nFiles
            load([files(p).folder '/' files(p).name]);
            if round(costFun(optParam)*1e4)<=round((bestCost+chi2inv(0.95,1))*1e4)
                params(p,:)=optParam;
            elseif round(costFun(optParam))<=round((bestCost+chi2inv(0.95,1)))
                params(p,:)=optParam;
                disp('This parameter set is above the limit, likely due to numerical differences in the simulation')
            else
                disp('ERROR!!')
            end
            fprintf('Done loading %i of %i files\n',p,nFiles)
        end
        params(~any(~isnan(params),2),:)=[]; % Removes dummy parameter sets
        if ~exist(fileDir,'dir')
            mkdir(fileDir)
        end
    else
        load([fileDir '/allParams.mat'])
    end
    
    time=unique(expData{'Time',:});
end
%% Get uncertainty and plot
if choice == 1 % Estimation data
    [ boundry, bestSim] = MinMax(model, params, bestparam, expData, experiments(toEstimateOn,:));
    Plot_Exocytosis_Uncertainty (bestSim, boundry, expData)
elseif choice == 2 % Validation
    [ boundry, bestSim] = MinMax(model, params, bestparam, expData, experiments(toValidateOn,:));
    Plot_Exocytosis_Uncertainty (bestSim, boundry, expData)
    
elseif choice == 3 % Prediction
    exoExperiments=experiments({'EPI_ATP', 'CL_ATP', 'CL_Ca'},:);
    [boundry, bestSim] = MinMax(model, params, bestparam, expData, exoExperiments);
    boundry30 = MinMax(model, params, bestparam, expData, exoExperiments, 12*60, 0.3);
    boundry60 = MinMax(model, params, bestparam, expData, exoExperiments, 12*60, 0.6);
    for k = 1:length(exoExperiments.Properties.RowNames)
        subplot(3,2,k+1), Plot_Subplot(boundry, bestSim('Time',:), exoExperiments.Properties.RowNames(k), [], [])
        subplot(3,2,k+1), hold on, Plot_Subplot(boundry30, bestSim('Time',:), exoExperiments.Properties.RowNames(k), [], [66, 194, 245]/255)
        subplot(3,2,k+1), hold on, Plot_Subplot(boundry60, bestSim('Time',:), exoExperiments.Properties.RowNames(k), [], [0 0 1])
        axis([0 12 0 24.66])
        fprintf('----%s----\n',exoExperiments.Properties.RowNames{k})
        fprintf('Peak, Normal: %.2f\n',max(boundry{k,'Max'}));
        fprintf('Peak, Inhib30: %.2f\n',max(boundry30{k,'Max'}));
        fprintf('Peak, Inhib60: %.2f\n',max(boundry60{k,'Max'}));
        fprintf('Peak down, 30: %.2f\n',(1-max(boundry30{k,'Max'})/max(boundry{k,'Max'}))*100);
        fprintf('Peak down, 60: %.2f\n',(1-max(boundry60{k,'Max'})/max(boundry{k,'Max'}))*100);
    end
    
    load('fig5a-data.mat')
    subplot(3,2,1)
    b=bar([1 2 3 4],siRNA{'Mean',:});
    hold on
    
    b.FaceColor = 'flat';
    b.CData=[207 207 207; 207 207 207; 77 77 77; 77 77 77]/255;
    errorbar([1 2 3 4],siRNA{'Mean',:},siRNA{'SEM',:},'marker','none','linestyle','none','color','k');
    set(gca, 'xticklabel', {'Control','EPI', 'Control','EPI'}, 'XTickLabelRotation', 0)
    set(gca, 'ytick',[0 50 100],'box','off')
    ylabel({'Adiponectin release'; '(μg/g protein)'})
    
    load('fig5e-data.mat')
    subplot(3,2,k+2)
    b=bar([1 2 3 4],release_iono{'Mean',:});
    hold on
    
    b.FaceColor = 'flat';
    b.CData=[207 207 207; 207 207 207; 207 207 207; 207 207 207]/255;
    errorbar([1 2 3 4],release_iono{'Mean',:},release_iono{'SEM',:},'marker','none','linestyle','none','color','k');
    set(gca, 'xticklabel', {'Control','FSK/IBMX', 'Iono.','FSK/IBMX/Iono.'}, 'XTickLabelRotation', 25)
    set(gca, 'ytick',[1 3 5],'box','off')
    ylabel({'Adiponectin release'; '(fold increase)'})
    
    set(gcf,'Position',[554 358 906 765]) %[1000 676 560 662]
elseif choice==4 % Insight
    
    mechanismExperiments=experiments({'EPI_ATP', 'CL_ATP', 'CL_Ca', 'noCa_ATP'},:);
    [boundry, bestSim] = MinMax(model, params, bestparam, expData, mechanismExperiments, 30*60);
    boundry{'Time','Max'}=boundry{'Time','Max'}/60;
    figure(4)
    m=size(boundry.MaxStates,3)+1;
    n=size(boundry.MaxStates,1)-1;
    k=1;
    
    boundry.MaxStates(abs(boundry.MaxStates)<1e-16)=0;
    boundry.MinStates(abs(boundry.MinStates)<1e-16)=0;
    boundry.MaxStates(:,:,2)=  boundry.MaxStates(:,:,2)*1e3;
    boundry.MinStates(:,:,2)=  boundry.MinStates(:,:,2)*1e3;
    
    ylabels = {{'\DeltaC/\Deltat','(fF/s)'}, {'Active \beta-receptors','(a.u)'}, {'cAMP increase','over basal (\muM)'}, {'Releasable pool','(a.u)'}, {'Adiponectin release','(fold over CL at t=15)'}};
    for j=1:n
        subplot(m,n,k)
        fill([boundry{'Time','Max'} fliplr(boundry{'Time','Max'})],...
            [boundry.Max(j,:), fliplr(boundry.Min(j,:))],[0.95,0.65,0],'EdgeColor',[0.95,0.65,0],'linewidth',2)
        title(boundry.Properties.RowNames{j},'Interpreter','none')
        k=k+1;
        axis([-.5 30.5 -1 30])
        box off
        xlabel('Time (Min)')
    end
    subplot(m,n,1)
    ylabel(ylabels{1})
    
    for i = 1:m-1
        subplot(m,n,k)
        ylabel(ylabels{i+1})
        hold on
        for j=1:n
            subplot(m,n,k)
            fill([boundry{'Time','Max'} fliplr(boundry{'Time','Max'})],...
                [boundry.MaxStates(j,:,i), fliplr(boundry.MinStates(j,:,i))],[ .8 .8 .8],'EdgeColor','[ .8 .8 .8]','linewidth',2)
            
            if i==1
                ylim([-5 100])
            elseif i==2
                ylim([-1 250])
            elseif i==3
                ylim([-0.1 3.5])
            elseif i==4 && j~=3
                ylim([-0.1 1.4])
            end
            xlim([-.5 30.5])
            
            box off
            xlabel('Time (Min)')
            k=k+1;
        end
    end
    
    CLAdi= table([6.09e-5     6.62e-5] , [0 1.24e-5], [15 30],'variablenames',{'MeanValues','SEMValues','Time'});
    CLAdi{:,1:end-1}=CLAdi{:,1:end-1}/CLAdi.MeanValues(CLAdi.Time==15);
    subplot(m,n,18)
    hold on
    errorbar(CLAdi.Time,CLAdi.MeanValues, CLAdi.SEMValues,'ko','linewidth',2,'MarkerFaceColor','auto')
    
    u = boundry.MaxStates(1:2,end,2)';
    l = boundry.MinStates(1:2,end,2)';
    mu = [2.24e-02, 1.95e-03 ]*1e3;
    s = [6.72e-04, 5.10e-04 ]*1e3;
    subplot(m,n,9)
    hold on
    errorbar(30, mu(1),s(1),'ko','linewidth',2,'MarkerFaceColor','auto')
    subplot(m,n,10)
    hold on
    errorbar(30, mu(2),s(2),'ko','linewidth',2,'MarkerFaceColor','auto')
    
    fprintf('Epi stimuli. Model prediction: %.2f - %.2f, experimental data: %.2f±%.2f\n',l(1), u(1), mu(1),s(1))
    fprintf('CL  stimuli. Model prediction: %.2f - %.2f, experimental data: %.2f±%.2f\n',l(2), u(2), mu(2),s(2))
    
    figure(5)
    subplot(1,2,1)
    set(gcf,'OuterPosition',[100,100,500,400])
    bar(30,mu(1), 'w', 'BarWidth',2)
    hold on
    errorbar(30, mu(1),s(1),'ko','linewidth',2,'MarkerFaceColor','auto', 'CapSize',18)
    title('EPI')
    set(gca, 'YLimSpec', 'Tight');
    box off
    
    subplot(1,2,2)
    bar(30,mu(2), 'w','BarWidth',2)
    hold on
    errorbar(30, mu(2),s(2),'ko','linewidth',2,'MarkerFaceColor','auto','CapSize',18)
    title('CL')
    set(gca, 'YLimSpec', 'Tight');
    box off
end

cd('..')

%% Plotting functions
function [] = Plot_Exocytosis_Uncertainty(sim, boundry, expData)
design = sim.Properties.RowNames';
design(strcmp(design,'Time'))=[];
extra=7;
for e=design
    switch e{:}
        case 'Ca_ATP',            figure(1), subplot(4,2,3)
        case 'Ca_noATP',          figure(1), subplot(4,2,4)
        case 'noCa_ATP',          figure(1), subplot(4,2,5)
        case 'noCa_noATP',        figure(1), subplot(4,2,6)
        case 'highCa_ATP',        figure(1), subplot(4,2,7)
            
        case 'Ca_ATP_nocAMP',     figure(2), subplot(3,1,1)
        case 'noCa_ATP_nocAMP',   figure(2), subplot(3,1,2)
        case 'noCa_noATP_nocAMP', figure(2), subplot(3,1,3)
            
        case 'EPI_ATP',           figure(1), subplot(4,2,1)
        case 'CL_Ca',            figure(1), subplot(4,2,2)
            
        case 'CL_ATP'
            figure(1);
            hold on
            shadeColor=[0.95,0.65,0];
            Plotcolor=[1, 0.2, 0];
            simTime=sim{'Time','Measures'}./60;
            subplot(2,1,1)
            hold on
            title('Simulation')
            xx=[simTime fliplr(simTime) ];
            yy=[boundry{e,'Min'} fliplr(boundry{e,'Max'})];
            h=fill(xx,yy,shadeColor);
            set(h, 'EdgeColor','None');
            plot(simTime,sim{e,'Measures'},'color', Plotcolor , 'LineWidth' , 1.5 ) %Plot the best parameter set
            xlabel('Time (Min)')
            ylabel('\DeltaC/\Deltat (fF/s)')
            axis([-0.2 10 -2.5 20])
            box off
            
            subplot(2,1,2)
            errorbar(expData{'Time',e}./60,expData{'Mean',e},expData{'SEM',e},'ko', 'MarkerFaceColor','auto', 'MarkerFaceColor','auto')
            title('Experimental data','Interpreter', 'none')
            xlabel('Time (Min)')
            ylabel('\DeltaC/\Deltat (fF/s)')
            axis([-0.2 10 -2.5 20])
            box off
            sgtitle('External stimulus, CL no Ca')
            
            figure(15)
            tInd=ismember(simTime,expData{'Time',e}./60);
            y=sim{e,'Measures'};
            
            errorbar(expData{'Time',e}./60,expData{'Mean',e},expData{'SEM',e},'ko', 'MarkerFaceColor','auto', 'MarkerFaceColor','auto')
            hold on
            plot(simTime(tInd),y(tInd))
            
            cost= nansum((expData{'Mean',e}-y(tInd)).^2./expData{'SEM',e}.^2);
            dgf=length(expData{'Mean',e});
            sprintf('Prediction, CL_ATP, Cost: %.1f, dgf: %i, Chi2: %.1f', cost, dgf, chi2inv(0.95, dgf))
            
            
        otherwise,                figure(extra), extra=extra+1;
    end
    Plot_Subplot(boundry, sim, e, expData)
end

figure(1), SetSameAxes(1,[0,2],[],[],0)
end

function []=Plot_Subplot(boundry, bestSim,e, data, shadeColor)
hold on
if nargin<5 || isempty(shadeColor)
    shadeColor=[0.95,0.65,0];
    alpha=1;
else
    alpha=1;
end
bestColor=[1, 0.2, 0];
simTime=bestSim{'Time','Measures'}./60;
if any(strcmp(boundry.Properties.RowNames,e))
    xx=[simTime fliplr(simTime) ];
    yy=[boundry{e,'Min'} fliplr(boundry{e,'Max'})];
    h=fill(xx,yy,shadeColor);
    set(h, 'EdgeColor','None');
    set(h,'facealpha',alpha)
end
if any(strcmp(bestSim.Properties.RowNames,e))
    plot(simTime,bestSim{e,'Measures'},'color', bestColor , 'LineWidth' , 1.5 ) %Plot the best parameter set
end
if ~isempty(data) && any(strcmp(data.Properties.VariableNames,e))
    errorbar(data{'Time',e}./60,data{'Mean',e},data{'SEM',e},'ko', 'MarkerFaceColor','auto')
end

title(e,'Interpreter', 'none')
xlabel('Time (Min)')
ylabel('\DeltaC/\Deltat (fF/s)')
axis tight
box off
end

%% Estimating uncertainty
function [ boundry, sim] = MinMax(model, params, bestparam, expData, experiments, tend, betaReduction)
if nargin<6
    tend=60*12;
end
if nargin<7
    betaReduction=[];
end
params=[params; bestparam];
for i = 1 : size(params,1)
    try
        intialconditions = SimulateSteadyState(params(i,:), model, 10000, betaReduction);
        [sim,~,~, maxcAMP]=SimulateExperiments(params(i,:), 0:1:tend, intialconditions, model, expData, experiments); % Simulates the experiments.
        
        sim_valuesStates=sim.States;
        idxTime = ismember(sim{'Time','Measures'}, 15*60);
        
        if any(idxTime)
            idxAdi = 4;
            idxCL = strcmp(sim.Properties.RowNames,'CL_ATP');
            sim_valuesStates(1:end-1,:,idxAdi)=sim_valuesStates(1:end-1,:,idxAdi)./sim_valuesStates(idxCL,idxTime,idxAdi);
        end
        if i == 1
            minValues=sim.Measures;
            maxValues=minValues;
            
            minValuesStates=sim_valuesStates;
            maxValuesStates=minValuesStates;
        end
        if maxcAMP<inf
            sim_values=sim.Measures;
            minValues(minValues>sim_values)=sim_values(minValues>sim_values);
            maxValues(maxValues<sim_values)=sim_values(maxValues<sim_values);
            
            
            minValuesStates(minValuesStates>sim_valuesStates)=sim_valuesStates(minValuesStates>sim_valuesStates);
            maxValuesStates(maxValuesStates<sim_valuesStates)=sim_valuesStates(maxValuesStates<sim_valuesStates);
        end
        
    catch err
        disp('error in sim');
    end
    fprintf('Done with %i of %i parameter sets\n',i,size(params,1))
end
boundry=table(maxValues, minValues,maxValuesStates, minValuesStates, 'VariableNames',{'Max','Min','MaxStates','MinStates'},'RowNames',sim.Properties.RowNames);
end
