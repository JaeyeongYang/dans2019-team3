cd('C:\Users\Seungwoo\Documents\GitHub\dans2019-ccslab')
Data = xlsread('behav.xlsx');
% Data : 3943 * 14 variables
gain = Data(:,6);
loss = Data(:,7);
gamble = Data(:,14);
RT = Data(:,11);
Matrix_Total = zeros(16,16);
Matrix_gamble = zeros(16,16);
Matrix_RT = zeros(16,16);
for i=1:length(gain)
    index2 = ((gain(i)-9)+1)/2;
    index1 = loss(i)-4;
    
    Matrix_gamble(index1,index2) = Matrix_gamble(index1,index2) + gamble(i);
    Matrix_Total(index1,index2) = Matrix_Total(index1,index2) + 1;
    Matrix_RT(index1,index2) = Matrix_RT(index1,index2) + RT(i);
end

%% Figure1,B
FigureB = Matrix_gamble ./ Matrix_Total;
fig = figure; fig.Position = [0 0 1000 1000];
imagesc(FigureB);
colormap('jet');
g=gca;
g.XDir = 'normal';
g.YDir = 'rev';
g.CLim = [0 1]; % color rescaling
c=colorbar;
xlabel('Potential Gain ($)');
ylabel('Potential Loss ($)');
title('Probability of acceptance');
% Save
SaveRoot=['C:\Users\Seungwoo\Documents\SNU-DANS'];
FileName=['Probabiliyt of acceptance'];
if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
cd(SaveRoot);
fig.PaperPositionMode = 'auto';
print([FileName '.png'],'-dpng','-r0')
print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');

%% Figure1,C
FigureC = Matrix_RT ./ Matrix_Total;
fig = figure; fig.Position = [0 0 1000 1000];
imagesc(FigureC);
colormap('jet');
g=gca;
g.XDir = 'normal';
g.YDir = 'rev';
g.CLim = [1.1 1.8]; % color rescaling
c=colorbar;
xlabel('Potential Gain ($)');
ylabel('Potential Loss ($)');
title('Response Time');
%save
SaveRoot=['C:\Users\Seungwoo\Documents\SNU-DANS'];
FileName=['Response Time'];
if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
cd(SaveRoot);
fig.PaperPositionMode = 'auto';
print([FileName '.png'],'-dpng','-r0')
print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');


%%
clear all
cd('C:\Users\Seungwoo\Dropbox\DANS_project')
list=dir('parameter_loss_aversion.tsv');
Data = tdfread([list(1).name]);

clear x y
for i=1:16
    lamda=Data.lambda(i);
    rho = Data.rho(i);
    x = -20:0.1:40;
    for k=1:length(x)
        if x(k) >= 0
            y(i,k) = x(k)^rho;
        elseif x(k) < 0
            y(i,k) = -lamda * (-x(k))^rho;
        end
    end
end

fig=figure;
shadedErrorBar(x,mean(y), std(y))
save('Individual_Rho_lamda.mat','x','y');
fig=figure;
plot(x,y,'k-');

%% hBaseDM_Individual data
clear all
cd('E:\Document')
list=dir('parameter_loss_aversion.tsv');
Data = tdfread([list(1).name]);

clear x y
for i=1:17
    lamda=Data.lambda(i);
    rho = Data.rho(i);
    x = -40:0.1:40;
    for k=1:length(x)
        if x(k) >= 0
            y(i,k) = x(k)^rho;
        elseif x(k) < 0
            y(i,k) = -lamda * (-x(k))^rho;
        end
    end
end

fig=figure; hold on
fig.Position = [0 0 800 800];
for i=1:16
    plot(x,y(i,:),'-');
end
p1 = plot(x,y(17,:),'r'); p1.LineWidth=2;
plot(x,x,'k-');

SaveRoot=['E:\Document\DANS-ANN'];
FileName=['Individual loss aversion'];
SaveFigure(gcf,SaveRoot, FileName);

%% Low rho, Individual, Fig1 B,C
for t=1:1
    clear all;
    list=dir('sub*');
    
    k=1;
    clear Data1 Data2 Data3
    for n=1:1
        i=5;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=15;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=14;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=9;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=4;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    %
    gain = data.gain;
    loss = data.loss;
    gamble = data.gamble;
    RT = data.response_time;
    
    Matrix_Total = zeros(16,16);
    Matrix_gamble = zeros(16,16);
    Matrix_RT = zeros(16,16);
    for i=1:length(gain)
        index2 = ((gain(i)-9)+1)/2
        index1 = loss(i)-4
        
        Matrix_gamble(index1,index2) = Matrix_gamble(index1,index2) + gamble(i);
        Matrix_Total(index1,index2) = Matrix_Total(index1,index2) + 1;
        Matrix_RT(index1,index2) = Matrix_RT(index1,index2) + RT(i);
    end
    
    % Figure1,B
    FigureB = Matrix_gamble ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureB);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0 1]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Probability of acceptance');
    % Save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Probabiliyt of acceptance_low rho'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
    
    % Figure1,C
    FigureC = Matrix_RT ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureC);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [1 2]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Response Time');
    %save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Response Time_low rho'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
end

%% high rho, Individual, Fig1 B,C
for t=1:1
    clear all;
    list=dir('sub*');
    
    k=1;
    clear Data1 Data2 Data3
    for n=1:1
        i=11;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=1;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=8;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=16;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=10;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    %
    gain = data.gain;
    loss = data.loss;
    gamble = data.gamble;
    RT = data.response_time;
    
    Matrix_Total = zeros(16,16);
    Matrix_gamble = zeros(16,16);
    Matrix_RT = zeros(16,16);
    for i=1:length(gain)
        index2 = ((gain(i)-9)+1)/2
        index1 = loss(i)-4
        
        Matrix_gamble(index1,index2) = Matrix_gamble(index1,index2) + gamble(i);
        Matrix_Total(index1,index2) = Matrix_Total(index1,index2) + 1;
        Matrix_RT(index1,index2) = Matrix_RT(index1,index2) + RT(i);
    end
    
    % Figure1,B
    FigureB = Matrix_gamble ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureB);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0 1]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Probability of acceptance');
    % Save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Probabiliyt of acceptance_high rho'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
    
    % Figure1,C
    FigureC = Matrix_RT ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureC);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0.85 1.8]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Response Time');
    %save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Response Time_high rho'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
end

%% middle rho, Individual, Fig1 B,C
for t=1:1
    clear all;
    list=dir('sub*');
    
    k=1;
    clear Data1 Data2 Data3
    for n=1:1
        i=2;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=3;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=6;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=7;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=12;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end

        clear Data1 Data2 Data3
    for n=1:1
        i=13;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    %
    gain = data.gain;
    loss = data.loss;
    gamble = data.gamble;
    RT = data.response_time;
    
    Matrix_Total = zeros(16,16);
    Matrix_gamble = zeros(16,16);
    Matrix_RT = zeros(16,16);
    for i=1:length(gain)
        index2 = ((gain(i)-9)+1)/2
        index1 = loss(i)-4
        
        Matrix_gamble(index1,index2) = Matrix_gamble(index1,index2) + gamble(i);
        Matrix_Total(index1,index2) = Matrix_Total(index1,index2) + 1;
        Matrix_RT(index1,index2) = Matrix_RT(index1,index2) + RT(i);
    end
    
    % Figure1,B
    FigureB = Matrix_gamble ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureB);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0 1]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Probability of acceptance');
    % Save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Probabiliyt of acceptance_middle rho'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
    
    % Figure1,C
    FigureC = Matrix_RT ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureC);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0.85 1.8]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Response Time');
    %save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Response Time_middle rho'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
end


%% Low lambda, Individual, Fig1 B,C
for t=1:1
    clear all;
    list=dir('sub*');
    
    k=1;
    clear Data1 Data2 Data3
    for n=1:1
        i=1;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=11;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=7;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=10;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=8;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    %
    gain = data.gain;
    loss = data.loss;
    gamble = data.gamble;
    RT = data.response_time;
    
    Matrix_Total = zeros(16,16);
    Matrix_gamble = zeros(16,16);
    Matrix_RT = zeros(16,16);
    for i=1:length(gain)
        index2 = ((gain(i)-9)+1)/2
        index1 = loss(i)-4
        
        Matrix_gamble(index1,index2) = Matrix_gamble(index1,index2) + gamble(i);
        Matrix_Total(index1,index2) = Matrix_Total(index1,index2) + 1;
        Matrix_RT(index1,index2) = Matrix_RT(index1,index2) + RT(i);
    end
    
    % Figure1,B
    FigureB = Matrix_gamble ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureB);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0 1]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Probability of acceptance');
    % Save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Probabiliyt of acceptance_low lambda'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
    
    % Figure1,C
    FigureC = Matrix_RT ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureC);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [1 2]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Response Time');
    %save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Response Time_low lambda'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
end

%% high lambda, Individual, Fig1 B,C
for t=1:1
    clear all;
    list=dir('sub*');
    
    k=1;
    clear Data1 Data2 Data3
    for n=1:1
        i=2;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=3;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=9;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=13;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=15;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    %
    gain = data.gain;
    loss = data.loss;
    gamble = data.gamble;
    RT = data.response_time;
    
    Matrix_Total = zeros(16,16);
    Matrix_gamble = zeros(16,16);
    Matrix_RT = zeros(16,16);
    for i=1:length(gain)
        index2 = ((gain(i)-9)+1)/2
        index1 = loss(i)-4
        
        Matrix_gamble(index1,index2) = Matrix_gamble(index1,index2) + gamble(i);
        Matrix_Total(index1,index2) = Matrix_Total(index1,index2) + 1;
        Matrix_RT(index1,index2) = Matrix_RT(index1,index2) + RT(i);
    end
    
    % Figure1,B
    FigureB = Matrix_gamble ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureB);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0 1]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Probability of acceptance');
    % Save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Probabiliyt of acceptance_high lambda'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
    
    % Figure1,C
    FigureC = Matrix_RT ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureC);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0.9 2]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Response Time');
    %save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Response Time_high lambda'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
end

%% middle lambda, Individual, Fig1 B,C
for t=1:1
    clear all;
    list=dir('sub*');
    
    k=1;
    clear Data1 Data2 Data3
    for n=1:1
        i=4;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=5;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=6;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=12;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    clear Data1 Data2 Data3
    for n=1:1
        i=14;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end

        clear Data1 Data2 Data3
    for n=1:1
        i=16;
        File1=list(3*i-2).name;
        File2=list(3*i-1).name;
        File3=list(3*i).name;
        
        Data1 = tdfread([File1]);
        Data2 = tdfread([File2]);
        Data3 = tdfread([File3]);
        
        for m=1:length(Data1.respcat)
            if (Data1.respcat(m) ~= -1)
                data.gain(k) = Data1.gain(m);
                data.loss(k) = Data1.loss(m);
                data.gamble(k) = Data1.respcat(m);
                data.response_time(k) = Data1.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data2.respcat)
            if (Data2.respcat(m) ~= -1)
                data.gain(k) = Data2.gain(m);
                data.loss(k) = Data2.loss(m);
                data.gamble(k) = Data2.respcat(m);
                data.response_time(k) = Data2.response_time(m);
                k=k+1;
            end
        end
        for m=1:length(Data3.respcat)
            if (Data3.respcat(m) ~= -1)
                data.gain(k) = Data3.gain(m);
                data.loss(k) = Data3.loss(m);
                data.gamble(k) = Data3.respcat(m);
                data.response_time(k) = Data3.response_time(m);
                k=k+1;
            end
        end
    end
    
    %
    gain = data.gain;
    loss = data.loss;
    gamble = data.gamble;
    RT = data.response_time;
    
    Matrix_Total = zeros(16,16);
    Matrix_gamble = zeros(16,16);
    Matrix_RT = zeros(16,16);
    for i=1:length(gain)
        index2 = ((gain(i)-9)+1)/2
        index1 = loss(i)-4
        
        Matrix_gamble(index1,index2) = Matrix_gamble(index1,index2) + gamble(i);
        Matrix_Total(index1,index2) = Matrix_Total(index1,index2) + 1;
        Matrix_RT(index1,index2) = Matrix_RT(index1,index2) + RT(i);
    end
    
    % Figure1,B
    FigureB = Matrix_gamble ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureB);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0 1]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Probability of acceptance');
    % Save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Probabiliyt of acceptance_middle lambda'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
    
    % Figure1,C
    FigureC = Matrix_RT ./ Matrix_Total;
    fig = figure; fig.Position = [0 0 1000 1000];
    imagesc(FigureC);
    colormap('jet');
    g=gca;
    g.XDir = 'normal';
    g.YDir = 'rev';
    g.CLim = [0.9 2]; % color rescaling
    c=colorbar;
    xlabel('Potential Gain ($)');
    ylabel('Potential Loss ($)');
    title('Response Time');
    %save
    SaveRoot=['E:\Document\DANS-ANN'];
    FileName=['Response Time_middle lambda'];
    if ~exist([SaveRoot],'dir') mkdir([SaveRoot]); end
    cd(SaveRoot);
    fig.PaperPositionMode = 'auto';
    print([FileName '.png'],'-dpng','-r0')
    print( gcf, '-painters', '-r300', [FileName '.pdf'], '-dpdf');
end
