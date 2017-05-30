%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     %
% Stochastic Tn-TM bindign simulation %
%                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: simulation_stoch_2.m, First Created on October 2013
%
% DESCRIPTION: Simulation of a two-compoent reaction with a compoent being
%              a lattice
%
% Tested with Matlab R2011b
% Christopher Solis, South Dakota State University, updated July 2016
%
% BEFORE WE BEGIN

% BEGIN
clear all;
clc;
PathName = pwd;
%Constants
TIME = 40;
VOLUME = 10000;
NN = 100;

k_on = 220;
k_off = 2.6;
sigma = 3.7; % when there is cooperativity sigma>1;
Totbins = round(logspace(1,3.5,NN));
ReservoirN = round(logspace(0,3.7,NN)); 
%
% Generate random filament lengths following an exponential distribution
Pn = rand(Totbins(end),1);
n_mu = 40;
ni = round(-n_mu.*log(Pn));

%% Matrix
Ym = [];
Free_Bm = [];
InteractConts = [];

for N = 1:NN
AB_save = [];
Free_B = [];
NUMBINS = Totbins(N);
ReservoirB = ReservoirN(N); % number of free particles of A. 
ReservoirB_t = ReservoirB;
A = zeros(1,NUMBINS); % filament length
AB_save(1,:) = A;

% Select the location of the divisions in the filament template NUMBINS
DIVS1 = 1;
DIVS2 = 1;
Count = 1;
if ni(1) >NUMBINS   
DIVS1 = 1;
else
while DIVS2(end)<NUMBINS  
DIVS1 = DIVS2;
DIVS2 = [DIVS2;sum(ni(1:Count))];   
Count = Count+1;
end
end

for i=1:TIME;
    
    r=rand(1,NUMBINS);
    for j=1:NUMBINS
        if A(1,j) == 0 % JUMP ON
            if r(1,j)<k_on*(ReservoirB_t/VOLUME)*((NUMBINS-sum(A))/VOLUME);
                A(1,j) = 1;
            end
        else % % JUMP OFF
            num_contacts = 0;
            if j > 1% not left edge
                if ismember(j,DIVS1)
                else
                num_contacts = num_contacts + A(j-1);
                end
            end
            if j< NUMBINS % not right edge
                if ismember(j,DIVS1-1)
                else
                num_contacts = num_contacts + A(j+1);
                end
            end            
            k_off_eff = k_off*((1/sigma)^num_contacts)*(sum(A)./VOLUME);
            if r(1,j)<k_off_eff 
                A(1,j) = 0;
            end        
        end
    end
    % recycled variables
    ReservoirB_t = ReservoirB - sum(A); % Current cocnetration of A
    
    AB_save(i,:) = A;
    Free_B(i) = ReservoirB_t;
end
%saved variables
    VALUES = round(0.05*TIME);
    Ym(N) = mean(sum(AB_save(TIME-VALUES:TIME,:),2))/NUMBINS;
    Free_Bm(N) = mean(Free_B(TIME-VALUES:TIME))/VOLUME;
    CC = [];
    CC = bwconncomp(AB_save(TIME,:));
    LL = []; % Length determination!
    for k = 1:length(CC.PixelIdxList);
        LL(k) = length(CC.PixelIdxList{k});
        
    end
    Length(N) = mean(LL(LL>1));
    InteractConts(N,:) = [sum(LL(LL==1))./NUMBINS sum(LL(LL==2))./NUMBINS sum(LL(LL>2))./NUMBINS];
    Percent = N*100/NN
end
%%
UNIT = 560000000;
X_DAT = ReservoirN/UNIT;
Y_DAT = Ym;
%plot(Free_Bm,Ym,'.k')
%semilogx(X_DAT,Y_DAT,'.k');

HILL = 1;
if HILL
X_DAT2 = logspace(log10(min(X_DAT)),log10(max(X_DAT)),100);    
Options = statset('Robust','on','MaxIter',100);
modelFunBind =  @(p,x) p(1) + p(2)*((1+2*p(3).*x - sqrt((1+2*p(3).*x).^2 - 4*p(3).^2.*x.^2))./(2.*p(3).*x));
startingVals = [0.2, 1.9, 1,1E-5];
[coefEsts,Resid,Jacobian] = nlinfit(X_DAT(5:end),Y_DAT(5:end), modelFunBind, startingVals,Options);
CurvData = modelFunBind(coefEsts,X_DAT2);
CoeffErr = nlparci(coefEsts,Resid,'jacobian',Jacobian,'alpha',0.05);
Kd = 1./coefEsts(3); % as M 
Kd_Err = 1./(CoeffErr(3,2)-CoeffErr(3));
Kd2 = Kd*1E6; % as uM
Kd2Err = Kd_Err*1E6;

fig1 = figure(1);
plot(X_DAT,Y_DAT,'.k','Linewidth',2,'MarkerSize',20);% DATA
hold on
semilogx(X_DAT2,CurvData,'-k','Linewidth',2);% CURVEFIT
%plot(((0.5*coefEsts(2)*coefEsts(4)-1.5*coefEsts(1)*coefEsts(4))/(0.5*coefEsts(2)+1.5*coefEsts(1)))^(1/coefEsts(3)),coefEsts(2)/2,'*g')
%legend('location', 'NorthWest');
%ylim([0 1]);
%plot(X_DAT,c2(X_DAT))
axis([0 1E-6 -0.01 1]);
hold off
set(gca,'Linewidth',2,'FontSize',18);
xlabel('Total RU concetration (M)');
ylabel('\theta_{Tn-Tpm}');
text(100.0*min(X_DAT(1)),0.15*max(CurvData),[sprintf('K_d = %0.3f',Kd2),'\pm',num2str(Kd2Err,'%0.3f'),' \muM'], 'FontName', 'Helvetica','FontSize',18,'BackgroundColor',[1 1 1]);
end
axis square
%% Plot length
fig2 = figure(2);
DO = 1;
if DO
plot(X_DAT,Length./max(Length),'.k','MarkerSize',20);
%set(gca,'Xscale','log');
xlabel('Total RU concetration (M)');
ylabel('Normalized lenght'); % normalized length of continour ligand-bound regions
xlim([1E-8,1E-5]);
%ylim([0 1])
%ylim([-0.01 1]);
axis square
end
set(gca,'Linewidth',2,'FontSize',18);
%% Types of bindign interactions
% Isolated Tn-Tpm bindign events versus contigous Tn-Tpm binding events
fig3 =figure(3);
plot(X_DAT,InteractConts(:,3),'.','Color',[0 0 0],'MarkerSize',20)% loooks better as blue though
hold on
plot(X_DAT,InteractConts(:,1),'o','Color',[0.5 0.5 0.5],'LineWidth',1.5,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1])% loooks better as blue though

plot(X_DAT,InteractConts(:,2),'.','Color',[0.5 0.5 0.5],'MarkerSize',20);% loooks better as continous red though

%set(gca,'XScale','log');
axis([0 1E-6 0 1]);
hold off
xlabel('Total Ru concetration (M)');
ylabel('Probability of interactions'); % normalized length of continour ligand-bound regions
LEG = legend('Contigous, two sides','Contigous, one side','Isolated');
set(LEG,'Box','off','Location','NorthWest');
set(LEG,'Location','NorthWest');
set(LEG,'FontSize',18);
axis square

% Areas under the curve
A1 = trapz(X_DAT,InteractConts(:,1));
A2 = trapz(X_DAT,InteractConts(:,2));
A3 = trapz(X_DAT,InteractConts(:,3));
A1p = 100*A1./sum([A1 A2 A3]); % Percent of 1-side interactions
A2p = 100*A2./sum([A1 A2 A3]); % Percent of 0-side interactions
A3p = 100*A3./sum([A1 A2 A3]); % Percent of 2-side interactions
A_weight = 1*A1p/100 + 0*A2p/100 + 2*A3p/100; % Average number of interactions
set(gca,'Linewidth',2,'FontSize',18);
%% Binding Constnat

%Teoretical
%Kb_0 = k_on./k_off

%Esperimental 
%A_eq = (NUMBINS - sum(mean(A_save(TIME-VALUES:TIME,:))))/VOLUME;
%B_eq = (mean(Free_B(TIME-VALUES:TIME)))/VOLUME; 
%AB_eq = (sum(mean(A_save(TIME-VALUES:TIME,:))))/VOLUME
%Kb = AB_eq/(A_eq*B_eq)
%Kd = 1/Kb


%% Save data into data file
%
Choice = questdlg('Save data?', 'Save3','Yes','No','No');
switch Choice
    %
    case 'Yes'  
        Msg = msgbox('Select destination','Destination Folder','help');
        uiwait(Msg);
        folder_name = uigetdir(PathName,'Select Folder'); % select folder
         print(fig1,'-depsc2', '-painters', [folder_name,'/','_BindCurv','.eps']); 
         print(fig2,'-depsc2', '-painters', [folder_name,'/','_Length','.eps']); 
         print(fig3,'-depsc2', '-painters', [folder_name,'/','_Ratios','.eps']); 
         disp('Data saved');
    case 'No'
        disp('Data was not saved');
end
