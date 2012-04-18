function [T, Y, E, P] = main_remodeling_sim(fobs,fnet,MaxIter)
% 
% Dynamics of the ethanolamine glycerophospholipid lipid remodeling network
%
% Copyright Mar 2012, Lu Zhang, all rights reserved.
% Boston College
%
% Usage: 
% main_remodeling_sim('../demo/DATA/18_0-18_1.txt','../demo/DATA/18_0-18_1.sif',1000);
%
% input
%
% fobs	lipid pulse-chase mass spec data
% fnet	remodeling network
% MaxIter	maximum number of iterations to update parameter values, default 1000
%
% output
%
% T	collocation time points
% Y	predicted values at collocation points
% E	recorded error for each iteration
% P	recorded parameter values for each iteration
%

%%
if nargin < 3
	MaxIter = 1000;
end
is_plot = 'T';

% Use medium-scale algorithm for lsqlin 
options = optimset('LargeScale','off','Display','off');

% Generate Output Directory
[pathstr, prefix] = fileparts(fobs);
if ~isempty(pathstr); 
    outdir = [regexprep(pathstr, 'DATA', 'results') filesep prefix];
else
    outdir = [prefix];
end
if ~exist(outdir)
    mkdir(outdir);
end

% Input flux network
Network = importdata(fnet);
Network = unique(Network,'rows');

% Input pulse-chase measurement
temp = importdata (fobs,'\t',1);
Species = temp.textdata(2:size(temp.textdata,1),1);
Time = str2num(char(temp.textdata(1,3:2:size(temp.textdata,2))))';
Data = temp.data(:,1:2:size(temp.data,2))';
STD = temp.data(:,2:2:size(temp.data,2))';


% B-spline Algorithm Parameter Settings
num_Spline_i = 5;
num_Spline_f = 21;
beta = 10; % lamda = beta^2;
dt = 1/2;
thu = 0.0001;
du = 0.01;

num_para = length(unique(Network(:,2)));
num_sp = length(Species);
num_time = length(Time);
t_start = Time(1);
t_end = Time(num_time);
num_col = (t_end - t_start)/dt +1;
T_col = t_start : dt : t_end;

%% Start B-spline algorithm
A1 = construct_BSpline_matrix( Time , num_Spline_i-1 , t_start , t_end , 'B' );
Matrix_A1 = zeros(num_time*num_sp,num_Spline_i*num_sp);
for i = 1 : num_sp
    Matrix_A1((i-1)*num_time+1:i*num_time,(i-1)*num_Spline_i+1:i*num_Spline_i) = A1;
end

A2 = construct_BSpline_matrix( T_col , num_Spline_i-1 , t_start , t_end , 'B');
Matrix_A2 = zeros(num_col*num_sp,num_Spline_i*num_sp);
for i = 1 : num_sp
    Matrix_A2((i-1)*num_col+1:i*num_col,(i-1)*num_Spline_i+1:i*num_Spline_i) = A2;
end
Coeff = lsqlin(Matrix_A1,Data(:),-Matrix_A2,zeros(num_col*num_sp,1),[],[],[],[],[],options);

Fit_col = reshape(Matrix_A2*Coeff,num_col,num_sp);

A_Time_B = construct_BSpline_matrix(Time , num_Spline_f-1 , t_start , t_end ,'B');
A_col_dB = construct_BSpline_matrix(T_col, num_Spline_f-1 , t_start , t_end ,'dB');
A_col_B = construct_BSpline_matrix(T_col, num_Spline_f-1 , t_start , t_end ,'B');

Matrix_Time = zeros(num_time*num_sp,num_Spline_f*num_sp);
for i = 1 : num_sp
    Matrix_Time((i-1)*num_time+1:i*num_time,(i-1)*num_Spline_f+1:i*num_Spline_f) = A_Time_B.*beta;
end

Matrix_col = zeros(num_col*num_sp,num_Spline_f*num_sp);
for i = 1 : num_sp
    Matrix_col((i-1)*num_col+1:i*num_col,(i-1)*num_Spline_f+1:i*num_Spline_f) = A_col_dB;
end

Z = [zeros(num_col*num_sp,1); Data(:).*beta];
A = [zeros(num_col*num_sp,num_para) Matrix_col ; zeros(num_time*num_sp,num_para) Matrix_Time];

% Iterate to optimize parameter values
num_iter = 1;
last_err = 0;
X0 = [];
E = [];
P = [];
while (1)
    A(1:num_col*num_sp,1:num_para) = zeros(num_col*num_sp,num_para);
    for alpha = 1 : num_sp
        for t = 1 : num_col
            for k = 1: num_para
                inx = find(Network(:,1)==alpha & Network(:,2)==k);
                if (~isempty(inx))
                   A((alpha-1)*num_col+t,k) = A((alpha-1)*num_col+t,k) + Fit_col(t,alpha)*length(inx);
                end
                inx = find(Network(:,3)==alpha & Network(:,2)==k);
                if (~isempty(inx))
                    A((alpha-1)*num_col+t,k) = A((alpha-1)*num_col+t,k) - sum(Fit_col(t,Network(inx,1)));
                end
            end
        end
    end
    [X,resnorm] = lsqlin(A,Z,[],[],[],[],[zeros(num_para,1) ],[],X0);
    Coeff = reshape (X(num_para+1:length(X)),num_Spline_f,num_sp);
    Fit_col = A_col_B * Coeff;
    [num_iter resnorm resnorm-last_err]
    E(num_iter) = resnorm;
    P(num_iter,:) = X(1:num_para);
    if (abs(resnorm-last_err)<thu)
        break;
    elseif (num_iter>MaxIter & resnorm<min(E(num_iter-3:num_iter-1)))  %m-3 or m-5
        break;
    elseif (num_iter>MaxIter*2)
        break;
    end
    num_iter = num_iter + 1;
    X0 = X;
    last_err = resnorm;
end

%% Output results

% LOG file
X2 = X(1:num_para);
[B inx] = sort(X2,'descend');

fid = fopen([outdir filesep 'log_' prefix '.txt'],'wt');
fprintf(fid,'%s\n\nLipid\t%s\nNetwork\t%s\nnum_iter\t%d\nerror\t%1.2f\n\n',date,fobs,fnet,num_iter,resnorm);
for j = 1 : num_para
    i = inx(j);
    a = find(Network(:,2)==i);
    if strcmp(Species{Network(a(1),1)}(1:4),Species{Network(a(1),3)}(1:4))
        fprintf(fid,'%d sn2 %s -> %s %1.4f\n',i,Species{Network(a(1),1)}(6:9),Species{Network(a(1),3)}(6:9),X2(i));
        fprintf('%d sn2 %s -> %s %1.4f\n',i,Species{Network(a(1),1)}(6:9),Species{Network(a(1),3)}(6:9),X2(i))
    else
        fprintf(fid,'%d sn1 %s -> %s %1.4f\n',i,Species{Network(a(1),1)}(1:4),Species{Network(a(1),3)}(1:4),X2(i));
        fprintf('%d sn1 %s -> %s %1.4f\n',i,Species{Network(a(1),1)}(1:4),Species{Network(a(1),3)}(1:4),X2(i))
    end
end
if (num_iter>MaxIter*2)
    fprintf(fid,'\nWarning: algorithm may not converge!\n');
end
fclose(fid);

tempfilename = ['..' filesep 'tmp' filesep 'tp16c19154_9bed_4a92_9441_4b248cb094ee'];
save tempfilename X2 num_para num_sp Network;
[T,Y] = ode45(@mysim,[t_start:du:t_end],Fit_col(1,:)); 

% PLOT simulation results
COLOR = repmat({[0 0 1],[0 0.4980 0],[0.8471 0.1608 0],[0 0.7490 0.7490],[0.7490 0 0.7490],[0.7490 0.7490 0],[0.3137 0.3176 0.3137]},1,4);
marker = [repmat({'o'},7,1); repmat({'s'},7,1); repmat({'d'},7,1);repmat({'+'},7,1)];

if (is_plot=='T')
    figure();
    for i = 1:num_sp
        errorbar(repmat(Time',1,1),Data(:,i),STD(:,i),'LineStyle','none','Marker',marker{i},'Color',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerFaceColor',COLOR{i}); hold on
    end
    xlim([t_start t_end]);
    xlabel('Chase time (h)','FontName','Arial','FontSize',12);
    ylabel('Relative abundance (%)','FontName','Arial','FontSize',12);
    plot(T,Y);hold on
    legend_handle=legend(Species,'Location','NorthEastOutside','FontName','Arial','FontSize',11);
    yl = get(gca,'Ylim');
    set(gca,'Ylim',[-4 yl(2)]);
    set(legend_handle, 'Box', 'off')
    print('-f1','-djpeg','-r300',[outdir filesep 'sim_' prefix '.jpg']);
end

fprintf('\nOutput to Directory: %s\n', outdir)

