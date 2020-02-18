function [prsSol, lsqSol, xdata, ydata] = multistage_noise_model_optimization(dataFilename,its,runParallel,maxIter,descr,init)
% [prsSol, lsqSol, xdata, ydata] = multistage_noise_model_optimization(dataFilename,its,runParallel,maxIter,descr,init)
%
% Main function that fits multistage noise model (Weber, Shea-Brown, &
% Rieke 2020) to data
%
% Inputs
%   dataFilename: filename of data to load for fitting; must contain
%                 variables 'xdata' and 'ydata', vectors of equal length
%                 that contain the inputs and responses, respectively
%   its: iterations to run (e.g., [1 3 4] will run 3 different iterations
%        using those numbers as labels)
%   runParallel: 0 (non-parallel) or 1 (parallel, requires Parallel
%                Computing Toolbox); default: 0
%   maxIter: maximum number of iterations for optimization; default: 5000
%   descr: string to record description of dataset or these iterations
%   init: optional vector of initial conditions for optimization; default:
%         random initial conditions on each iteration   
%
% Outputs
%   prsSol: solution of optimization. Entries are:
%               upstream noise
%               multiplicative noise
%               downstream noise strength
%               nonlinearity param 1: vertical stretch
%               nonlinearity param 2: horizontal stretch
%               nonlinearity param 3: horizontal offset
%               nonlinearity param 4: vertical offset
%               downstream noise probability
%            Note that upstream noise is in the same units as xdata (not
%            z-scored), and nonlinearity params are in
%            "reparameterized" units (i.e., to be plotted with
%            nl_sr_reparam function).
%   lsqSol: least squares fit of nonlinearity to data. Nonlinearity
%           parameters are in original (not reparameterized) units (i.e.,
%           to be plotted with nl_sr function).
%   xdata: input data used for fitting
%   ydata: response data used for fitting

if ~exist('runParallel','var') || isempty(runParallel) 
    runParallel = 0;
end
if ~exist('descr','var') || isempty(descr)
    descr = '';
end
if ~exist('maxIter','var') || isempty(maxIter)
    maxIter = 5000;
end

% addpath('/home/username/chebfun-master');  % add chebfun package to path if necessary

%% if running parallel, open pool
if runParallel
    myCluster = parcluster('local');
    poolObj = parpool('local',myCluster.NumWorkers);
end

%% load data files and find best nonlinearity with least squares
load([dataFilename '.mat'])  % data file must contain variables 'xdata' and 'ydata'
best = 1e10;
LsqOptions = optimset('TolX',1e-16,'TolFun',1e-16,'MaxFunEvals',10000,'MaxIter',1000,'display','off');
lsqSol = [.001 .1 0 0];
nl = @nl_sr; 
if isrow(xdata)
    xdata = xdata';
    ydata = ydata';
end
for itLSQ = 1:3
    [xCurr, resnorm] = lsqcurvefit(nl,lsqSol,xdata,ydata,[],[],LsqOptions);
    
    if resnorm<best
        lsqSol = xCurr;
        best = resnorm;
    end
end


%% set up hyperparameters for optimization
MnmOptions = [];
MnmOptions.tolFun = 1e-12;
MnmOptions.tolX = 1e-12;
MnmOptions.maxIter = maxIter;
MnmOptions.maxFunEvals = maxIter;

for it = its  % for the desired number of iterations (different initial conditions)
    
    fnameSave = ['./savedResults/mnm_optim_results_' dataFilename '_' num2str(it)];
    fid = fopen([fnameSave '.txt'],'w');
    fclose(fid);
    
    if runParallel
        f = @(prs) -calculate_mnm_likelihood_parallel([prs(1:4) prs(5)/prs(4) prs(6)/prs(4) prs(7:8)],xdata,ydata,fnameSave); % note: optimization performs better when parameters are normalized in this way
    else
        f = @(prs) -calculate_mnm_likelihood([prs(1:4) prs(5)/prs(4) prs(6)/prs(4) prs(7:8)],xdata,ydata,fnameSave);
    end
    
    % upper and lower bounds on parameters
    prLB = [ 0.01 * std(xdata)   5e-2   1e-3   1e-2   1e-4*230/std(xdata)   -200   0   0];
    prUB = [    2 * std(xdata)   10     15     100    230/std(xdata)        200    5   1];
    
    
    %%% find set of initial conditions that gives non-zero likelihood
    lsqSolTrans = [lsqSol(1) lsqSol(2)*lsqSol(1) lsqSol(3)*lsqSol(1) lsqSol(4)];  
    if ~exist('init','var') || isempty(init)
        init = [prLB(1:3) + .5*(prUB(1:3)-prLB(1:3)) + (rand-.5)*(prUB(1:3)-prLB(1:3)) lsqSolTrans + .4*lsqSolTrans.*(2*rand(size(lsqSolTrans))-1)  rand]; % perturb around lsqSol
        init(init<prLB) = prLB(init<prLB); 
        init(init>prUB) = prUB(init>prUB);
        fInit = f(init);
        while fInit == Inf || ~isreal(fInit)
            initTemp = [prLB(1:3) + .5*(prUB(1:3)-prLB(1:3)) + (rand-.5)*(prUB(1:3)-prLB(1:3)) lsqSolTrans + .4*lsqSolTrans.*(2*rand(size(lsqSolTrans))-1)  rand]; 
            initTemp(initTemp<prLB) = prLB(initTemp<prLB);
            initTemp(initTemp>prUB) = prUB(initTemp>prUB);
            init = initTemp;
            fInit = f(init);
        end
    end
    
    %%% save initial conditions
    save([fnameSave '.mat'],'descr','xdata','ydata','lsqSol','init','fInit')
    
    %%% run optimization
    disp('start fit')
    tic
    [prsSol, fval, exitflag, output, history, ~] = run_fit_fminsearchbnd(f,init,prLB,prUB,MnmOptions);
    timeToRun = toc;
    disp(['time to run: ' num2str(timeToRun) ' sec'])
    
    %%% save output
    save([fnameSave '.mat'],'xdata','ydata','lsqSol','timeToRun','prsSol', 'fval', 'exitflag', 'output','init','MnmOptions','history','prLB','prUB','descr')
    clear init
    
end

if runParallel
    delete(poolObj)
end