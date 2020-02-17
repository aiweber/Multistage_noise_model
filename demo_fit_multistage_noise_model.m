% This code fits a multiply stochastic model to neural responses. (See
% Weber, Shea-Brown, & Rieke: "Identification of multiple noise sources
% improves estimation of neural responses across stimulus conditions." bioRxiv 2020)

%% load parameters
dataFilename = 'test_data_downMixture'; % test_data_downMixture or test_data_allGaussian
iterationIDs = 1:2;                     % will run optimization for two different initial conditions, saved as "_1" and "_2"
runParallel = 1;                        % 0 will run without opening up a parallel pool of workers; 1 runs in parallel
maxIter = 10000;                        % maximum number of iterations to run optimization
descr = 'example dataset with downstream noise drawn from mixture model';  % optional description for this set of runs
init = [];                              % option to specify initial conditions of optimization

%% run optimization
[prsSol, lsqSol, xdata, ydata] = multistage_noise_model_optimization(dataFilename,iterationIDs,runParallel,maxIter,descr,init);

%% plot results
figure; hold on;
plot(xdata,ydata,'.','handlevisibility','off');
xVec = min(xdata):max(xdata);
plot(xVec,nl_sr_reparam(prsSol(4:7),xVec),'linewidth',2)
plot(xVec,nl_sr(lsqSol,xVec));
legend({'multistage noise model'; 'least squares'},'location','nw')
xlabel('inputs')
ylabel('responses')