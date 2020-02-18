function [xsol, fval, exitflag, output, history, options] = run_fit_fminsearchbnd(f,init,prLB,prUB,optionsStruct)
% [xsol, fval, exitflag, output, history, options] = run_fit_fminsearchbnd(f,init,prLB,prUB,optionsStruct)
%
% Wrapper for fminsearchbnd
%
% Inputs
%   f: function to optimize
%   init: initial conditions
%   prLB: lower bounds on parameters
%   prUB: upper bounds on parameters
%   optionsStruct: structure of options for optimization (fields tolFun,
%                  tolX, maxIter, maxFunEvals)
%
% Outputs
%   xsol, fval, exitflag, output: see fminsearch documentation
%   history: full history of parameters searched
%   option: structure of options used for optimization
%   

history.x = [];
history.fval = [];


if isfield(optionsStruct,'tolFun')
    tolFun = optionsStruct.tolFun;
else
    tolFun = 1e-6;
end
if isfield(optionsStruct,'tolX')
    tolX = optionsStruct.tolX;
else
    tolX = 1e-6;
end
if isfield(optionsStruct,'maxIter')
    maxIter = optionsStruct.maxIter;
else
    maxIter = 500;
end
if isfield(optionsStruct,'maxFunEvals')
    maxFunEvals = optionsStruct.maxFunEvals;
else
    maxFunEvals = 500;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% call optimization
    options = optimset('OutputFcn',@outfunfminsearchbnd,...
        'Display','iter','tolfun',tolFun,'tolx',tolX,'maxFunEvals',maxFunEvals,'maxIter',maxIter);
    [xsol, fval, exitflag, output] = fminsearchbnd(f,init,prLB,prUB,options);
   

    
function stop = outfunfminsearchbnd(x,optimValues,state)
stop = false;

switch state
    case 'iter'
        % Concatenate current point and objective function
        % value with history. x must be a row vector.
        history.fval = [history.fval; optimValues.fval];
        history.x = [history.x; x];
end
end

end