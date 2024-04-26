function [x,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,x0,options)
% Copyright 2023, Gurobi Optimization, LLC
%
%INTLINPROG A mixed integer programming (MIP) example using the
%   Gurobi MATLAB interface
%
%   This example is based on the intlinprog interface defined in the
%   MATLAB Optimization Toolbox. The Optimization Toolbox
%   is a registered trademark of The Math Works, Inc.
%
%   x = INTLINPROG(f,intcon,A,b) solves the MIP problem:
%
%     minimize     f'*x
%     subject to    A*x <= b,
%                     x(j) integer, where j is in the vector
%                                   intcon of integer constraints.
%
%   For large problems, you can pass A as a sparse matrix and b as a
%   sparse vector.
%
%   x = INTLINPROG(f,intcon,A,b,Aeq,beq) solves the MIP problem:
%
%     minimize     f'*x
%     subject to    A*x <= b,
%                 Aeq*x == beq,
%                     x(j) integer, where j is in the vector
%                                   intcon of integer constraints.
%
%   For large problems, you can pass Aeq as a sparse matrix and beq as a
%   sparse vector. You can set A=[] and b=[] if no inequalities exist.
%
%   x = INTLINPROG(f,intcon,A,b,Aeq,beq,lb,ub) solves the MIP problem:
%
%     minimize     f'*x
%     subject to    A*x <= b,
%                 Aeq*x == beq,
%           lb <=     x <= ub,
%                     x(j) integer, where j is in the vector
%                                   intcon of integer constraints.
%
%   You can set lb(j) = -inf, if x(j) has no lower bound, and ub(j) = inf,
%   if x(j) has no upper bound. You can set Aeq=[] and beq=[] if no
%   equalities exist.
%
%   x = INTLINPROG(f,intcon,A,b,Aeq,beq,lb,ub,X0) solves the problem above
%   with MIP start set to X0.
%
%   You can set lb=[] or ub=[] if no bounds exist.
%
%   x = INTLINPROG(f,intcon,A,b,Aeq,beq,lb,ub,x0,OPTIONS) solves the
%   problem above given the specified OPTIONS. Only a subset of possible
%   options have any effect:
%
%     OPTIONS.Display               'off' or 'none' disables output,
%     OPTIONS.MaxTime               time limit in seconds,
%     OPTIONS.MaxFeasiblePoints     MIP feasible solution limit,
%     OPTIONS.RelativeGapTolerance  relative MIP optimality gap,
%     OPTIONS.AbsoluteGapTolerance  absolute MIP optimality gap.
%
%   x = INTLINPROG(PROBLEM) solves PROBLEM, which is a structure that must
%   have solver name 'intlinprog' in PROBLEM.solver. You can also specify
%   any of the input arguments above using fields PROBLEM.f, PROBLEM.A, ...
%
%   [x,fval] = INTLINPROG(f,intcon,A,b) returns the objective value at the
%   solution. That is, fval = f'*x.
%
%   [x,fval,exitflag] = INTLINPROG(f,intcon,A,b) returns an exitflag
%   containing the status of the optimization. The values for exitflag and
%   the corresponding status codes are:
%
%      2  stopped prematurely, integer feasible point found
%      1  converged to a solution
%      0  stopped prematurely, no integer feasible point found
%     -2  no feasible point found
%     -3  problem is unbounded
%
%   [x,fval,exitflag,OUTPUT] = INTLINPROG(f,intcon,A,b) returns information
%   about the optimization. OUTPUT is a structure with the following fields:
%
%     OUTPUT.message          Gurobi status code
%     OUTPUT.relativegap      relative MIP optimality gap
%     OUTPUT.absolutegap      absolute MIP optimality gap
%     OUTPUT.numnodes         number of branch-and-cut nodes explored
%     OUTPUT.constrviolation  maximum violation for constraints and bounds
%

% Initialize missing arguments
if nargin == 1
    if isa(f,'struct') && isfield(f,'solver') && strcmpi(f.solver,'intlinprog')
        [f,intcon,A,b,Aeq,beq,lb,ub,x0,options] = probstruct2args(f);
    else
        error('PROBLEM should be a structure with valid fields');
    end
elseif nargin < 4 || nargin > 10
    error('INTLINPROG: the number of input arguments is wrong');
elseif nargin < 10
    options = struct();
    if nargin == 9
        if isa(x0,'struct') || isa(x0,'optim.options.SolverOptions')
            options = x0; % x0 was omitted and options were passed instead
            x0 = [];
        end
    else
        x0 = [];
        if nargin < 8
            ub = [];
            if nargin < 7
                lb = [];
                if nargin < 6
                    beq = [];
                    if nargin < 5
                        Aeq = [];
                    end
                end
            end
        end
    end
end

%build Gurobi model
model.obj = f;
model.A = [sparse(A); sparse(Aeq)]; % A must be sparse
n = size(model.A, 2);
model.vtype = repmat('C', n, 1);
model.vtype(intcon) = 'I';
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b(:); beq(:)]); % rhs must be dense
if ~isempty(x0)
    model.start = x0;
end
if ~isempty(lb)
    model.lb = lb;
else
    model.lb = -inf(n,1); % default lb for MATLAB is -inf
end
if ~isempty(ub)
    model.ub = ub;
end

% Extract relevant Gurobi parameters from (subset of) options
params = struct();

if isfield(options,'Display') || isa(options,'optim.options.SolverOptions')
    if any(strcmp(options.Display,{'off','none'}))
        params.OutputFlag = 0;
    end
end

if isfield(options,'MaxTime') || isa(options,'optim.options.SolverOptions')
    params.TimeLimit = options.MaxTime;
end

if isfield(options,'MaxFeasiblePoints') ...
        || isa(options,'optim.options.SolverOptions')
    params.SolutionLimit = options.MaxFeasiblePoints;
end

if isfield(options,'RelativeGapTolerance') ...
        || isa(options,'optim.options.SolverOptions')
    params.MIPGap = options.RelativeGapTolerance;
end

if isfield(options,'AbsoluteGapTolerance') ...
        || isa(options,'optim.options.SolverOptions')
    params.MIPGapAbs = options.AbsoluteGapTolerance;
end

params.Method = 2;
params.Crossover = 0;
params.BarConvTol = 1e-3;
params.NodeMethod = 2;
params.Presolve = 0;
params.Aggregate = 0;
% %params.NormAdjust = 0;
% %params.PreSparsify = 1;
params.Threads = 128;
params.BarHomogeneous = 1;
%params.Heuristics = 0;
% Solve model with Gurobi
result = gurobi(model, params);

% Resolve model if status is INF_OR_UNBD
if strcmp(result.status,'INF_OR_UNBD')
    params.DualReductions = 0;
    warning('Infeasible or unbounded, resolve without dual reductions to determine...');
    result = gurobi(model,params);
end

% Collect results
x = [];
output.message = result.status;
output.relativegap = [];
output.absolutegap = [];
output.numnodes = result.nodecount;
output.constrviolation = [];

if isfield(result,'x')
    x = result.x;
    if nargout > 3
        slack = model.A*x-model.rhs;
        violA = slack(1:size(A,1));
        violAeq = norm(slack((size(A,1)+1):end),inf);
        viollb = model.lb(:)-x;
        violub = 0;
        if isfield(model,'ub')
            violub = x-model.ub(:);
        end
        output.constrviolation = max([0; violA; violAeq; viollb; violub]);
    end
end

fval = [];

if isfield(result,'objval')
    fval = result.objval;
    if nargout > 3 && numel(intcon) > 0
        U = fval;
        L = result.objbound;
        output.relativegap = 100*(U-L)/(abs(U)+1);
        output.absolutegap = U-L;
    end
end

if strcmp(result.status, 'OPTIMAL')
    exitflag = 1;
elseif strcmp(result.status, 'INFEASIBLE') ...
        || strcmp(result.status, 'CUTOFF')
    exitflag = -2;
elseif strcmp(result.status, 'UNBOUNDED')
    exitflag = -3;
elseif isfield(result, 'x')
    exitflag = 2;
else
    exitflag = 0;
end

% Local Functions =========================================================

function [f,intcon,A,b,Aeq,beq,lb,ub,x0,options] = probstruct2args(s)
%PROBSTRUCT2ARGS Get problem structure fields ([] is returned when missing)

f = getstructfield(s,'f');
intcon = getstructfield(s,'intcon');
A = getstructfield(s,'Aineq');
b = getstructfield(s,'bineq');
Aeq = getstructfield(s,'Aeq');
beq = getstructfield(s,'beq');
lb = getstructfield(s,'lb');
ub = getstructfield(s,'ub');
x0 = getstructfield(s,'x0');
options = getstructfield(s,'options');

function f = getstructfield(s,field)
%GETSTRUCTFIELD Get structure field ([] is returned when missing)

if isfield(s,field)
    f = getfield(s,field);
else
    f = [];
end
