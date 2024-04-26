%   gurobi()
%
%    gurobi ( model, params )
%    gurobi ( model )
%
%    This function optimizes the given model. The algorithm used for the
%    optimization depends on the model type (simplex or barrier for a
%    continuous model; branch-and-cut for a MIP model). Upon successful
%    completion it will return a struct variable containing solution
%    information.
%
%    Please consult Variables and Constraintssection in the reference manual
%    for a discussion of some of the practical issues associated with
%    solving a precisely defined mathematical model using finite-precision
%    floating-point arithmetic.
%
%    Arguments:
%
%    model: The model struct must contain a valid Gurobi model. See the
%    model argument section for more information.
%
%    params: The params struct, when provided, contains a list of modified
%    Gurobi parameters. See the params argument section for more
%    information.
%
%    Example usage:
%     result = gurobi(model, params);
%     if strcmp(result.status, 'OPTIMAL');
%       fprintf('Optimal objective: %e\n', result.objval);
%       disp(result.x)
%     else
%       fprintf('Optimization returned status: %s\n', result.status);
%     end
%
%    All details on input and output arguments, and general
%    information on the Gurobi Optimizer Matlab Interface, are given in
%    <a href="matlab:web(fullfile(fileparts(which('gurobi')), 'html', 'index.html'))">included html documentation</a>, and on-line at the <a href="matlab:web('https://www.gurobi.com/documentation/10.0/refman/matlab_api_overview.html')">Gurobi Documentation</a> page.
%
% Copyright 2023, Gurobi Optimization, LLC
%
