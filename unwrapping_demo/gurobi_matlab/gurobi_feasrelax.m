%   gurobi_feasrelax()
%
%    gurobi_feasrelax ( model, relaxobjtype, minrelax, penalties, params )
%    gurobi_feasrelax ( model, relaxobjtype, minrelax, penalties )
%
%    This function computes a feasibility relaxation for the input model
%    argument. The feasibility relaxation is a model that, when solved,
%    minimizes the amount by which the solution violates the bounds and
%    linear constraints of the original model. You must provide a penalty to
%    associate with relaxing each individual bound or constraint (through
%    the penalties argument). These penalties are interpreted in different
%    ways, depending on the value of the relaxobjtype argument.
%
%    Arguments:
%
%    model: The model struct must contain a valid Gurobi model. See the
%    model argument section for more information.
%
%    relaxobjtype: The approach used to impose penalties on violations.
%    If you specify relaxobjtype=0, the objective for the feasibility
%    relaxation is to minimize the sum of the weighted magnitudes of the
%    bound and constraint violations.
%    If you specify relaxobjtype=1, the objective for the feasibility
%    relaxation is to minimize the weighted sum of the squares of the bound
%    and constraint violations.
%    If you specify relaxobjtype=2, the objective for the feasibility
%    relaxation is to minimize the weighted count of bound and constraint
%    violations.
%    In all cases, the weights are taken from penalties.lb, penalties.ub and
%    penalties.rhs. You can provide the special penalty value Inf to
%    indicate that the corresponding bound or constraint cannot be relaxed.
%
%    minrelax: The minrelax argument is a boolean that controls the type of
%    feasibility relaxation that is created. If minrelax=false, optimizing
%    the returned model gives a solution that minimizes the cost of the
%    violation. If minrelax=true, optimizing the returned model finds a
%    solution that minimizes the original objective, but only from among
%    those solutions that minimize the cost of the violation. Note that
%    gurobi_feasrelax must solve an optimization problem to find the minimum
%    possible relaxation when minrelax=true, which can be quite expensive.
%
%    penalties: The penalties argument is a struct array, having the
%    following optional fields (default: all Inf):
%    lb Penalty for violating each lower bound.
%    ub Penalty for violating each upper bound.
%    rhs Penalty for violating each constraint.
%
%    To give an example, if a constraint with penalties.rhs value p is
%    violated by 2.0, it would contribute 2*p to the feasibility relaxation
%    objective for relaxobjtype=0, 2*2*p for relaxobjtype=1, and p for
%    relaxobjtype=2.
%
%    params: The params struct, when provided, contains a list of modified
%    Gurobi parameters. See the params argument section for more
%    information.
%
%    Return value:
%
%    A struct containing two fields:
%    result.model, a struct variable, as described in the model argument
%    section.
%    result.feasobj, a scalar. If minrelax==true this is the relaxation
%    problem objective value, 0.0 otherwise.
%
%    Example usage:
%    model = gurobi_read('stein9.mps');
%    penalties.lb = ones(length(model.lb),1);
%    penalties.ub = ones(length(model.ub),1);
%    penalties.rhs = ones(length(model.rhs),1);
%    feasrelaxresult = gurobi_feasrelax(model, 0, false, penalties);
%
%    All details on input and output arguments, and general
%    information on the Gurobi Optimizer Matlab Interface, are given in
%    <a href="matlab:web(fullfile(fileparts(which('gurobi')), 'html', 'index.html'))">included html documentation</a>, and on-line at the <a href="matlab:web('https://www.gurobi.com/documentation/10.0/refman/matlab_api_overview.html')">Gurobi Documentation</a> page.
%
% Copyright 2023, Gurobi Optimization, LLC
%
