%   gurobi_relax()
%
%    gurobi_relax ( model, params )
%    gurobi_relax ( model )
%
%    Create the relaxation of a MIP model. Transforms integer variables into
%    continuous variables, and removes SOS and general constraints.
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
%    Return value:
%
%    A model struct variable, as described in the model parameter section.
%
%    Example usage:
%    model = gurobi_read('stein9.mps');
%    relaxed = gurobi_relax(model);
%
%    All details on input and output arguments, and general
%    information on the Gurobi Optimizer Matlab Interface, are given in
%    <a href="matlab:web(fullfile(fileparts(which('gurobi')), 'html', 'index.html'))">included html documentation</a>, and on-line at the <a href="matlab:web('https://www.gurobi.com/documentation/10.0/refman/matlab_api_overview.html')">Gurobi Documentation</a> page.
%
% Copyright 2023, Gurobi Optimization, LLC
%
