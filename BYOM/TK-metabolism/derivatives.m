%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. The
% origninal script is linked to the script files for Calanus2016, fit body residues with
% two-compartment TK. The current script is linked to Raths et al. 2023 (elimination resistant) As input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Original Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>

% * Modified script authors (receptor binding): Annika Mangold-Döring & Johannes Raths
% * Date: January 2023

% * Current version authors (biotransform): Johannes Raths
% * Date: January 2026

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

C_p   = X(1); % state 1 (internal parent concentrations) at previous time point
C_m  = X(2); % state 2 (1st BTP concentration) at previous time point
C_tot = X(3); % state 3 (total tissue concentration - parent + BTP) %JR20160113 redundant but tracked

% These concentrations are all expressed on volume basis. Make sure to
% work with molar concentrations to account for mass differences between
% parent and BTP (the BTP gains or looses funcional groups, which alters
% the molar weight).

% metabolism and biotransformation may be used interchangable in this
% script

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ke    = par.ke(1);    % elimination rate constant, d-1
ku    = par.ku(1);    % uptake rate constant, L/kg/d

%implementing tranformation rates as in Fu & Rösch 2018
% only phase one metabolism considered
km  = par.km(1);   % biotransformation rate [d-1]
kem = par.kem(1);  % metabolite elimination rate [d-1]


%% Extract correct exposure for THIS time point
% Allow for external concentrations to change over time, either
% continuously, or in steps, or as a static renewal with first-order
% disappearance. For constant exposure, the code in this section is skipped
% (and could also be removed).

if glo.timevar(1) == 1 % if we are warned that we have a time-varying concentration ...
    c = read_scen(-1,c,t,glo); % use read_scen to derive actual exposure concentration
    % Note: this has glo as input to the function to save time!
    % the -1 lets read_scen know we are calling from derivatives (so need one c)
end

%% Calculate the derivatives
% This is the actual model, specified as a system of two ODEs:

% Calculating the change of ligand concentration in membrane protein
% fraction with either (1) the Michaelis-Menten kinetics or (2) the second
% order kinetics. 

% Pre-allocate (helps catch missing assignments in new model options)

    dC_m = km*C_p - kem*C_m; %added another compartment containing the metabolite
    dC_p  = ku*c - ke*C_p - km*C_p; %added biotransformation as a second elimination pathway
    dC_tot = dC_p + dC_m; %total chemical concentration (parent + 1stBTP
    

dX = [dC_p;dC_m;dC_tot]; % collect all derivatives in one vector dX,
