%% BYOM, byom_calanus_2016_onecomp.m 
%
% * Author: Tjalling Jager 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* An organism is exposed to a chemical in its surrounding
% medium. The animal accumulates the chemical according to standard
% one-compartment first-order kinetics.  
%
% *This script:* One-compartment TK of C2-naphthalene in _Calanus finmarchicus_.
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 
%
%
% September 2022 
% Modifications by AMD to include (irrevirsible) receptor binding of an 
% antagonist. Here with the example of the antagonist thiacloprid (THI) 
% binding to the nicotinic-acetylcholine receptor (nAChR) as observed 
% in Gammarus pulex. in Raths et al. 2023

% January 2026
% Modifications by JR to include biotransformation. Receptor binding was
% removed from this version of the script.


%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.
clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
% set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 1; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood



%%% test run with dataset G_CIT_21 (G.pulex, citalopram, 21°C) from Raths et al. (2023, speed it up)

% Internal concentrations of parent compound in target organism in [µmol/kg] 
% equals input data for state 1.
DATA{1} = [0.5    1	1   1
0	0	0	0
0.025	0.510	0.700	0.550
0.050	0.920	0.640	0.700
0.100	0.670	0.600	0.660
0.200	0.610	0.650	0.750
0.400	0.970	0.590	0.550
0.800	0.630	0.780	0.650
1.400	0.460	0.430	0.430
1.425	0.080	0.030	0.110
1.450	0.050	0.050	0.050
1.500	0.010	0.010	0.010
1.600	0.000	0.000	0.000

];

% Internal concentrations of the primary transformation product (phase one metabolite) in [µmol/kg] 
% equals input data for state 2
% !!! important: comment this vector for state 2 out if no % biotransformation is used (a vector of only zeros will cause an error)
% DATA{2} = [0.50    1	1   1
% 0       0       0       0
% 0.25	0.00	0.00	0.00
% 0.50	0.00	0.00	0.00
% 1.00	0.00	0.00	0.00
% 2.00	0.00	0.00	0.00
% 4.00	0.00	0.00	0.00
% 8.00	0.00	0.00	0.00
% 14.00	0.00	0.00	0.00
% 14.25	0.00	0.00	0.00
% 14.50	0.00	0.00	0.00
% 15.00	0.00	0.00	0.00
% 16.00	0.00	0.00	0.00
% 18.00	0.00	0.00	0.00
% 22.00	0.00	0.00	0.00
% 28      0.00	0.00	0.00
% ];

% If needed, weights for individual measurements can be defined
% For this, uncommend the following line and specify your weights

% W{1} = 21 * ones(size(DATA{1})-1); % i.e. each point is a pooled sample of 21 animals

% In this data set, exposure was time-varying and reported as a series of
% concentrations over time. Here, the scenario is used as a linear forcing
% series (which has an analytical solution, and is thus much faster than
% the ODE version). Double time entries can be used, to account for replicate and sampling variancs.
% Unit: [µmol/L]
Cw1 = [0.5        1
0	19.7
0.2	17.2
0.4	14.7
0.8	14.2
1.4	13.5
1.4001	0
1.6	0

];

make_scen(4,Cw1); % prepare as linear-forcing function interpolation (can use glo.use_ode = 0)  

% Create a table with nicer labels for the legends
Scenario = [1]; 
Label = {'Clomazone concentration [mg/kg]'};
glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = Scenario; % scenarios (concentrations or identifiers)
X0mat(2,:) = 0;      % initial values state 1 (parent internal concentrations)
X0mat(3,:) = 0;      % initial values state 2 (1st BTP internal concentration)
X0mat(4,:) = 0;      % initial values state 3 (total chemical internal concentration)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval scale];
par.ke          = [60           1   0.01        200   2];  % elimination rate constant, d-1
par.ku          = [3         1  0.01        10    2];  % uptake rate constant, L/kg/d

par.km          = [0     0   0           10       1];  % 1stBTP biotransformation rate constant [d-1]
par.kem         = [0     0   0           10       1];  % 1stBTP elimination rate constant [d-1]


%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = ['Concentration parent (mg/kg)'];
glo.ylab{2} = ['Concentration metabolite (mg/kg)'];
glo.ylab{3} = ['Total chemical concentration (mg/kg)'];

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.
% par_out = calc_optim(par,opt_optim); % start the optimisation
% calc_and_plot(par_out,opt_plot); % call the plotting routine again to plot fits with CIs

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 

opt_optim.type = 4; % optimisation method 1) simplex, 4) parameter-space explorer
opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 0; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.bw    = 0; % plot in black and white
opt_plot.cn    = 0; % if set to 1, connect model line to points (only for bw=1)
opt_plot.annot = 1; % annotations in sub-plot: text box with parameter estimates or overall legend
glo.useode     = 1; % use the simplified analytical solution in simplefun.m (0) or the ODE solution in derivatives (1)
glo.stiff      = 1; % 0 = use ode45 with very strict tolerances (default, non-stiff mode); 1 = use ode15s for experiments with long time (i.e. more than one week) or imbalanced paramters

glo.R_mod = 1;  % default for print

opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
% no plotting here; we'll immediately plot with CIs below

%% Plot results with confidence intervals
% The following code can be used to make a standard plot (the same as for
% the fits), but with confidence intervals. Options for confidence bounds
% on model curves can be set using opt_conf (see prelim_checks).
% 
% Use opt_conf.type to tell calc_conf which sample to use: 
% -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% the model curves from the sample 
% 2) parameter sets from a joint likelihood region using the shooting 
% method (limited sets can be used), which will yield (asymptotically) 95% 
% CIs on predictions
% 3) as option 2, but using the parameter-space explorer

opt_conf.type    = 3; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 0; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time

out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs
