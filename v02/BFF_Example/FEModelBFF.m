% The following code generates a finite element based structural model of
% the BFF Aircraft, given a set of parameter values. The parameters 
% represent the structural properties of various parts of the aircraft and
% the mass / inertia distribution. The list of parameters for BFF FE model
% can be found in the reference.
%
% For a given set of parameters, the code generates a FE based structural
% grid of the aircraft as shown in Figure 4 in the reference. It, then, 
% generates a spatial domain mass and stiffness matrices of corresponding
% to the heave, bend and twist states of the gridpoints. These spatial 
% domain matrices are used to generate modal domain data including modal
% frequencies and mode shapes. A state space description, with force as an 
% input and acceleration as output at the gridpoints is also generated.
%
% Reference:
%
% Gupta, A., Moreno, C. P., Pfifer, H., Balas, G. J.,  "Updating a Finite
% Element Based Structural Model of a Small Flexible Aircraft", AIAA
% (Available at: 'http://www.aem.umn.edu/~AeroServoElastic/Papers/2015/
% GuptaEtAl_15AIAA_UpdatingAFiniteElementBasedStructuralModelOfSmallFlexAC
% .pdf')


% Defining optimal parameter values
Params.MuFl = 0.6;
Params.E = 56;
Params.G = 234;
Params.MuW = 0.62;
Params.Mass=[0.10, 0.00, 0.60, 0.20, 0.00, 0.14, 0.19, 0.08];
Params.Inertia=[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000;
                0.060, 0.000, 0.000, 0.000, 0.006, 0.000, 0.002, 0.002]';

         
% Defining the modal damping coefficients 
% The first 6 are calculated from the GVT, the others are default values
Zeta = [1.55 1.06 2.06 2.33 2.85 2.55 5*ones(1,6)]/100;


% Number of modes
NumModes = length(Zeta);


% Load aircraft data structure
load AircraftData_v02


%% Generating the structural model of BFF


% Obtain griddata 
Griddata = BFFParameters2FEM(Params);


% Obtain spatial matrices and state space description
[MB,CB,KB,FEModelSS] = GetSpatial(Griddata,Zeta);


% Obtain modes and global aircraft mass and inertia from the FE Model
[Modes,AircraftData] = GetModal(MB,KB,Griddata,Zeta,AircraftData);

