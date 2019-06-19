% BFFParameters2FEM - Generates the FE model description of the BFF 
%   Aircraft
%
% [Griddata] = BFFParameters2FEM(Params)
%
% This function takes in specific parameter values of the FE model
% and returns a complete (geometric and material property) description
% of the FE model for the Body Freedom Flutter aircraft.
% 
%
% Inputs:
%
% -Params: Structure of parameter values used to model / optimize the BFF
%   Aircraft. It has following fields:
%   Params.MuFl: Mass per unit length of fuselage elements (Kg/m).
%   Params.E: Young's Modulus of the elements constituting the wings (Gpa).
%   Params.G: Torsional Modulus of the elements constituting the
%       wings (Mpa).
%   Params.MuW: Mass per unit length of wing elements (Kg/m).
%   Params.Mass: An eight element vector representing the point masses at 
%       gridpoints 1-8. As the mass distribution on the aircraft is 
%       assumed to be symmetric, the point masses at other gridpoints are 
%       automatically known from this vector. For example: Point mass 
%       at gridpoint-13 is equal to point mass at gridpoint-8
%       (Kg).
%   Params.Inertia: An 8-by-2 element matrix where each row represents the
%       point [roll inertia, pitch inertia] at gridpoints 1-8. As the mass
%       distribution on the aircraft is assumed to be symmetric, the
%       point inertias at other gridpoints are automatically known 
%       from this vector. For example: Point inertia at gridpoint-13 
%       is equal to point pitch inertia at gridpoint-8 (Kg-m^2).
% 
%
% Outputs:
% 
% -Griddata: Structure containing matrices that describe the geometry and
%   material properties of the FE model, describing the FE model 
%   completely. The fields of Griddata are:
%   Griddata.Np: Number of gridpoints in the FE model.
%   Griddata.Nel: Number of elements in the FE model.
%   Griddata.Cons: Np-by-3 matrix where the ith row represents the
%       presence of constraints on heave, bend and roll of the ith
%       gridpoint respectively. Presence of constraint is denoted by a 1
%       and absence by 0.
%   Griddata.GridP: Np-by-2 matrix where the ith row is [X, Y]. (X,Y) is
%       the coordinates of the ith grid point about the nose of the
%       aircraft (m).
%   Griddata.MatProp: Ns-by-5 matrix where the ith row contains the
%       section properties: Young's modulus, rigidity modulus, width,
%       height, mass/length of a particular material. Here, Ns is the
%       number of different materials used in the FE model.
%   Griddata.El: Nel-by-3 matrix where the ith row is [P1, P2, Index]. P1
%       and P2 are the end gridpoints of the ith element. Index specifies 
%       the material properties of the element by specifying the row 
%       number of MatProp matrix corresponding to the material of the 
%       ith element.
%   Griddata.MPt: Np-by-1 vector where the ith element contains the 
%       additional point mass at the ith gridpoint of the FE model (Kg).
%   Griddata.IPt: Np-by-2 matrix where ith row is [Iz, Ix]. Iz is point
%       roll inertia at ith gridpoint and Ix is point pitch inertia at ith 
%       gridpoint (Kg-m^2).


function [Griddata] = BFFParameters2FEM(Params)
% To simulate free-free condition, it is assumed that the finite element
% model is restricted in motion at the last gridpoint. 


%% Defining the gridpoints


% Number of total gridpoints (including the constrained one)
Np = 18;


% Number of elements
Nel = 21;


% Boundary conditions (Assign 1 if degree of freedom is restrained, 0 
% if not)
Cons = zeros(Np,3);  
Cons(Np,1:3) = 1;       % The last gridpoint is restricted in motion


% Defining gridpoint locations with respect to point of intersection of 
% spars and center chord
fac = 2.54/100;         % Conversion factor from inches to m
GridP = zeros(Np,2);

% Center body
GridP(1,:) = [0,0.2005];
GridP(2,:) = [0,0.4151];
GridP(3,:) = [0,-0.1392];

% Left wing 
GridP(4,:) = 15.65*[cosd(22),sind(22)]*fac;
GridP(5,:) = 20.4*[cosd(22), sind(22)]*fac;
GridP(6,:) = 35.65*[cosd(22), sind(22)]*fac;
GridP(7,:) = 50.03*[cosd(22), sind(22)]*fac;
GridP(8,:) = 65.65*[cosd(22), sind(22)]*fac;
GridP(14,:) = GridP(4,:) - [0 0.1016];
GridP(15,:) = GridP(4,:) + [0 0.0635];

% Mirror left wing to get grid points on the right wing
GridP(9:13,:) = [-GridP(4:8,1),GridP(4:8,2)];
GridP(16:17,:) = [-GridP(14:15,1),GridP(14:15,2)];


% Constrained gridpoint
GridP(Np,:) = [0,1];

% Shifting origin to the nose of the aircraft
GridP = GridP - repmat([0,-0.3903],[Np 1]);


%% Defining the structural elements


% Material properties matrix
% Each row contains the structural properties of a specific material used
% to model the aircraft. A row contains the Young's modulus, rigidity 
% modulus, width, height, mass/length a particular material. 
MatProp = [5e15, 0.6e15, 0.0792, 0.0106, Params.MuFl;
   Params.E*1e9, Params.G*1e6, 0.0762, 0.0065, Params.MuW;
   1e-6, 1e-4, 1e-5, 1e-5, 1e-2;
   5e15, 0.6e15, 0.08, 0.01, 0.01]; % [N/m^2 N/m^2 m m kg/m]


% Element information matrix
% Nel-by-3 matrix where the ith row is [P1, P2, Index]. P1 and P2 are the 
% end gridpoints of the ith element. Index specifies the material
% properties of the element by specifying the row number of Mat_prop 
% matrix corresponding to the material of the ith element.
El = [
    1 4 1
    4 5 2;
    5 6 2;
    6 7 2;
    7 8 2;
    1 9 1;
    9 10 2;
    10 11 2;
    11 12 2;
    12 13 2;
    1 3 1;
    1 2 1;
    1 Np 3
    3 9 1
    3 4 1
    2 9 1
    2 4 1
    17 9 4
    9 16 4
    4 14 4
    15 4 4
    ];


% Determine # of elements associated with each grid point
ElCount = zeros(Np,1);
for i=1:Np
    idx = find(El(:,1:2)==i);
    ElCount(i) = numel(idx);
end


%% Defining mass distribution 


% Divide the point mass and point inertias by the number of elements
% sharing the grid points. Note that symmetry is assumed, so that the
% point masses and point inertias are only defined for one wing
for i=1:size(Params.Mass,2)
   Params.Mass(i) = Params.Mass(i)/ElCount(i);
   Params.Inertia(i,:) = Params.Inertia(i,:)/ElCount(i);   
end


% Additional point mass matrix.
% This matrix contain the value of point mass at a gridpoint divided by 
% the number of elements sharing the grid point.
MPt = zeros(Np,1);

MPt(1:3) = Params.Mass(1:3);            % Center Body
MPt(4:8) = Params.Mass(4:8);            % Left Wing
MPt(9:13)= Params.Mass(4:8);            % Right Wing


% Additional point inertias matrix
% This matrix contains the additional point inertia [Iz, Ix] at a
% gridpoint divided by the number of elements sharing the gridpoint.
IPt=zeros(Np,2);

IPt(1:3,:) = Params.Inertia(1:3,:);     % Center Body
IPt(4:8,:) = Params.Inertia(4:8,:);     % Left Wing
IPt(9:13,:)= Params.Inertia(4:8,:);     % Right Wing


%% Final Output


Griddata.Np = Np;                   % Number of gridpoints
Griddata.Nel = Nel;                 % Number of elements
Griddata.Cons = Cons;               % Constraint matrix
Griddata.GridP = GridP;             % Gridpoint coordinates
Griddata.MatProp = MatProp;         % Material property matrix
Griddata.El = El;                   % Element information matrix
Griddata.MPt = MPt;                 % Point mass matrix
Griddata.IPt = IPt;                 % Point inertia matrix


end