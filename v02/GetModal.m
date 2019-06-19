% GetModal: Generates modal data from spatial mass and stiffness matrices
%
% [Modes,AircraftData,Griddata] = GetModal(MB,KB,Griddata,Zeta,AircraftData)
% 
% This function takes in the spatial domain mass, damping and stiffness 
% matrices, the  and the modal damping ratios and outputs the modal domain
% data.
% 
%
% Inputs:
%
% -MB: A (3*(Np-1))-by-(3*(Np-1)) mass matrix of the aircraft (spatial).
%   The first three rows correspond to the heave, bend and twist states at
%   first node, the next three rows for the second node and so on for Np-1
%   nodes. Here, Np-1 is the number of free nodes.
% 
% -KB: A (3*(Np-1))-by-(3*(Np-1)) stiffness matrix of the aircraft 
%   (spatial). The first three rows correspond to the heave, bend and twist
%   states at first node, the next three rows for the second node and so
%   on for Np-1 nodes. Here, Np-1 is the number of free nodes.
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
%       aircraft(m).
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
%
% -Zeta: An m element vector containing the modal damping coefficients of
%   the modes of the aircraft. Here, m is equal to the number of modes
%   required to be output by the FE model. %
%
% -AircraftData (Optional): A structure containing the high level data
%   about the aircraft. It contains the following fields:
%   AircraftData.mAC: Total mass of the aircraft (Kg).
%   AircraftData.IB: Total inertias of the aircraft (Kg-m^2).
%   AircraftData.s: Wing area of the aircraft (m^2).
%   AircraftData.c: Reference chord of the aircraft (m).
%   AircraftData.cg: Position of the center of gravity of the aircraft (m).
%   AircraftData.b: Reference span of the aircraft (m).
%
%
% Outputs:
% 
% -Modes: A structure containing the matrices describing the final
%   structural model of the aircraft. It contains modal and spatial level
%   matrices like: 
%   Modes.Phif: An (Np-1)*3-by-m matrix where the ith column contains the
%       mode shape of ith mode. The first three rows correspond to the
%       heave, bend and twist at first node, the next three rows for the
%       second node and so on for Np-1 free nodes.
%   Modes.Phib: An (Np-1)*3-by-6 matrix where the columns contains the 
%       rigid mode shape The columns correspond to x, y, z displacement
%       and x, y, z rotation in the body fixed frame while the rows
%       correspond to the heave, bend and pitch of Np-1 free nodes. The
%       first three rows correspond to the heave, bend and twist at first
%       node, the next three rows for the 2nd node and so on for Np-1 free
%       nodes.
%   Modes.Omegaf: An m element matrix containing the calculated modal 
%       frequencies for the first m modes (rad/s).
%   Modes.Zeta: An m element matrix containing the modal damping ratios
%       for the first m modes.
%   Modes.Mff: An m-by-m modal mass matrix of the aircraft. Here, m is 
%       the number of modes considered.
%   Modes.Bff: An m-by-m modal damping matrix of the aircraft. Here, m
%       is the number of modes considered.
%   Modes.Kff: An m-by-mx modal stiffness matrix of the aircraft. Here,
%       m is the number of modes considered.
%   Modes.MB: An (3*(Np-1))-by-(3*(Np-1)) mass matrix of the aircraft 
%       (spatial). The first three rows correspond to the heave, bend and
%       twist states at first node, the next three rows for the 2nd node
%       and so on for Np-1 free nodes. Here, Np-1 is the number of free
%       nodes.
%   Modes.KB: An (3*(Np-1))-by-(3*(Np-1)) stiffness matrix of the aircraft 
%       (spatial). The first three rows correspond to the heave, bend and
%       twist states at first node, the next three rows for the 2nd node
%       and so on for Np-1 free nodes. Here, Np-1 is the number of free
%       nodes.
%
% -AircraftData: A structure containing the high level data about the
%   aircraft. If the optional input, AircraftData is not provided, the
%   output AircraftData has the following fields:
%   AircraftData.mAC: Total mass of the aircraft (Kg).
%   AircraftData.IB: Total inertias of the aircraft (Kg-m^2).
%   AircraftData.cg: Position of the center of gravity of the aircraft (m).
%   If AircraftData is provided as an input, the function
%   overwrites these fields with values obtained from the FE model.
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
%       the coordinates of the ith grid point about the center of gravity
%       of the aircraft(m).
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


function [Modes,AircraftData,Griddata] = GetModal(MB,KB,Griddata,Zeta,varargin)


%% Eigenvalue decomposition of the spatial matrices


[V,D] = eig(KB,MB);     % Eigenvalue decomposition of the spatial matrices
Freq = (diag(sqrt(D))); % Obtaining frequencies (Rad/s)
[FreqSort,idx]=sort(Freq);  % Sorting the frequencies
VSort = V(:,idx);   % Sorting the eigenvalues in same order as frequencies


%% Generating flexible modal data


Phif = VSort(:,4:length(Zeta)+3);           % Flexible modes
Omegaf = FreqSort(4:length(Zeta)+3);        % Modal frequencies
Mff = Phif'*MB*Phif;                        % Modal mass matrix
Kff = Phif'*KB*Phif;                        % Modal stiffness matrix
Bff = 2*diag(Zeta'.*Omegaf)*Mff;            % Modal Damping matrix


%% Changing the origin on the grid to the CG location


% Obtaining CG location
CGloc = GetCGloc(Griddata);

% Changing the origin to the CG location
Griddata.GridP = Griddata.GridP - repmat([0 CGloc],[Griddata.Np 1]);


%% Generating rigid modal data


% Displacement matrix from body fixed coordinates to structural coordinates
d_xb = 1;
d_yb = -1;
d_zb = -1;
dd = [d_xb,d_yb,d_zb];

% Rotation matrix from body fixed coordinates to structural coordinates
th_xb = 90;     %[deg]
th_yb = -90;    %[deg]
th_zb = -90;    %[deg]

Rx = [1 0 0; 0 cosd(th_xb) -sind(th_xb); 0 sind(th_xb) cosd(th_xb)];    
Ry = [cosd(th_yb) 0 sind(th_yb);0 1 0; -sind(th_yb) 0 cosd(th_yb)];
Rz = [cosd(th_zb) -sind(th_zb) 0; sind(th_zb) cosd(th_zb) 0; 0 0 1];

% Creating the rigid body modes
Nf = Griddata.Np-1;     % Number of free gridpoints
Phib = zeros(Nf,6);     % Rigid body mode matrix

% Modified gridpoint matrix with Z-coordinate included
Grid_p2 = [Griddata.GridP(1:Nf,[2 1]) zeros(Nf,1)]';         

% Obtaining rigid modes
for ij = 1:6    
    switch ij                                                   
        case {1,2,3}
            vecrb(:,:,ij) = zeros(Nf,6);
            vecrb(:,ij,ij) = dd(ij)*ones(Nf,1);
            for ii = 5:-1:0
                Phib((1:Nf)*6-ii,ij) = vecrb(:,6-ii,ij);
            end
        case 4
            vecrb(:,:,ij) = zeros(Nf,6);
            vecrb(:,1:3,ij) = (Rx*Grid_p2)'-[Grid_p2(1,:)' zeros(Nf,1) zeros(Nf,1)];
            vecrb(:,ij,ij) = sind(th_xb)*ones(Nf,1);
            for ii = 5:-1:0
                Phib((1:Nf)*6-ii,ij) = vecrb(:,6-ii,ij);
            end
        case 5
            vecrb(:,:,ij) = zeros(Nf,6);
            vecrb(:,1:3,ij) = (Ry*Grid_p2)'-[zeros(Nf,1) Grid_p2(2,:)' zeros(Nf,1)];
            vecrb(:,ij,ij) = sind(th_yb)*ones(Nf,1);
            for ii = 5:-1:0
                Phib((1:Nf)*6-ii,ij) = vecrb(:,6-ii,ij);
            end
        case 6
            vecrb(:,:,ij) = zeros(Nf,6);
            vecrb(:,1:3,ij) = (Rz*Grid_p2)'-[zeros(Nf,1) zeros(Nf,1)  Grid_p2(3,:)'];
            vecrb(:,ij,ij) = sind(th_zb)*ones(Nf,1);
            for ii = 5:-1:0
                Phib((1:Nf)*6-ii,ij) = vecrb(:,6-ii,ij);
            end
    end
    
end

% Index to get the required mode shapes only
idx = [3:6:(Nf)*6;4:6:(Nf)*6;5:6:(Nf)*6];

% Final rigid mode shape vector
Phib = Phib(idx(:),:);                                          


%% Output modal data


Modes.Phif = Phif;              % Flexible Modes
Modes.Phib = Phib;              % Rigid modes
Modes.Omegaf = Omegaf;          % Modal frequencies
Modes.Zeta = Zeta;              % Modal Damping
Modes.Mff = Mff;                % Modal Mass Matrix
Modes.Bff = Bff;                % Modal Damping Matrix
Modes.Kff = Kff;                % Modal Stiffness Matrix
Modes.MB = MB;                  % Spatial Mass Matrix
Modes.KB = KB;                  % Spatial Stiffness Matrix


%% Output aircraft data


% Check if Aircraft_data is provided as input 
if length(varargin)==1
    AircraftData=varargin{1};
end


% Calculate the total mass of the aircraft and update it in Aircraft_data
AircraftData.mAC = sum(sum(MB((1:Griddata.Np-1)*3-2,(1:Griddata.Np-1)*3-2)));


% Calculate the inertia matrix of the aircraft
Mr = Phib'*MB*Phib;
Izb = Mr(4,4);      % Roll inertia   
Ixb = Mr(5,5);      % Pitch inertia
Iyb = Izb + Ixb;    % Yaw inertia

AircraftData.IB = diag([Izb,Ixb,Iyb]);    % Inertia matrix


% CG location of the aircraft

AircraftData.cg = CGloc;
end
