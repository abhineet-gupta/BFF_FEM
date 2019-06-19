% GetCGloc: Calculates the CG location of the aircraft
%
% [CGloc] = GetCGloc(Griddata)
% 
% This function calculates the center of gravity of the aircraft for a
% given finite element griddata.
%
%
% Input:
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
%
% Output:
% 
% -CGloc: A variable which specifies the z coordinate of the centre of
%   gravity of the aircraft (m).


function [CGloc] = GetCGloc(Griddata)


%% Defining variables


mass_moment_x = 0;      % Variable for moment of mass in x direction
mass_moment_z = 0;      % Variable for moment of mass in z direction
mass_total = 0;         % Variable for total mass


%% Importing gridata


Np = Griddata.Np;               % Number of gridpoints
Nel = Griddata.Nel;             % Number of elements
GridP = Griddata.GridP;         % Gridpoint coordinates
MatProp = Griddata.MatProp;     % Material property matrix
El = Griddata.El;               % Element information matrix
MPt = Griddata.MPt;             % Point mass matrix



%% For distributed mass


for el_num = 1:Nel
    % Length of each element
    lengthi = norm(GridP(El(el_num,1),:)-GridP(El(el_num,2),:));      
    
    % Mid point location of each element
    locationi = (GridP(El(el_num,1),:) + GridP(El(el_num,2),:))/2;    
    
    % Mass per unit length of each element
    mu = MatProp(El(el_num,3),5);                                      
    
    
    % Summing up the mass moment in x and z direction
    mass_moment_x = mass_moment_x + (lengthi * locationi(1) * mu);          
    mass_moment_z = mass_moment_z + (lengthi * locationi(2) * mu);      
    
    % Summing up the total mass
    mass_total = mass_total + (lengthi * mu);                           
end

%% For point masses


for gridp=1:Np
    % Location of each gridpoin
    locationi = GridP(gridp,:);
    
    % Number of elements which share the gridpoint
    num_shared = length(find(El(:,1:2)==gridp));        
   
    
    % Summing up the mass moment in x and z direction
    mass_moment_x = mass_moment_x + (locationi(1) * num_shared * MPt(gridp));          
    mass_moment_z = mass_moment_z + (locationi(2) * num_shared * MPt(gridp));          
    
    % Summing up the total mass
    mass_total = mass_total + (num_shared * MPt(gridp));                               
end


%% Output the CG location


% [x,z] coordinates of CG location
CGloc_temp=[mass_moment_x mass_moment_z]/mass_total;        

% Z coordinate of CG location
CGloc = CGloc_temp(2);                                      


end