% GetSpatial: Generates the spatial domain mass and stiffness matrices
%
% [MB,CB,KB,FEModelSS] = GetSpatial(Griddata,Zeta)
% 
% This function takes in the FE description of the aircraft structure and 
% generates the spatial level mass, damping and stiffness matrices
% corresponding the to heave, bend and twist states of the gridpoints and
% the structural model in state space format.
% 
%
% Inputs: 
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
%
% -Zeta: An m element vector containing the modal damping coefficients of
%   the modes of the aircraft. Here, m is equal to the number of modes
%   required to be output by the FE model. 
%
%
% Outputs:
%
% -MB: A (3*(Np-1))-by-(3*(Np-1)) mass matrix of the aircraft (spatial).
%   The first three rows correspond to the heave, bend and twist states at
%   first node, the next three rows for the second node and so on for Np-1
%   nodes. Here, Np-1 is the number of free nodes.
%
% -CB: A (3*(Np-1))-by-(3*(Np-1)) damping matrix of the aircraft (spatial).
%   The first three rows correspond to the heave, bend and twist states at 
%   first node, the next three rows for the second node and so on for Np-1
%   nodes. Here, Np-1 is the number of free nodes.
% 
% -KB: A (3*(Np-1))-by-(3*(Np-1)) stiffness matrix of the aircraft 
%   (spatial). The first three rows correspond to the heave, bend and twist
%   states at first node, the next three rows for the second node and so
%   on for Np-1 nodes. Here, Np-1 is the number of free nodes.
%
% -FEModelSS: A state space system representation of the structural model
%   of the aircraft. It has 3*(Np-1) inputs: corresponding to force,
%   bending moment and twisting moments at (Np-1) free points, 3*(Np-1)
%   outputs corresponding to acceleration in heave, bend and twist at
%   (Np-1) free points and 6*(Np-1) states corresponding to heave, bend,
%   twist and their corresponding velocities at (Np-1) points. Here, Np-1
%   is the number of free nodes.


function [MB,CB,KB,FEModelSS] = GetSpatial(Griddata,Zeta)


%% Obtaining FE model description from Griddata


Np = Griddata.Np;               % Number of gridpoints
Nel = Griddata.Nel;             % Number of elements
Cons = Griddata.Cons;           % Constraint matrix
MatProp = Griddata.MatProp;     % Material property matrix
El = Griddata.El;               % Element information matrix
MPt = Griddata.MPt;             % Point mass matrix
IPt = Griddata.IPt;             % Point inertia matrix
GridP = Griddata.GridP;         % Gridpoint coordinates


%% Calculating Number of equations


ne=Cons;  nec=0;
for i=1:Np
   for j=1:3
      if Cons(i,j)==0       % Check for existance of boundary condition
         nec=nec+1;
         ne(i,j)=nec+1;     % Save number of equation in array ne
      end
   end
end
nec1=nec+1;                 % Index 1 is going to be for trash	


%% Obtaining global mass and inertia matrix


% Initiating Array
KK=sparse(nec1,nec1);           % Global stiffness matrix
MM=sparse(nec1,nec1);           % Global mass matrix


% Constructing the global matrices in sparce form
for ik=1:Nel
    i00=El(ik,1:2);                     % Nodes 1-2 for element
    mi=[ne(i00(1),:) ne(i00(2),:)];     % Local incidence in global matrix
    E=MatProp(El(ik,3),1);				% Elastic constant (Young Modulus)
    G=MatProp(El(ik,3),2);      % Torsion constant (Modulus of Rigidity)
    b1=MatProp(El(ik,3),3);     % Section width
    h1=MatProp(El(ik,3),4);     % Section height
    mas=MatProp(El(ik,3),5);    % Mass per length
    
    % Obtaining element inertias
    Iz=1/12*b1*(h1^3);          % Cross section inertia
    Iy=1/12*h1*(b1^3);                          
    J=Iy+Iz;
    xi=1/12*mas*(b1^2+h1^2);    % Inertia per length
    xn=GridP(i00,1);			% Element coordinates
    zn=GridP(i00,2); 
   
   
    % Creating individual element matrices
    mv=MPt(i00,1);      % Temporary variable holding point masses
    iv1=IPt(i00,1);     % Temporary variable holding point inertias
    iv2=IPt(i00,2);
    le=sqrt((xn(2,1)-xn(1,1))^2+(zn(2,1)-zn(1,1))^2);   % Length of element
    cose=(xn(2,1)-xn(1,1))/le;
    sen=(zn(2,1)-zn(1,1))/le;
    tao=[1 0 0;0 cose -sen; 0 sen cose];    % Coordinate transformation 
    t = blkdiag(tao,tao);

    % Defining elemental stiffness matrix
    ke=[12*E*Iz/(le^3) 6*E*Iz/(le^2) 0 -12*E*Iz/(le^3) 6*E*Iz/(le^2) 0;
        6*E*Iz/(le^2) 4*E*Iz/le 0 -6*E*Iz/(le^2) 2*E*Iz/le 0;
        0 0 G*J/le 0 0 -G*J/le; 
        -12*E*Iz/(le^3) -6*E*Iz/(le^2) 0 12*E*Iz/(le^3) -6*E*Iz/(le^2) 0;
        6*E*Iz/(le^2) 2*E*Iz/le 0 -6*E*Iz/(le^2) 4*E*Iz/le 0;
        0 0 -G*J/le 0 0 G*J/le];
   
    % Defining temporary variables to be used in mass matrices
    aa=mas*le/420; bb =0; cc=xi*le/6;            
    maddt=zeros(6); maddt([1 4],[1 4])=diag(mv); 
    maddt([2 5],[2 5])=diag(iv1); maddt([3 6],[3 6])=diag(iv2);
    
    % Defining elemental mass matrix
    me=[156*aa 22*le*aa -7*bb/2 54*aa -13*le*aa -3*bb/2;
        22*le*aa 4*le^2*aa -le*bb/2 13*le*aa -3*le^2*aa -le*bb/3;
        -7*bb/2 -le*bb/2 2*cc -3*bb/2 le*bb/3 cc;
        54*aa 13*le*aa -3*bb/2 156*aa -22*le*aa -7*bb/2;
        -13*le*aa -3*le^2*aa le*bb/3 -22*le*aa 4*le^2*aa le*bb/2;
        -3*bb/2 -le*bb/3 cc -7*bb/2 le*bb/2 2*cc];
   
    % Assembling to generate global matrices
    KK(mi,mi)=KK(mi,mi)+t'*ke*t;
    MM(mi,mi)=MM(mi,mi)+t'*me*t+maddt;
end


% Number of degrees of freedom
nec2=2:nec1;


%% Generating final matrices


% Stiffness Matrix
KB = full(KK(nec2,nec2));

% Mass Matrix
MB = full(MM(nec2,nec2));


% Damping Matrix

% Eigenvalue decomposition of the spatial matrices
[V,D] = eig(KB,MB);         
Freq = (diag(sqrt(D)));     % Obtaining frequencies (Rad/s)
[FreqSort,idx]=sort(Freq);  % Sorting the frequencies
VSort = V(:,idx);           % Sorting the eigenvectors in same order

% Add additional damping coefficients 
Zeta = [Zeta , repmat(0.05,1,size(MB,1)-length(Zeta))];

% Generating global damping matrix
mff = VSort'*MB*VSort;                     % Modal mass matrix
bff = 2*diag(Zeta'.*FreqSort)*mff;         % Modal damping matrix
CB = VSort'\bff/VSort;                     % Global damping matrix


%% Obtaining state space model


Nst = size(MB);                           
A = [zeros(Nst) eye(Nst); -MB\KB -MB\CB];           % A matrix
B = [zeros(Nst); inv(MB)];                          % B matrix
C = [-MB\KB -MB\CB];                                % C matrix
D = inv(MB);                                        % D matrix

FEModelSS=ss(A,B,C,D);      % Structural model in state space format


end
