clear all

%This code is for Question 1, Part A, Part ii. Which is 4 elements with the
%left hand support being fixed and the right hand being pinned.

numberElements=4; %4 elements to analyse for
numberNodes=numberElements+1; %nodes=elements+1
dofNode=3; %the degrees of freedom per node is 3 for a frame analysis
nodesEl=2; %nodes per element is 2; i and j
dofEl=dofNode*nodesEl; %Degree of freedom per element is dof per node * nodes per element

%Discretisation process

elementNodes=[1 2; 2 3; 3 4; 4 5]; % the global noding system
nodeCoordinates=[0 0; 0 15; 0 30; 60 30; 60 15]; %Coordinates for the node positions
bcDof=[1,2,3,13,14]; %NOTE node 15 is missing as a pinned support is present, meaning rotation is ALLOWED

A=200; %cm^2 units
E=20*10^6 %N/cm^2 units
I=10^4*6; %cm^4 units

P=[0;0;0;150000;0;0;0;0;0;0;0;0;0;0;0]; %Horizontal force applied at node two

xx=nodeCoordinates(:,1); %Extracting the x coordinates
yy=nodeCoordinates(:,2); %Extracting the y coordinates
sysDof=dofNode*numberNodes; %Degrees of freedom for the whole system
noConstraints=length(bcDof); %The number of constraints used in a later loop is equal to the number of boundary conditions
elDofs=zeros(numberElements,dofEl); %Sizing a matrix to store the element degrees of freedom

K=zeros(sysDof,sysDof); %Sizing the global stiffness matrix


for e=1:numberElements
    elDofs(e,1)=(elementNodes(e,1)*3)-2; % dof U3i-2
    elDofs(e,2)=(elementNodes(e,1)*3)-1; %dof U3i-1
    elDofs(e,3)=(elementNodes(e,1)*3); %dof U3i
    elDofs(e,4)=(elementNodes(e,2)*3)-2; %dof U3j-2
    elDofs(e,5)=(elementNodes(e,2)*3)-1; %dof U3j-1
    elDofs(e,6)=(elementNodes(e,2)*3); %dof U3j
    index1=elementNodes(e,1); %Calls node i DOFs
    index2=elementNodes(e,2); %Calls node j DOFs
    xx1=xx(index1); %x_i coordinate
    xx2=xx(index2);%x_j coordinate
    yy1=yy(index1); %y_i coordinate
    yy2=yy(index2); %y_j coordinate
    L=((xx2-xx1)^2+(yy2-yy1)^2)^0.5; %Length of element e calculation
    c=(xx2-xx1)/L; %cos(pheta) equal
    s=(yy2-yy1)/L; %sin(pheta) equal
    ALI=A*L^2/I; %A and I constants for every element, length is specific
    EIL=E*I/L^3; %E and I constants for every element, length is specific
    
    
    Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1]; %Transformation matrix
    kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2]; %element stiffness before transformation
    ke=transpose(Te)*kebar*Te; %tranformed element stiffness
    
    K(3*e-2:3*e+3,3*e-2:3*e+3)=K(3*e-2:3*e+3,3*e-2:3*e+3)+ke; %Assembly procedure
end

%Setting the rows and columns to 0 for specific boundary conditions

for c=1:noConstraints
    i=bcDof(c);
    j=bcDof(c);
    K(i,:)=0;
    K(:,j)=0;
    P(i)=0;
end


activeDof=setdiff([1:sysDof],[bcDof]);  %produce a list of active DOFs
Kr = K(activeDof,activeDof);    %produce a reduced stiffness matrix
Fr = P(activeDof);  


U=Kr\Fr; %Calculating the displacements and rotations present

%Resizing the displacement and rotation matrix

uFinal=zeros(3*4+3,1)
uFinal(4:3*4,1)=U(1:size(U)-1,1);
uFinal(3*4+3,1)=U(3*4-2,1);

%POST PROCESSING

nodalForce = K*uFinal; %Forces and moments in each node
elementForce = zeros(numberElements,dofEl); %Sizing the element forces for each degree of freedom for every element

for e = 1:numberElements
    elDofs(e,1) = (elementNodes(e,1)*3)-2;      %dof U_3i-2
    elDofs(e,2) = (elementNodes(e,1)*3)-1;      %dof U_3i-1
    elDofs(e,3) = (elementNodes(e,1)*3);        %dof U_3i
    elDofs(e,4) = (elementNodes(e,2)*3)-2;      %dof U_3j-2
    elDofs(e,5) = (elementNodes(e,2)*3)-1;      %dof U_3j-1
    elDofs(e,6) = (elementNodes(e,2)*3);        %dof U_3j
    index1 = elementNodes(e,1);                 %calls node i DOFs
    index2 = elementNodes(e,2);                 %calls node j DOFs
    xx1 = xx(index1);                           %x_i coordinates
    xx2 = xx(index2);                           %x_j coordinates
    yy1 = yy(index1);                           %y_i coordinates
    yy2 = yy(index2);                           %y_j coordinates
    L = (((xx2-xx1)^2)+((yy2-yy1)^2))^(0.5); %element length calculation
    C(e) = (xx2-xx1)/L;                      %element cosine calculation
    S(e) = (yy2-yy1)/L;                      %element sine calculation
    c = C(e);                                   %renamed for convenience
    s = S(e);                                   %renamed for convenience
    cs = c*s;                                   %definition
    c2 = c^2;                                   %definition
    s2 = s^2;                                   %definition
    ALI=A*L^2/I;
    EIL=E*I/L^3;                               %definition


    Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1]; %Transformation matrix of a frame
    kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2]; %element stiffness matrix NOT transformed
    ke=transpose(Te)*kebar*Te %Transformed element stiffness matrix



    gDof = elDofs(e,:);
    cDofStart = gDof(1);
    cDofEnd = gDof(end);
    uCurrent = uFinal(cDofStart:cDofEnd);
    elementForce(e,:) = ke*uCurrent; %Element forces calculated for its' 6 DoF
end

Xdisp=zeros(size(xx)); %Sizing a matrix specifically for x displacements
Ydisp=zeros(size(yy)); %Sizing a matrix specifically for y displacements

Xdisp(2:numberNodes-1,1)=U(1:3:(numberNodes-2)*3,1); %Storing data in appropriate locations
Ydisp(2:numberNodes-1,1)=U(2:3:(numberNodes-2)*3+1,1); %Storing data in appropriate locations



%Plotting and labelling the orginal shape
plot(xx,yy)
hold on
axis([-10 max(xx)+10 -10 max(yy)+10])
xlabel('Horizontal Distance (cm)');
ylabel('Vertical Distance (cm)');
title('Orginal Shape');
hold off

%Plotting and labbeling the deformed shape
Xfinal=nodeCoordinates(:,1)+Xdisp(:,1); %The final x coordinate positions after loading
Yfinal=nodeCoordinates(:,2)+Ydisp(:,1); %The final y coordiante positions after loading
figure
plot(Xfinal, Yfinal)
hold on
axis([-10 max(Xfinal)+10 -10 max(Yfinal)+10])
xlabel('Horizontal Distance (cm)');
ylabel('Vertical Distance (cm)');
title('Deformed Shape');
hold off
    