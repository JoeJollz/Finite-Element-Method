clear all

%This code is for Question 1, Part B, Part ii. Which is 40 elements with
%the left hand support being FIXED, and the right hand support being PINNED

n=40; %Number of elements present
numberElements=n; %Variable names present to make clearer for the reader

%Discretisation process begins
a=floor((numberElements+3)/6)*2+1; %Number of node present in the left column, including the support and the top.
c=floor(((numberElements+3-a))/2); %Number of nodes present in the right column, including the support and the top nodes.
b=numberElements+3-a-c; %Number of nodes in the top beam, including both the top right and top left corners

numberNodes=numberElements+1;%automatic
dofNode=3; %degrees of freedom per node of a frame is 3
nodesEl=2; %nodes per element is 2; i and j
dofEl=dofNode*nodesEl; %degrees of freedom per element is dof per node * number of nodes per element

elementNodes=zeros(numberElements,2); %creating zeros ready for global node input


for e=1:numberElements   %defining all the global nodes
    elementNodes(e,1)=e;   %x coordinate node 
    elementNodes(e,2)=e+1;  %y coordinate node
end

%nodeCoordinates

nodeCoordinates=zeros(n+1,2);

nodeCoordinates(1:a,2)=0:30/(a-1):30;
nodeCoordinates(a+1:n+1-c,2)=30;
nodeCoordinates(a:n+2-c,1)=0:60/(b-1):60;
nodeCoordinates(n-c+2:n+1,2)=30:-15/(c-1):15;
nodeCoordinates(n-c+2:n+1,1)=60;
%Discretisation process ends. Check the coordinates in the nodeCoordinates
%dateset for optimal distribution between the two columns and single beam.

%Boundary conditions will be the 3 degrees of freedom for the first node,
%and the 2 degrees of freedom at the final node, as the right support is 
%now PINNED and therefore ALLOWING ROTATION
bcDof=[numberElements/numberElements;numberElements/numberElements+1;numberElements/numberElements+2;numberNodes*3-2;numberNodes*3-1];

A=200; %cm^2 units
E=20*10^6; %N/cm^2 units
I=(6*10^4); %cm^4 units

P=zeros(numberNodes*3,1); %sizing the externally applied forces vector
P(((a+1)/2)*3-2,1)=150000; %Applying the 150kN horizontal force at the middle node in the left column. Located at (0,15).

xx=nodeCoordinates(:,1); %Extracting and storing the x coordinate data
yy=nodeCoordinates(:,2); %Extracting and storing the y coorindate data
sysDof = dofNode*numberNodes; %calculating the total number of degrees of freedom in the whole system
noConstraints = length(bcDof); %number of constraints is equal to the number of boundary condtions
elDofs = zeros(numberElements,dofEl); %sizing the element degrees of freedome matrix

K=zeros(sysDof,sysDof); %Sizing the global stiffness matrix


for e=1:numberElements
    elDofs(e,1)=(elementNodes(e,1)*3)-2; %dof U3i-2
    elDofs(e,2)=(elementNodes(e,1)*3)-1; %dof U3i-1
    elDofs(e,3)=(elementNodes(e,1)*3); %dof U3i
    elDofs(e,4)=(elementNodes(e,2)*3)-2; %dof U3j-2
    elDofs(e,5)=(elementNodes(e,2)*3)-1; %dof U3j-1
    elDofs(e,6)=(elementNodes(e,2)*3); %dof U3j
    index1=elementNodes(e,1); %calls node i DOFs
    index2=elementNodes(e,2); %calls node j DOFs
    xx1=xx(index1); %x_i coordinates
    xx2=xx(index2); %x_j coordinates
    yy1=yy(index1); %y_i coordinates
    yy2=yy(index2); %y_j coordinates
    L=((xx2-xx1)^2+(yy2-yy1)^2)^0.5; %Calclating the length of the specific element
    c=(xx2-xx1)/L; %cos(pheta) equals
    s=(yy2-yy1)/L; %sin(pheta) equals
    ALI=A*L^2/I; %A and I constant, length is specific to each element
    EIL=E*I/L^3; %E and I constant for all elements, length is specific to each element.


    Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1]; %Transformation matrix
    kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2]; %element stiffness matrix pre transformation
    ke=transpose(Te)*kebar*Te; %element stiffness post transformation

    K(3*e-2:3*e+3,3*e-2:3*e+3)=K(3*e-2:3*e+3,3*e-2:3*e+3)+ke; %Assembly procedure
end

%Apply boundary conditions
for c = 1:noConstraints
    i = bcDof(c);
    j = bcDof(c);
    K(i,:) = 0; %replace stiffness terms in row i with zero
    K(:,j) = 0; %replace stiffness terms in row j with zero
    P(i) = 0;   %replace force terms in row j with zero
end



activeDof=setdiff([1:sysDof],[bcDof]);  %produce a list of active DOFs
Kr = K(activeDof,activeDof);    %produce a reduced stiffness matrix
Fr = P(activeDof);             %produce a reduced force vector

U = Kr\Fr;  %solve for stiffness

%post-processing to find nodal and element forces

%Resizing the displacement matrix
uFinal=zeros(3*n+3,1)
uFinal(4:3*n,1)=U(1:size(U)-1,1);
uFinal(3*n+3,1)=U(3*n-2,1);


nodalForce = K*uFinal; 
elementForce = zeros(numberElements,dofEl);

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


    Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1];
    kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2];                   %element force
    ke=transpose(Te)*kebar*Te



    gDof = elDofs(e,:);
    cDofStart = gDof(1);
    cDofEnd = gDof(end);
    uCurrent = uFinal(cDofStart:cDofEnd);
    elementForce(e,:) = ke*uCurrent; %forces for the 6 degrees of freedom present in the given element e
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
