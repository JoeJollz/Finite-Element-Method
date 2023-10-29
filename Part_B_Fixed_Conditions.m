clear all

%This code is for Question 1, Part B, Part i. Which is 40 elements with both
%supports being fixed

numberElements=40; %Number of elements to analyse for
n=numberElements; %Small notation, was oringally labelled numberElements to make clearer to reader

%Discretisation process begins 

a=floor((numberElements+3)/6)*2+1; %The number of nodes in the left column (include support and top node)
c=floor(((numberElements+3-a))/2); %Number of nodes in the right column (including support and top node)
b=numberElements+3-a-c; %number of nodes in the top beam (including the top left and right corners)

numberNodes=numberElements+1;%automatic
dofNode=3;
nodesEl=2;
dofEl=dofNode*nodesEl;

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
%Discretisation process ends. Open the nodeCoordinates data to see that
%The nodes are OPTIMALLY DISTRIBUTED 

%The boundary conditions are corresponding the to the degrees of freedom on
%the first and last node. Both supports are fixed, hence the 3
%degress of freedom on the first node are present and the 3 DoF on the last
%node are present
bcDof=[numberElements/numberElements;numberElements/numberElements+1;numberElements/numberElements+2;numberNodes*3-2;numberNodes*3-1;numberNodes*3];

A=200; %cm^2
E=20*10^6; %N/cm^2
I=(6*10^4); %cm^4

P=zeros(numberNodes*3,1); %External force application sizing
P(((a+1)/2)*3-2,1)=150000; %Horizontal point force applied in the middle of left column (newtons)

xx=nodeCoordinates(:,1); %extracting and storing the x coordinates
yy=nodeCoordinates(:,2); %extracting and storing the y coordinates
sysDof = dofNode*numberNodes; %Total number of Dofs in the system
noConstraints = length(bcDof); %Used in a later loop, equal to the number of boundary conditions present       
elDofs = zeros(numberElements,dofEl); %Sizing the degrees of freedom for an element

K=zeros(sysDof,sysDof); %Sizing the global stiffness matrix for the system.


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
    L=((xx2-xx1)^2+(yy2-yy1)^2)^0.5; %Element length calculation from pythagoras 
    c=(xx2-xx1)/L; %cos(pheta) equals
    s=(yy2-yy1)/L; %sin(pheta) equals
    ALI=A*L^2/I; %A and I are constants for each element, length is specific
    EIL=E*I/L^3; %E and I are constants for each element, length is specific


    Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1]; %Transformation matrix
    kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2]; %element stiffness pre transformation
    ke=transpose(Te)*kebar*Te; %element stiffness post ransformation

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

U = Kr\Fr;  %solve for displacements

%Post Processing to find element forces and nodal forces

uFinal=zeros(3*n+3,1) %Resizing the displacement matrix
uFinal(4:3*n,1)=U; 


nodalForce = K*uFinal; %Forces on nodes
elementForce = zeros(numberElements,dofEl); %sizing the element forces matrix for 6 degree of freedom per element

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
    elementForce(e,:) = ke*uCurrent; %element force specific to element e
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


    