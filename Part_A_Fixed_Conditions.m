clear all

%This code is for Question 1, Part A, Part i. Which is 4 elements with both
%supports being fixed.

numberElements=4;    %Number of elements to analyse for
numberNodes=numberElements+1;   %Nodes=elements+1
dofNode=3;                %3 degrees of freedom present, x, y, rotation
nodesEl=2;                %node i and node j
dofEl=dofNode*nodesEl;    %The degrees of freedom for a single element is equal to degrees of freedom per node * nodes per element

%Discretisation process begins

elementNodes=[1 2; 2 3; 3 4; 4 5]; %The global node labelling system. First elements (1 2), second element (2 3), etc
nodeCoordinates=[0 0; 0 15; 0 30; 60 30; 60 15]; %Coordinates for plotting the frame
bcDof=[1,2,3,13,14,15]; %Boundary conditions of the frame, these are found at node 1 and node 5. 

%Discretisation process complete

A=200; %units cm^2
E=20*10^6   %units N/cm^2
I=10^4*6;   %cm^4

P=[0;0;0;150000;0;0;0;0;0;0;0;0;0;0;0]; %force applied at node 2, x direction, hence degree of freedom 4.

xx=nodeCoordinates(:,1); %Storing the x coordinate data
yy=nodeCoordinates(:,2); %Storing the y coordinate data
sysDof=dofNode*numberNodes; %Total number of degrees of freedom for the system.
noConstraints=length(bcDof); %No of constraints will be equal to the number of boundary conditions.
elDofs=zeros(numberElements,dofEl); %Sizing a matrix ready to store all the degrees of freedom.

K=zeros(sysDof,sysDof); %Sizing the global stiffness matrix


for e=1:numberElements
    elDofs(e,1)=(elementNodes(e,1)*3)-2; %Calculating for the given element e; the degree of freedom U(3i-2)
    elDofs(e,2)=(elementNodes(e,1)*3)-1; %Calculating for the given element e; the degree of freedom U(3i-1)
    elDofs(e,3)=(elementNodes(e,1)*3); %Calculating for the given element e; the degree of freedom U(3i)
    elDofs(e,4)=(elementNodes(e,2)*3)-2; %Calculating for the given element e; the degree of freedom U(3j-2)
    elDofs(e,5)=(elementNodes(e,2)*3)-1; %Calculating for the given element e; the degree of freedom U(3j-1)
    elDofs(e,6)=(elementNodes(e,2)*3); %Calculating for the given element e; the degree of freedom U(3j)
    index1=elementNodes(e,1); %Retriving the index, which will be used to extract x and y coordination for the given element e.
    index2=elementNodes(e,2); %Retriving the index, which will be used to extract x and y coordination for the given element e.
    xx1=xx(index1); %Extracting the x coordinate relating to the i node of the element
    xx2=xx(index2); %Extracting the x coordinate relating to the j node of the element 
    yy1=yy(index1); %Extracting the y coordinate relating to the i node of the element
    yy2=yy(index2); %%Extracting the x coordinate relating to the j node of the element
    L=((xx2-xx1)^2+(yy2-yy1)^2)^0.5; %Pythagoras Theorem
    c=(xx2-xx1)/L; %cos(pheta) is equal to
    s=(yy2-yy1)/L; %sin(pheta) is equal to

    ALI=A*L^2/I; %The A and I are constant for all elements, The length is specific
    EIL=E*I/L^3; %The E and I are constant for all elements, The length is specific
    
    
    Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1]; %Transformation matrix
    kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2]; %kebar
    ke=transpose(Te)*kebar*Te; %transformed ke
    
    K(3*e-2:3*e+3,3*e-2:3*e+3)=K(3*e-2:3*e+3,3*e-2:3*e+3)+ke; %The assembly procedure
end

%Applying boundary conditions to the global stiffness matrix and and the
%externally applied load
noConstraints=6;
for c=1:noConstraints
    i=bcDof(c);
    j=bcDof(c);
    K(i,:)=0;
    K(:,j)=0;
    P(i)=0;
end


%Creating the reduced stiffness matrix and reduced force matrix
RefinedK=K(4:12,4:12);
RefinedP=P(4:12,1);

%Calculating the displacement and roation
U=RefinedK\RefinedP;

%Creating the appropriate size displacement and rotation matrix
uFinal=zeros(3*4+3,1)
uFinal(4:3*4,1)=U;

%POST_PROCESSING

nodalForce = K*uFinal; %Forces and moments in the nodes
elementForce = zeros(numberElements,dofEl); %sizing an matrix to fit 6 degrees of freedoms for all the elements

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


    Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1]; %Transformation matrix 
    kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2];                   %element force
    ke=transpose(Te)*kebar*Te %element stiffness matrix after being transformed



    gDof = elDofs(e,:); 
    cDofStart = gDof(1);
    cDofEnd = gDof(end);
    uCurrent = uFinal(cDofStart:cDofEnd);
    elementForce(e,:) = ke*uCurrent; %Storing the forces and moments for each i and j of each element
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

    