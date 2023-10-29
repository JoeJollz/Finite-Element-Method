clear all

%This code is for Question 1, Part C, BOTH i and ii. The user will be asked
%to select between i and ii through an input function.


numberElements=input('Enter the number of elements you would like to process: '); %User is requested to input the number of elements they wish to analysis for
Criteria=input('For both supports being fixed, input 1. For the left support fixed and right support pinned, input 2:' ); %Selecting which boundary conditions they would like to anaylsis for

if Criteria==1  %This following code is for both the supports being fixed.


    n=numberElements; %numberElements stored in simply varible

    %Discretisation process begins
    a=floor((n+3)/6)*2+1; %number of nodes in the left column (including support and top left corner node)
    c=floor(((n+3-a))/2); %number of nodes in the right column (including support and top right corner)
    b=n+3-a-c; %number of nodes in the top beam (including the top left and top right corner nodes

    numberNodes=numberElements+1; %nodes =elements +1
    dofNode=3; %degrees of freedom per node is 3 for a frame
    nodesEl=2; %Nodes per elements is 2; i and j
    dofEl=dofNode*nodesEl; %dof per element is dof per node (3) * number of nodes per element (2)

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
    %Discretisation process ends

    %boundary conditions for the fixed support means the 3 dofs for node 1
    %are present. And the 3 degrees of freedom for the final node is
    %present. As the x and y displacement is restricted along with the
    %rotation.
    bcDof=[numberElements/numberElements;numberElements/numberElements+1;numberElements/numberElements+2;numberNodes*3-2;numberNodes*3-1;numberNodes*3];
    
    A=200; %cm^2 units
    E=20*10^6; %N/cm^2
    I=(6*10^4); %cm^4

    P=zeros(numberNodes*3,1); %sizing the external forces matrix
    P(((a+1)/2)*3-2,1)=150000; %appling the horizontal force of 150kN to the node located at (0,15)

    xx=nodeCoordinates(:,1); %Extracting and storing the x coordinates
    yy=nodeCoordinates(:,2); %Extracting and storing the y coordinates
    sysDof = dofNode*numberNodes; %number of degrees of freedom for the whole system
    noConstraints = length(bcDof); %number of constraints present is equal to the number of boundary conditions
    elDofs = zeros(numberElements,dofEl); %sizing the degree of freedom for each element matrix

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
        L=((xx2-xx1)^2+(yy2-yy1)^2)^0.5; %Calculating length using Pythagoras
        c=(xx2-xx1)/L; %cos(pheta) equal
        s=(yy2-yy1)/L; %sin(pheta) equal
        ALI=A*L^2/I; %A and I are constants for each element, length is specific
        EIL=E*I/L^3; %E and I are constants for each element, length is specific

        Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1]; %Transformation matrix
        kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2]; %element stiffness matrix pre transformation
        ke=transpose(Te)*kebar*Te; %element stiffness matrix 

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

    %Re-sizing the displacement vector
    uFinal=zeros(3*n+3,1)
    uFinal(4:3*n,1)=U;


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
        elementForce(e,:) = ke*uCurrent;
    end

    
    
    
    
    
    %This code is used if the user has selected the left support to be
    %fixed and the right support to the pinned (Allowing rotation,
    %effecting the boundary condtions). 
    %The comments above also apply to this code.
    
else
    n=numberElements;

    a=floor((numberElements+3)/6)*2+1;
    c=floor(((numberElements+3-a))/2);
    b=numberElements+3-a-c;

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


    %Now changed, rotation is free at the final node due to support now
    %being pinned. So the very final global degree of freedom is not in the
    %bcDof matrix
    bcDof=[numberElements/numberElements;numberElements/numberElements+1;numberElements/numberElements+2;numberNodes*3-2;numberNodes*3-1];

    A=200; %Defining the cross-sectional areas for each element
    E=20*10^6; %Defining the youngs modulus of each element
    I=(6*10^4);

    P=zeros(numberNodes*3,1);
    P(((a+1)/2)*3-2,1)=150000;

    xx=nodeCoordinates(:,1);
    yy=nodeCoordinates(:,2);
    sysDof = dofNode*numberNodes;	        
    noConstraints = length(bcDof);          
    elDofs = zeros(numberElements,dofEl);

    L=zeros(1,numberElements);
    C=zeros(1,numberElements);
    S=zeros(1,numberElements);


    K=zeros(sysDof,sysDof);


    for e=1:numberElements
        elDofs(e,1)=(elementNodes(e,1)*3)-2;
        elDofs(e,2)=(elementNodes(e,1)*3)-1;
        elDofs(e,3)=(elementNodes(e,1)*3);
        elDofs(e,4)=(elementNodes(e,2)*3)-2;
        elDofs(e,5)=(elementNodes(e,2)*3)-1;
        elDofs(e,6)=(elementNodes(e,2)*3);
        index1=elementNodes(e,1);
        index2=elementNodes(e,2);
        xx1=xx(index1);
        xx2=xx(index2);
        yy1=yy(index1);
        yy2=yy(index2);
        L=((xx2-xx1)^2+(yy2-yy1)^2)^0.5;
        c=(xx2-xx1)/L;
        s=(yy2-yy1)/L;
        ALI=A*L^2/I;
        EIL=E*I/L^3;


        Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1];
        kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2];
        ke=transpose(Te)*kebar*Te;

        K(3*e-2:3*e+3,3*e-2:3*e+3)=K(3*e-2:3*e+3,3*e-2:3*e+3)+ke;
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
        s = S(e);                                   %renamed for convenience                                  %definition
        ALI=A*L^2/I;
        EIL=E*I/L^3;                               %definition


        Te=[c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; 0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1];
        kebar=EIL*[ALI 0 0 -ALI 0 0; 0 12 6*L 0 -12 6*L; 0 6*L 4*L^2 0 -6*L 2*L^2; -ALI 0 0 ALI 0 0; 0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2];                   %element force
        ke=transpose(Te)*kebar*Te



        gDof = elDofs(e,:);
        cDofStart = gDof(1);
        cDofEnd = gDof(end);
        uCurrent = uFinal(cDofStart:cDofEnd);
        elementForce(e,:) = ke*uCurrent;
    end
    
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
