clear
close all
tic
addpath(genpath('~/Documents/MATLAB/Images/altmany-export_fig-c062ae2'));

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% User Parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
FileBase = '/Users/fhilty/Desktop/Design-Tool-3D/Cont/G%d-%d.png'; %When doing batch runs this is the file base and destination, %d is replaced by IC's value
SaveBase = '/Users/fhilty/Desktop/Design-Tool-3D/Cont/G%d-%d'; %Base file name/location to save select workspace information when complete
CSVBase  = '/Users/fhilty/Desktop/Design-Tool-3D/Cont/G10.csv'; %Base file name/location to save information from 'Data' array when complete
Threshold = 100; %Threashhold value to differentiate phases in image
SegHeight = 3; %Height of a segment of image, determines the available space to search to find the best pathway as Image height divided by SegHeight
SegShift = 3; %Distance segments shift down to find the next search space for the next pathway as a segments height divided by SegShift

PerCut = 0.03; %Fraction of child pathways which are kept, used as cutoff value for ending the Genetic Algorithm part when less than PerCut of child pathways are kept
NumPath = 150; %Number of parent pathways making up the gene pool per segment
Mutation = 0.03;  %Mutation chance each primary node has when generating a child path in Genetic part
PD = 3; %Determine distance between primary nodes, primary nodes receive a random y-coordinate and are passed between parent and child
DPA = 5; %Displacements Per Attempt, the number of primary nodes randomly displaced during each attempt at finding a better path in Annealing Algorithm part
NumRun = 5; %Number of times to run algorithm, best result from all runs will be kept
SmCut = 0.95; %Cutoff value to determine if a smoothing operation should be kept
CPMax = 5000000; %Hard limit on how many attempted pathways to make

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% User Parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine segment boundaries
Eval = zeros(NumPath,3); % C=1 D=2 P=3
CEval = zeros(1,3); % C=1 D=2 P=3
Data = zeros(1,7); %FN,SR,IC,C,D,P,fa   This is saved at the end as a record
d=1;
for FN = [1] %if doing a batch of images this can be used as one variable in the file name
    for SR = [10] %Likewise this can be used as one variable in the file name
        for IC = [25,50,75,100] %This is the third possible variable in a file name for batch runs
            I = imread(sprintf(FileBase,SR,IC)); %Read in the desired image
            I = rgb2gray(I); %Convert image to grayscale
            [Iy,Ix] = size(I);

            %Find the volume percent of the microstructure with set Threshold
            s = 0;
            for m = 1:Iy
                for n = 1:Ix
                    if I(m,n)<Threshold
                        s = s + 1;
                    end
                end
            end
            VolPercent = s/(Ix*Iy);

            AT = 100*Ix; %Annealing Time, Number of times changes to a NewPath are attempted
            Segments = SegHeight*SegShift; %Determine number of segments
            NewY = round(Iy + (SegShift-1)*Iy/Segments); 
            Bounds(1,Segments+SegShift) = NewY; %Y coordinate for boundaries defining the search space for a pathway
            Bounds(1,1) = 1;
            for m = 1:Segments+SegShift-2
               Bounds(1,m+1) =  floor(m*Iy/Segments);
            end

            %Stitch top portion of image to bottom
            I(NewY,Ix) = 0; %Resize I to allow additional segments to be copied from top to bottom
            for m = 1:round((SegShift-1)*Iy/Segments)
                I(Iy+m,:) = I(m,:);
            end

            %Clean Image based on Threshold
            for m = 1:NewY
                for n = 1:Ix
                    if (I(m,n) >= Threshold) %Change > sign to < if secondary phase isnt dark phase
                        I(m,n) = 255; %White
                    else
                        I(m,n) = 0; %Black
                    end
                end
            end

            %Determine Primary Nodes
            if (rem(Ix,PD)==0 || rem(Ix,PD)==1)
                NPN = floor(Ix/PD);
            else
                NPN = floor(Ix/PD)+1;
            end
            PN = zeros(NPN+1,1); %Stores X-coordinate of each primary node
            for m = 0:NPN
                PN(m+1) = 1+m*PD;
            end
            if (PN(NPN+1) > Ix)
                PN(NPN+1) = Ix;
            end

            %The search routine to find a optimum pathway NumRun times and keeps best for each segment        
            FinPath = zeros(Ix,Segments); %The final pathways from the search
            FinEval = zeros(Segments,3); %final pathway evaluations, the fraction of pathway in secondary phase, normalized pathway lengths, and continuousness value
            for R = 1:NumRun %Loop for the full searching algorithm

                %Genetic Algorithm to find best pathway
                BestPath = zeros(Ix,Segments); %The best pathway from the gene pool from this run of the search
                BestEval = zeros(Segments,3); %Evaluations for the best pathway
                for S = 1:Segments

                    %Make first guess for all primary nodes
                    Nodes = zeros(Ix,NumPath); %Stores all Y-coordinates making up all pathways, Nodes' index is the X-coordinate 
                    for m = 1:NumPath
                        for n = 1:NPN+1
                             Nodes(PN(n),m) = randi([Bounds(1,S),Bounds(1,S+SegShift)]); %Assign random Y-coordinate within segments boundaries to each primary node
                        end
                    end

                    %Place guide nodes, Guide nodes are placed in a straight line between primary nodes
                    for m = 1:NumPath
                        for n = 1:NPN
                            for i = 1:PN(n+1)-PN(n)-1 %1 to the number of guide nodes between primary nodes
                                Nodes((n-1)*PD+i+1,m) = round(Nodes(PN(n),m)+i*((Nodes(PN(n+1),m)-Nodes(PN(n),m))/(PN(n+1)-PN(n)))); %Calculate position of each guide node
                            end
                        end
                    end

                    %Smooth lines between phase boundaries
                    Nodes = smoothPaths(NumPath,Threshold,Ix,Nodes,I);

                    %Evaluate Initial Guesses
                    Eval = evalPaths(NumPath,Ix,Nodes,I,Eval);

                    %Evolve system of paths with Genetic Algorithm until the percent cuttoff value is reached
                    PerAcc = 1; %Initialize PerAcc to enter while loop
                    CPAttempts = 1; %The number of child pathways attempted
                    CPAccepted = 0; %The number of child pathways accepted into gene pool
                    CPTotal = 1;
                    while(PerAcc > PerCut && CPTotal < CPMax)

                        %Mix two pathways
                        PA = 0; %Parent pathway A
                        PB = 0; %Parent pathway B
                        while PA == PB %Pick two different random pathways
                           PA = randi([1,NumPath]);
                           PB = randi([1,NumPath]);
                        end

                        %Mix primary nodes between the two parents
                        Child = zeros(Ix,1); %Stores child pathway node coordinates
                        for m = 1:NPN+1
                            if rand < Mutation %Mutation chance
                                Child(PN(m)) = randi([Bounds(1,S),Bounds(1,S+SegShift)]);
                            elseif rand > rand %Parent pathway A passes on gene
                                Child(PN(m)) = Nodes(PN(m),PA);
                            else %Parent pathway B passes on gene
                                Child(PN(m)) = Nodes(PN(m),PB);
                            end
                        end

                        %Place guide nodes for Child
                        for m = 1:NPN
                            for n = 1:PN(m+1)-PN(m)-1
                                Child((m-1)*PD+n+1) = round(Child(PN(m))+n*((Child(PN(m+1))-Child(PN(m)))/(PN(m+1)-PN(m))));
                            end
                        end

                        %Smooth line between phase boundaries
                        Child = smoothPaths(1,Threshold,Ix,Child,I);

                        %Evaluate child pathway
                        CEval = evalPaths(1,Ix,Child,I,CEval);

                        %Replace worst parent pathway with child pathway if it is better
                        if Eval(PB,3) >= Eval(PA,3) %Determine worst parent pathway
                            low = PA;
                        else
                            low = PB;
                        end
                        if CEval(1,1)/(CEval(1,2)) > Eval(low,1)/(Eval(low,2)) %If child pathway is better it replaces worst parent
                            Nodes(:,low) = Child(:);
                            Eval(low,:) = CEval(:);
                            CPAccepted = CPAccepted +1; %Incriment number of child pathways accepted into gene pool
                        end

                        %Determine the delta mean P periodically
                        if CPAccepted == 300 || CPAttempts > 100000 %When we have tried 100000 new child pathways or 300 have been accepted
                            PerAcc = CPAccepted/CPAttempts; %Update percent accepted to control while loop
                            CPAttempts = 1; %The number of child pathways attempted
                            CPAccepted = 0; %The number of child pathways accepted into gene pool

                            %Clean all pathways by eliminating zig-zag pattern from grains
                            Nodes = zigZag(NumPath,Threshold,SmCut,Ix,Nodes,I);

                            %Re-evaluate FinEval after removing zig-zag pattern
                            Eval = evalPaths(NumPath,Ix,Nodes,I,Eval);
                        end
                        CPAttempts = CPAttempts + 1;
                        CPTotal = CPTotal + 1;
                    end

                    %Smooth lines between phase boundaries
                    Nodes = smoothPaths(NumPath,Threshold,Ix,Nodes,I);

                    %Evaluate All Paths
                    Eval = evalPaths(NumPath,Ix,Nodes,I,Eval);

                    %Find best path for current segment
                    best = 1;
                    for m = 2:NumPath
                        if Eval(m,1)/(Eval(m,2)) > Eval(best,1)/(Eval(best,2))
                            best = m;
                        end
                    end
                    BestPath(:,S) = Nodes(:,best);
                    BestEval(S,:) = Eval(best,:);
                end

                %Annealing Algorithm to refine best pathways
                for S = 1:Segments

                    %Annealing loop
                    for T = 1:AT
                    %Neighbor creation
                        %Shift points slightly along best path to make Child path
                        sigma = 0.05*(Bounds(1,S+SegShift) - Bounds(1,S)); %determines displacement distribution
                        Child(:,1) = BestPath(:,S);
                        for m = 1:DPA
                            point = randi([1,Ix]); %Randomly pick a point from the line by its x coordinate
                            NewPoint = round(normrnd(BestPath(point,S),sigma)); %Displace that point
                            if NewPoint < Bounds(1,S) %Keep point within bounds
                                NewPoint = Bounds(1,S);
                            elseif NewPoint > Bounds(1,S+SegShift)
                                NewPoint = Bounds(1,S+SegShift);
                            end
                            Child(point,1) = NewPoint; %Make change to pathway
                        end

                        %Smooth Child between phase boundaries
                        Child = smoothPaths(1,Threshold,Ix,Child,I);

                        %Evaluate Child
                        CEval = evalPaths(1,Ix,Child,I,CEval);

                        %If Child has a higher continuousness keep it
                        if CEval(1,3) > BestEval(S,3)
                            BestPath(:,S) = Child(:,1);
                            BestEval(S,:) = CEval(1,:);
                        end
                    end
                     if  R == 1 || BestEval(S,1)/(BestEval(S,2)) > FinEval(S,1)/(FinEval(S,2)) %If the best pathway of this run is better than previous best then make it the Final Pathway
                        FinPath(:,S) = BestPath(:,S);
                        FinEval(S,:) = BestEval(S,:);            
                     end
                end
            end

            for R = 1:2

                %Eliminate zig-zag pattern in grains
                FinPath = zigZag(Segments,Threshold,SmCut,Ix,FinPath,I);

                %Smooth Final Pathways
                FinPath = smoothPaths(Segments,Threshold,Ix,FinPath,I);

                %Remove spikes
                for m = 1:Segments
                   for n = 1:Ix-1
                       if abs(FinPath(n,m) - FinPath(n+1,m)) > 5
                          if n < 6 %Near left edge   
                            FinPath(1:n-1,m) = FinPath(n,m);
                            Step = (FinPath(n+5,m) - FinPath(n,m))/5;
                            for i = n+1:n+4
                                FinPath(i,m) = round(FinPath(n,m)+((i-n)*Step));
                            end                
                          elseif n > Ix - 6 %Near right edge
                            FinPath(n+1:Ix,m) = FinPath(n,m);
                            Step = (FinPath(n,m) - FinPath(n-5,m))/5;
                            for i = n-4:n-1
                                FinPath(i,m) = round(FinPath(n,m)+((i-n+5)*Step));
                            end  
                          else %In interior
                            Step = (FinPath(n+5,m) - FinPath(n-5,m))/10;
                            for i = n-4:n+4
                                FinPath(i,m) = round(FinPath(n-5,m)+((i-n+5)*Step));
                            end  
                          end
                       end
                   end
                end
            end

            %Evaluate new FinPath
            FinEval = evalPaths(Segments,Ix,FinPath,I,FinEval);

            %Populate Data array with results
            Data(d,1) = FN;
            Data(d,2) = SR;
            Data(d,3) = IC;
            Data(d,4) = mean(FinEval(:,1));
            Data(d,5) = mean(FinEval(:,2));
            Data(d,6) = mean(FinEval(:,3));
            Data(d,7) = VolPercent;
            d=d+1;
            
            %Display final results
            disp(fprintf(FileBase,SR,IC))
            disp('Percent Dark Phase');
            disp(VolPercent);
            disp('Ave    C         D         P');
            disp(mean(FinEval));
            disp(' ')

            %Plot final pathways
            imshow(I)
            hold on
            for m = 1:Segments
                plot(FinPath(:,m),'linewidth',4)       
            end

            %Save control parameters and final results from matlab workspade and save figure
    %         save(sprintf(SaveBase,SR,IC),'Threshold','SegHeight','SegShift','PerCut','NumPath','Mutation','PD','DPA','NumRun','SmCut','FinPath','FinEval')
    %         export_fig(sprintf(FileBase, IC), '-jpg');

            toc
        end
    end
end
csvwrite(CSVBase,Data)

%Smooth lines between phase boundaries
function Nodes = smoothPaths(NumPath,Threshold,Ix,Nodes,I)
    for m = 1:NumPath
        PhBnd = [1 Nodes(1,m)]; %Set first point as first phase boundary
        Phase = logical(I(Nodes(1,m),1) <= Threshold); %Determine phase of first point
        for n = 2:Ix %Loop through all x points of a path
            if ((n <= Ix-2) && (Phase ~= logical(I(Nodes(n,m),n) <= Threshold)) && (Phase ~= logical(I(Nodes(n+1,m),n+1)) <= Threshold) && (Phase ~= logical(I(Nodes(n+2,m),n+2) <= Threshold))) %Smooth between this phase boundary and previous one
                Step = (Nodes(n,m) - PhBnd(1,2))/(n - PhBnd(1,1));
                for i = PhBnd(1,1)+1:n-1
                    Nodes(i,m) = round(Nodes(PhBnd(1,1),m)+((i-PhBnd(1,1))*Step));
                end
                Phase = logical(I(Nodes(n,m),n) <= Threshold);
                PhBnd = [n Nodes(n,m)]; 
            end
            if (n == Ix) %When at the right border of image smooth to previous phase boundary
                Step = (Nodes(n,m) - PhBnd(1,2))/(n - PhBnd(1,1));
                for i = PhBnd(1,1)+1:n-1
                    Nodes(i,m) = round(Nodes(PhBnd(1,1),m)+((i-PhBnd(1,1))*Step));
                end
                Phase = logical(I(Nodes(n,m),n) <= Threshold);
                PhBnd = [n Nodes(n,m)];
            end
        end
    end
end

%Evaluate Paths
function Eval = evalPaths(NumPath,Ix,Nodes,I,Eval)
    for m = 1:NumPath
        PLS = 0; %Path Length in Secondary Phase
        TPL = 0; %Total Path Length
        for n = 1:Ix-1
            dl = sqrt(1+(Nodes(n+1,m)-Nodes(n,m))^2); %Length of path segment between nodes
            if I(Nodes(n,m),n) == 0 && I(Nodes(n+1,m),n+1) == 0 %Both points in secondary phase
                PLS = PLS + dl;
            elseif I(Nodes(n,m),n) ~= I(Nodes(n+1,m),n+1) %Points in different phases
                PLS = PLS + 0.5*dl;
            end 
            TPL = TPL + dl;
        end
        Eval(m,1) = PLS/TPL; %Percentage of pathway in secondary phase
        Eval(m,2) = TPL/Ix;  %Normalized pathway legnth
        Eval(m,3) = Eval(m,1)/Eval(m,2); %Cointinousness of pathway
    end
end

%Eliminate Zig-Zag
function Nodes = zigZag(NumPath,Threshold,SmCut,Ix,Nodes,I)
    for m = 1:NumPath
        Child(:,1) = Nodes(:,m);
        PhBnd = zeros(6,2); %X and Y coordinates for 6 phase boundary points
        PhBnd(1,1) = 1;
        PhBnd(1,2) = Child(1,1); %Set first point as first phase boundary
        Phase = logical(I(Child(1,1),1) <= Threshold); %Determine phase of first point
        for n = 2:Ix %Loop through all x points of a path
            if Phase ~= logical(I(Child(n,1),n) <= Threshold) %Find next phase boundary
                b = n; %Incriment to track x-coordinate of current point
                PB = 1; %Phase boundaries kept so far
                while PB < 6 && b < Ix %Find up to the next 5 phase boundaries
                    if (Phase ~= logical(I(Child(b,1),b) <= Threshold)) 
                        PB = PB + 1;
                        PhBnd(PB,1) = b;
                        PhBnd(PB,2) = Child(b,1);
                        Phase = logical(I(Child(b,1),b) <= Threshold); %Update phase
                    end
                    b = b + 1;
                end
                j = 0; %Incriments which phase boundary is being considered
                while j < PB %Starts with first and last phase boundaries and moves on to closer phase boundaries to first 
                    dy = PhBnd(PB-j,2) - PhBnd(1,2); %vertical distance between first and current phase boundary
                    dx = PhBnd(PB-j,1) - PhBnd(1,1); %horizontal distance between first and current phase boundary
                    for i = PhBnd(1,1)+1:PhBnd(PB-j,1)-1 %Create straight line between first and current phase boundary
                        Child(i,1) = round(PhBnd(1,2)+((i-PhBnd(1,1))*dy/dx));
                    end
                    PLS = 0; %Path Length in Secondary Phase
                    TPL = 0; %Total Path Length
                    for i = 0:dx-1 %incriment through x-coordinates
                        cx = PhBnd(1,1)+i; %Current X point
                        cy = round(PhBnd(1,2)+(i*dy/dx)); %Current Y point
                        nx = cx+1; %Next X point
                        ny = round(PhBnd(1,2)+((i+1)*dy/dx)); %Next Y point
                        dl = sqrt(1+(ny-cy)^2);
                        if I(cy,cx) == 0 && I(ny,nx) == 0 %Both points in secondary phase
                            PLS = PLS + dl;
                        elseif I(cy,cx) ~= I(ny,nx) %Points in different phases
                            PLS = PLS + 0.5*dl;
                        end
                        TPL = TPL + dl;
                    end
                    j = j + 1; %Update j
                    if PLS/TPL > SmCut %If new path segment is in secondary phase more than the smoothing cutoff value keep it
                        Nodes(:,m) = Child(:,1);
                        PhBnd(1,1) = PhBnd(PB-j+1,1); %Make current phase boundary the new starting phase boundary
                        PhBnd(1,2) = PhBnd(PB-j+1,2);
                        Phase = logical(I(PhBnd(1,2),PhBnd(1,1)) <= Threshold);
                        j = PB; %When a smoothing operation is successful stop while loop
                    end
                    Child(:,1) = Nodes(:,m); %Removes the changes to the pathway if a smoothing operation was unsuccessful
                end
            end
        end
    end
end