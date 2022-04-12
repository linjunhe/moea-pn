function SPEA2PN2_U(Global)
% <algorithm> <S>
% Strength Pareto evolutionary algorithm 2

%------------------------------- Reference --------------------------------
% E. Zitzler, M. Laumanns, and L. Thiele, SPEA2: Improving the strength
% Pareto evolutionary algorithm, Proceedings of the Fifth Conference on
% Evolutionary Methods for Design, Optimization and Control with
% Applications to Industrial Problems, 2001, 95-100.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
    PopObj  = Population.objs;
    Zmin    = min(Population.objs,[],1);

    % articulate preference
    Pref = {[1, 1],...
            [1, 1]};
    
    % initialize transformation functions and transform
    TransformFunc = initTransformFunc(Pref);
    PopObj = nonlinearTransform(PopObj, TransformFunc, Zmin, Global.gen/Global.maxgen);

    Fitness    = CalFitness(PopObj);
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,Fitness);
        Offspring  = GA(Population(MatingPool));
        Zmin       = min([Zmin;Offspring.objs],[],1);
        [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Global.N,TransformFunc,Zmin,Global.gen/Global.maxgen);
    end
end

function [Population,Fitness] = EnvironmentalSelection(Population,N,TransformFunc,Zmin,t)
% The environmental selection of SPEA2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj  = Population.objs;
    
    % nonlinear transform
    PopObj = nonlinearTransform(PopObj, TransformFunc,Zmin,t);
    
    %% Calculate the fitness of each solution
    Fitness = CalFitness(PopObj);

    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(PopObj(Next,:),sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end

function Fitness = CalFitness(PopObj)
% Calculate the fitness of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);

    %% Detect the dominance relation between each two solutions
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end
    
    %% Calculate S(i)
    S = sum(Dominate,2);
    
    %% Calculate R(i)
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end
    
    %% Calculate D(i)
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);
    
    %% Calculate the fitnesses
    Fitness = R + D';
end