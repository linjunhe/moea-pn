function SMSEMOAPN_1(Global)
% <algorithm> <SMSEMOA>
% S metric selection based evolutionary multiobjective optimization with
% Pref1
% algorithm
% r --- 1.1 --- r of reference point

%------------------------------- Reference --------------------------------
% M. Emmerich, N. Beume, and B. Naujoks, An EMO algorithm using the
% hypervolume measure as selection criterion, Proceedings of the
% International Conference on Evolutionary Multi-Criterion Optimization,
% 2005, 62-76.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    r = Global.ParameterSet(1.1);
    
    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NDSort(Population.objs,inf);

    % articulate preference
    Pref = {[1, 10],...
            [1, 1]};
    
    Zmin = min(Population.objs);
    % initialize transformation functions and transform
    TransformFunc = initTransformFunc(Pref);
    
    %% Optimization
    while Global.NotTermination(Population)
        for i = 1 : Global.N
            drawnow();
            Offspring = GAhalf(Population(randperm(end,2)));
            Zmin = min(Zmin, Offspring.obj);
            [Population,FrontNo] = Reduce([Population,Offspring],FrontNo,r,TransformFunc,Zmin,Global.gen/Global.maxgen);
        end
    end
end

function [Population,FrontNo] = Reduce(Population,FrontNo,r,TransformFunc,Zmin,t)
% Delete one solution from the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Identify the solutions in the last front
    FrontNo   = UpdateFront(Population.objs,FrontNo);
    LastFront = find(FrontNo==max(FrontNo));
    PopObj    = Population(LastFront).objs;
    [N,M]     = size(PopObj);
    
    % nonlinear transform
    PopObj = nonlinearTransform(PopObj,TransformFunc,Zmin,t);

    %% Calculate the contribution of hypervolume of each solution
    deltaS = inf(1,N);
    ref = r*ones(1,M);

    if N > 1
        deltaS = CalHVC(PopObj,ref,N);
    end
    
    %% Delete the worst solution from the last front
    [~,worst] = min(deltaS);
    FrontNo   = UpdateFront(Population.objs,FrontNo,LastFront(worst));
    Population(LastFront(worst)) = [];
end