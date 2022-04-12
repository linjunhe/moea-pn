function PopObj = nonlinearTransform(varargin)

    PopObj = varargin{1};
    [N,M]  = size(PopObj);
    t = varargin{4};

    if size(varargin{3},1) == 1
        % normalize population with current ND solutions
        TransformFunc = varargin{2};
        Zmin = varargin{3};
        FrontNo = NDSort(PopObj,N+1);
        Zmax = max(PopObj(FrontNo==1,:),[],1);
        PopObj1 = (PopObj - Zmin)./(Zmax - Zmin + 1e-6);
    else
        % normalize population with PF
        TransformFunc = varargin{2};
        PF = varargin{3};
        fmin   = min(PF,[],1);
        fmax   = max(PF,[],1);
        PopObj1 = (PopObj-fmin)./(fmax-fmin);
    end

    PopObj2 = PopObj1;
    % transform
    for m = 1:M
        x = PopObj1(:,m);
        y = TransformFunc{m}(x);
        % linear mapping for solutions outside [0,1]^M
        y(x>1|x<0) = x(x>1|x<0);
        PopObj2(:,m) = y;
    end

    PopObj = (1-t)*PopObj1 + t*PopObj2;

end