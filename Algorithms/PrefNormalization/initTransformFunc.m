function TransformFunc = initTransformFunc(Pref)
    M = length(Pref);
    TransformFunc = cell(1, M);
    for m = 1:M
        pref = Pref{m}; 
        if length(pref) == 2
            %% beta distribution type
            a = pref(1);
            b = pref(2);
            cdf_func = @(x) betacdf(x,a,b);
        else
            %% collected data type
            data = pref;
            data = (data-min(data))./(max(data)-min(data));
            % calculate ecdf
            [cdff,cdfx] = ecdf(data);
            % calculate breakpoints
            xj = cdfx(2:end);
            Fj = (cdff(1:end-1)+cdff(2:end))/2;
            % expand breakpoints
            n = length(xj);
            xj_expanded = [xj(1)-Fj(1)*(xj(2)-xj(1))/((Fj(2)-Fj(1)));
                  xj;
                  xj(n)+(1-Fj(n))*((xj(n)-xj(n-1))/(Fj(n)-Fj(n-1)))];
            Fj_expanded = [0; Fj; 1];
            % linearly interpolate
            cdf_func = @(x) interp1(xj_expanded,Fj_expanded,x,'linear','extrap');
        end
        TransformFunc{m} = cdf_func;
    end
end