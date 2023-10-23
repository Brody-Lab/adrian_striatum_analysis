function fit = do_fit(xs,autocorr,params)
    arguments
        xs
        autocorr
        params.max_skip_refrac (1,1) {mustBeInteger,mustBeNonnegative} = 0;
        params.init (1,3) = [0.1 0.1 1]; % good for autocorrelation timescale when data is in seconds
        params.verbose (1,1) logical = false;
    end
    if params.max_skip_refrac>0
        [~,idx] = min(diff(autocorr));%find(diff(autocorr)<0,1); 
        idx=min(idx,1+params.max_skip_refrac);        
        if idx>1 && params.verbose
            warning('Removed the 1st %d datapoint(s) from fitting due to a refractory period.',idx-1);
        end
    else
        idx=1;
    end
    f = @(b,xdata)b(1)*(exp(-xdata/b(2))+b(3));
    options = optimoptions('lsqcurvefit','Display','off');    
    b=lsqcurvefit(f,params.init,xs(idx:end),autocorr(idx:end),[],[],options);
    fit = struct('A',b(1),'tau',b(2),'B',b(3),'f',f,'b',b);
end