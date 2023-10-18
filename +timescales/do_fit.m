function fit = do_fit(xs,autocorr,params)
    arguments
        xs
        autocorr
        params.skip_refrac (1,1) logical  = true;
        params.init (1,3) = [0.1 0.1 1]; % good for autocorrelation timescale when data is in seconds
    end
    if params.skip_refrac
        [~,idx] = min(diff(autocorr));
        if idx>1
            warning('Removed the 1st %d datapoint(s) from fitting due to a refractory period.',idx-1);
        end
    else
        idx=1;
    end
    f = @(b,xdata)b(1)*(exp(-xdata/b(2))+b(3));
    b=lsqcurvefit(f,params.init,xs(idx:end),autocorr(idx:end));
    fit = struct('A',b(1),'tau',b(2),'B',b(3),'f',f,'b',b);
end