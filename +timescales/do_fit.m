function fit = do_fit(xs_ms,diags)
    init = [0.1 100 0.1]; % A, tau (ms), B
    f = @(b,xdata)b(1)*(exp(-xdata/b(2))+b(3));
    b=lsqcurvefit(f,init,xs_ms,diags);
    fit = struct('A',b(1),'tau',b(2),'B',b(3),'f',f,'b',b);
end