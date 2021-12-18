%pass two vector arguments:  xilist, list of temporal frequencies, and fhatlist,
%list of f. coefficients at those frequencies
%** FREQUENCIES MUST BE EVENLY SPACED, from (-N/2)/duration to ((N-1)/2)/duration **

%returns two vectors: tlist, list of times evenly spaced at dt, from [-a to a]
%and flist, list of function values at those times



function [tlist,flist]=inv_f_trans_on_vector(xilist,fhatlist)

    if size(fhatlist,1) > size(fhatlist,2)
        fhatlist = transpose(fhatlist);
    end

    N=length(xilist);
    klist=0:N-1;
    jlist=0:N-1;

    one_over_duration=xilist(2)-xilist(1);
    duration=1/one_over_duration;
    dt=duration/N;

    tlist=-duration/2 + dt*[0:N-1];

    gtildelist=fhatlist.*exp(-pi*1i*jlist);
    ftildelist=ifft(gtildelist);

    flist=1/dt * exp(pi*1i*N/2) * exp(-pi*1i*klist) .*ftildelist;

end

    