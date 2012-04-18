function bout=B(t,s,h,dr)
% 
% Dynamics of the ethanolamine glycerophospholipid lipid remodeling network
%
% Copyright Mar 2012, Lu Zhang, all rights reserved.
% Boston College


if (dr=='B') % B(t)
    if (t>=s-2*h && t<s-h)
        bout=((t-s+2*h)/h)^3/6;
    elseif (t>=s-h && t<s)
        bout=1.0/6+(t-s+h)/h/2+((t-s+h)/h)^2/2-((t-s+h)/h)^3/2;
    elseif (t>=s && t<s+h)
        bout=1.0/6+(s+h-t)/h/2+((s+h-t)/h)^2/2-((s+h-t)/h)^3/2;
    elseif (t>=s+h && t<s+2*h)
        bout=((s+2*h-t)/h)^3/6;
    else
        bout=0;
    end;
    
elseif (dr=='dB') % dB(t)/dt
    if (t>=s-2*h && t<s-h)
        bout=(t-s+2*h)^2/h^3/2;
    elseif (t>=s-h && t<s)
        bout=1/(2*h)+(t-s+h)/h^2-3*(t-s+h)^2/h^3/2;
    elseif (t>=s && t<s+h)
        bout=(-1)/(2*h)-(s+h-t)/h^2+3*(s+h-t)^2/h^3/2;
    elseif (t>=s+h && t<s+2*h)
        bout=(-1)*(s+2*h-t)^2/h^3/2;
    else
        bout=0;
    end;
end

return
