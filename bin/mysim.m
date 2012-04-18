function dy = mysim(t,y)
%
% Dynamics of the ethanolamine glycerophospholipid lipid remodeling network
%
% Copyright Mar 2012, Lu Zhang, all rights reserved.
% Boston College

tempfilename = ['..' filesep 'tmp' filesep 'tp16c19154_9bed_4a92_9441_4b248cb094ee'];

load(tempfilename);

dy = zeros(num_sp,1);
for alpha = 1 : num_sp
    for k = 1 : num_para
        inx = find(Network(:,1)==alpha & Network(:,2)==k);
        if (~isempty(inx))
            dy(alpha) = dy(alpha) - X2(k)*y(alpha)*length(inx);
        end
        inx = find(Network(:,3)==alpha & Network(:,2)==k);
        if (~isempty(inx))
            dy(alpha) = dy(alpha) + X2(k)*sum(y(Network(inx,1)));
        end
    end
end
