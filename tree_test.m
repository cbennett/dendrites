function outte  = tree_test(n,ind,sum,junctures,levels,isInhibited,inputspikeste,dim,mat,syn,conductancemat,outte)
vst = 0.5;
bias = 0.01;
% vectors to clear each time
gates = zeros(junctures,1);
current = zeros(sum,1);
active = zeros(sum,syn);
branchcons = zeros(syn,1);
for k = 1:levels
    if k == 1 % synapses affected directly by retina
        first = 2^(levels);
        jxnret = first/2;
        for f = 1:first
            branchcons(:) = mat(n,f,:);
            active(f,:) = sign(branchcons);
            for dd = 1:dim
                if inputspikeste(ind,dd) == 1
                    activex = find(branchcons > 0);
                    length = size(activex,1);
                    if length > 0
                        for ii = 1:length
                            if abs(branchcons(activex(ii))) == dd
                                current(f) = current(f) + conductancemat(f,activex(ii),n)*vst;
                            end
                        end
                    end
                    activei = find(branchcons < 0);
                    length2 = size(activei,1);
                    if length2 > 0
                        for ii = 1:length2
                            if abs(branchcons(activei(ii))) == dd
                                current(f) = current(f) - (conductancemat(f,activei(ii),n)*vst);
                            end
                        end
                    end
                end
            end
        end
        % compute gates
        for ff = 1:jxnret
            if current(2*ff) > 0 && current(2*ff-1) > 0 && isInhibited(ff,n) ~= 1
                gates(ff) = gates(ff) + 1;
            end
        end
    elseif k < levels % forward prop
        jxnit = 0;
        for i  = 1:(k-1)
            jxnit = jxnit + (2^(levels-i+1))/2;
        end
        jxnfirst = jxnit + 1;
        jxnlev = 2^(levels-k+1)/2;
        jxnlast = jxnfirst + jxnlev -1;
        first = jxnfirst*2 - 1;
        last = jxnlast*2 + 1;
        offset = 2^levels;
        for ff = first:last
            branchcons(:) = mat(n,ff,:);
            %active(ff,:) = sign(branchcons);
            for dd = 1:dim
                if inputspikeste(ind,dd) == 1
                    activex = find(branchcons > 0);
                    length = size(activex,1);
                    if length > 0
                        for ii = 1:length
                            if abs(branchcons(activex(ii))) == dd
                                current(ff) = current(ff) + conductancemat(ff,activex(ii),n)*vst;
                            end
                        end
                    end
                    activei = find(branchcons < 0);
                    length2 = size(activei,1);
                    if length2 > 0
                        for ii = 1:length2
                            if abs(branchcons(activei(ii))) == dd
                                current(ff) = current(ff) - (conductancemat(ff,activei(ii),n)*vst);
                            end
                        end
                    end
                end
            end
            if gates(ff - offset) > 0
               current(ff) = current(ff) + bias;
            end
        end
        for gg = jxnfirst:jxnlast
        index = 2*gg;
            if current(index) > 0 && current(index-1) > 0 && isInhibited(gg,n) ~= 1
                gates(gg) = gates(gg) + 1;
            end
        end    
    else
        branchcons(:) = mat(n,sum,:);
        for dd = 1:dim
            if inputspikeste(ind,dd) == 1
                activex = find(branchcons > 0);
                length = size(activex,1);
                if length > 0 
                    for ii = 1:length
                        if abs(branchcons(activex(ii))) == dd
                            current(sum) = current(sum) + conductancemat(sum,activex(ii),n)*vst;
                        end
                    end
                end
                activei = find(branchcons < 0);
                length2 = size(activei,1);
                if length2 > 0
                    for ii = 1:length2
                        if abs(branchcons(activei(ii))) == dd
                            current(sum) = current(sum) - (conductancemat(sum,activei(ii),n)*vst);
                        end
                    end
                end
            end
        end
        if gates(junctures - 1) > 0
            current(sum) = current(sum) + bias;
        end
        if  current(sum) > 0 && isInhibited(junctures,n) ~= 1
            gates(junctures) = 1;
            outte(n,ind) = outte(n,ind) + 1;
        end
    end
    
end