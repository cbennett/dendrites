function [conductancemat,outtr]  = tree_prop(n,c,ind,sum,junctures,levels,levelsback,isInhibited,inputspikes,dim,mat,syn,conductancemat,outtr,assigntable)
vst = 0.5;
vplus = 3.5;
vminus = 6.5;
vt1 = 3.1;
vt2 = 6.0;
writelev = 64;
gmax = 6e-3;
gmin = 1e-6;
bias = 0.01;
% vectors to clear each time
current = zeros(sum,1);
gates = zeros(junctures,1);
deltamat = zeros(sum,syn);
branchcons = zeros(syn,1);
active = -1*ones(sum,syn);
for k = 1:levels
    if k == 1 % synapses affected directly by retina
        first = 2^(levels);
        jxnret = first/2;
        %potentiate synapses and produce output spikes
        for f = 1:first
            branchcons(:) = mat(n,f,:);
            for dd = 1:dim
                if inputspikes(ind,dd) == 1
                    activex = find(branchcons > 0);
                    length = size(activex,1);
                    if length > 0
                        for ii = 1:length
                            if abs(branchcons(activex(ii))) == dd
                               active(f,activex(ii)) = 1;
                               current(f) = current(f) + conductancemat(f,activex(ii),n)*vst;
                            end
                        end
                    end
                    activei = find(branchcons < 0);
                    length2 = size(activei,1);
                    if length2 > 0 
                        for ii = 1:length2
                            if abs(branchcons(activei(ii))) == dd
                               active(f,activei(ii)) = 1;
                               current(f) = current(f) - (conductancemat(f,activei(ii),n)*vst);
                            end
                        end
                    end
                end
            end
        end
        % compute gates
        for ff = 1:jxnret
            if current(2*ff) > 0  && current(2*ff-1) > 0 && isInhibited(ff,n) ~= 1
                gates(ff) = gates(ff) + 1;
            end
        end
        if levels == levelsback
            for ff = 1:jxnret
                if gates(ff) == 1
                    deltamat((2*ff-1),:) = active(2*ff-1,:); % pos branch
                    deltamat(2*ff,:) =  active(2*ff,:); % neg branch
                    %else
                    %    deltamat((2*ff-1),:) = -1*active(ff,:); % pos branch
                    %    deltamat(2*ff,:) = 1*active(ff,:); % neg branch
                end
            end
        end
    elseif k< levels
        jxnit = 0;
        for i  = 1:(k-1)
            jxnit = jxnit + (2^(levels-i+1))/2;
        end
        jxnfirst = jxnit + 1;
        jxnlev = 2^(levels-k+1)/2;
        jxnlast = jxnfirst + jxnlev -1;
        % forward prop & compute gates
        first = jxnfirst*2 - 1;
        last = jxnlast*2 + 1;
        offset = 2^(levels);
        for ff = first:last
            branchcons(:) = mat(n,ff,:);
            for dd = 1:dim
                if inputspikes(ind,dd) == 1
                    activex = find(branchcons > 0);
                    length = size(activex,1);
                    if length > 0
                       for ii = 1:length
                           if abs(branchcons(activex(ii))) == dd
                              active(ff,activex(ii)) == 1;
                              current(ff) = current(ff) + conductancemat(ff,activex(ii),n)*vst;
                           end
                       end
                    end
                    activei = find(branchcons < 0);
                    length2 = size(activei,1);
                    if length2 > 0
                       for ii = 1:length2
                           if abs(branchcons(activei(ii))) == dd
                              active(ff,activei(ii)) == 1;
                              current(ff) = current(ff) - conductancemat(ff,activei(ii),n)*vst;
                           end
                       end
                    end
                end
            end
            if gates(ff - offset) > 0
               current(ff) = current(ff) + bias;
            end
            %if current(ff) > 0
            %   active(ff,:) = sign(branchcons);
            %end
        end
        for gg = jxnfirst:jxnlast
            index = 2*gg;
            if current(index) > 0 && current(index-1) > 0 && isInhibited(gg,n) ~= 1
                gates(gg) = gates(gg) + 1;
            end
        end
        if k >= (levels - levelsback) %supervised!
            for ff = jxnfirst:jxnlast
                if gates(ff) == 1
                    deltamat((2*ff-1),:) = active(2*ff-1,:); % pos branch
                    deltamat(2*ff,:) = active(2*ff,:); % neg branch
                    %else
                    %    deltamat((2*ff-1),:) = -1*active(ff,:); % pos branch
                    %    deltamat(2*ff,:) = 1*active(ff,:); % neg branch
                end
            end
        end
    else % Ultimate branch:look one gate back and adjust weights
        branchcons(:) = mat(n,sum,:);
        for dd = 1:dim
            if inputspikes(ind,dd) == 1
                activex = find(branchcons > 0);
                length = size(activex,1);
                if length > 0
                    for ii = 1:length
                        if abs(branchcons(activex(ii))) == dd
                            active(sum,activex(ii)) = 1;
                            current(sum) = current(sum) + conductancemat(sum,activex(ii),n)*vst;
                        end
                    end
                end
                activei = find(branchcons < 0);
                length2 = size(activei,1);
                if length2 > 0
                   for ii = 1:length2
                       if abs(branchcons(activei(ii))) == dd
                          active(sum,activei(ii)) = 1;
                          current(sum) = current(sum) - (conductancemat(sum,activei(ii),n)*vst);
                       end
                   end
                end
            end
        end
        if gates(junctures-1) == 1
           current(sum) = current(sum) + bias;
           deltamat(sum,:) = active(sum,:);
        end 
        if current(sum) > 0 &&  isInhibited(junctures,n) ~= 1
            %active(sum,:) = sign(branchcons);
            gates(junctures) = 1;
            outtr(n,ind) = outtr(n,ind) + 1;
        end
        if outtr(n,ind) == 1 && assigntable(c,n) ~= 1 % FALSE POS
            %HL case > dpunish the synapses which caused spikes
            for f = 1:sum
                X2 = find(deltamat(f,:) == 1);
                len = size(X2,2);
                if len > 0
                    for i = 1:len
                        conductancemat(f,X2(i),n) = condevolvein(conductancemat(f,X2(i),n), vminus,vt1,vt2,writelev,gmax,gmin,1);
                    end
                end
                X3 = find(deltamat(f,:) == -1);
                len = size(X3,2);
                if len > 0
                    for i = 1:len
                        conductancemat(f,X3(i),n) = condevolvein(conductancemat(f,X3(i),n), vplus,vt1,vt2,writelev,gmax,gmin,1);
                    end
                end
            end
        elseif outtr(n,ind) == 0 && assigntable(c,n)== 1 % FALSE NEG
        %    %if outtr(n,ind) == 0 && assigntable(c,n)== 1
           % no spike & WAS expected : weakly strengthen path
            for f = 1:sum
                X2 = find(deltamat(f,:) == 1);
                len = size(X2,2);
                if len > 0
                    for i = 1:len
                        conductancemat(f,X2(i),n) = condevolvein(conductancemat(f,X2(i),n), vplus,vt1,vt2,writelev,gmax,gmin,1);
                    end
                end
                X3 = find(deltamat(f,:) == -1);
               len = size(X3,2);
                if len > 0
                    for i = 1:len
                        conductancemat(f,X3(i),n) = condevolvein(conductancemat(f,X3(i),n), vminus,vt1,vt2,writelev,gmax,gmin,1);
                    end
                end
            end
        elseif outtr(n,ind) == 1 && assigntable(c,n) == 1 % REWARD
        %if outtr(n,ind) == 1 && assigntable(c,n) == 1 % REWARD
            %if outtr(n,ind) == 1 && assigntable(c,n) == 1 % corr + anti-corr
            % most important step! strongly potentiate correlation . Also weaken all ones not contr (slightly)
            for f = 1:sum
                X2 = find(deltamat(f,:) == 1);
                len = size(X2,2);
                if len > 0
                    for i = 1:len
                        conductancemat(f,X2(i),n) = condevolvein(conductancemat(f,X2(i),n), vplus,vt1,vt2,writelev,gmax,gmin,2);
                    end
                end
                X3 = find(deltamat(f,:) == -1);
                len = size(X3,2);
                if len > 0
                    for i = 1:len
                        conductancemat(f,X3(i),n) = condevolvein(conductancemat(f,X3(i),n), vplus,vt1,vt2,writelev,gmax,gmin,2);
                    end
                end
            end
            
        end
    end
end
%deltamat
%gates
%current
end