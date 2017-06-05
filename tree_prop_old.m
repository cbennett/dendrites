function [conductancemat,outtr]  = tree_prop_old(n,c,ind,sum,junctures,levels,levelsback,andnumb,isInhibited,inputspikes,dim,mat,syn,conductancemat,outtr,assigntable)
		   %nanosyn stuff
		    vst = 0.5;
			vplus = 3.5;
			vminus = 6.5;
			vt1 = 3.1;
			vt2 = 6.0;
			writelev = 256;
			gmax = 6e-3;
			gmin = 1e-6;
            bias = 0.1;
             % vectors to clear each time
		    coinc = zeros(sum,1);
            current = zeros(sum,1);
            currentdiff = zeros(junctures,1);
            spikes = zeros(sum,1);
            gates = zeros(junctures,1);
            deltamat = zeros(sum,syn);
            active = zeros(sum,syn);
            for k = 1:levels
                if k == 1 % synapses affected directly by retina
                    first = 2^(levels);
                    jxnret = first/2;
                    %potentiate synapses and produce output spikes
                    for f = 1:first
                        branchcons = zeros(syn,1);
                        branchcons(:) = mat(n,f,:);
                        for dd = 1:dim
                            if inputspikes(ind,dd) == 1
                                connected = find(branchcons,dd);
                                %active(f,:) = connected;
                                length = size(connected);
                                for ii = 1:length
                                    active(f,connected(ii)) = 1;
                                    if connected(ii) == dd
                                        %conductancemat(n,f,i(ii)) = condevolvein(conductancemat(n,f,i(ii)), vin,vt1,vt2,deltag,gmax,gmin);
                                        current(f) = current(f) + conductancemat(f,connected(ii),n)*vst;
                                        coinc(f) = coinc(f) + 1;
                                    end
                                end
                            end
                        end
                        if coinc(f) > (andnumb-1)
                            spikes(f) = spikes(f) + 1;
                        end
                    end
                    % compute gates
                    for ff = 1:jxnret
                        currentdiff(ff) = current(2*ff)- current(2*ff-1);
                        if spikes(2*ff) ==1 || spikes(2*ff-1) == 1  && currentdiff(ff) > 0 && isInhibited(ff,n) ~= 1
                            gates(ff) = gates(ff) + 1;
                        end
                    end
                    if levels == levelsback
                        for ff = 1:jxnret
                            if gates(ff) == 1
                                deltamat((2*ff-1),:) = 1*active(ff,:); % pos branch
                                deltamat(2*ff,:) = -1*active(ff,:); % neg branch
                            else
                                deltamat((2*ff-1),:) = -1*active(ff,:); % pos branch
                                deltamat(2*ff,:) = 1*active(ff,:); % neg branch
                            end
                        end
                    end
                elseif k< levels
                    if k == 2
                        jxnfirst = jxnret + 1;
                        jxnlast = jxnret + 2^(levels-k+1)/2;
                    else
                        jxnfirst = jxnlast + 1;
                        jxnlast = jxnfirst + 2^(levels-k+1)/2-1;
                    end
                    % forward prop & compute gates
                    first = jxnfirst*2 - 1;
                    last = jxnlast*2 + 1;
                    for ff = first:last
                        %potentiate synapses and produce output spikes
                        for dd = 1:dim
                              if inputspikes(ind,dd) == 1
                                 X = mat(n,ff,:);
                                 connected = find(X,dd);
                                 length = size(connected);
                                 for ii = 1:length
                                    active(ff,connected(ii)) = 1;
                                    if connected(ii) == dd
                                         current(ff) = current(ff) + conductancemat(ff,connected(ii),n)*vst;
                                         coinc(ff) = coinc(ff) + 1;
                                    end
                                 end
                              end 
                        end
                        if coinc(ff) > (andnumb-1)
                            spikes(ff) = spikes(ff) + 1;
                        end  
                    end
                    for gg = jxnfirst:jxnlast
                        index = 2*gg;
                        if spikes(index) ==1 || spikes(index-1) == 1 && (current(index)-current(index-1) > 0) && isInhibited(gg,n) ~= 1
                            gates(gg) = gates(gg) + 1;
                        end
                    end
                    if k >= (levels - levelsback) %supervised!
                        for ff = jxnfirst:jxnlast
                            if gates(ff) == 1
                                deltamat((2*ff-1),:) = 1*active(ff,:); % pos branch
                                deltamat(2*ff,:) = -1*active(ff,:); % neg branch
                            else
                                deltamat((2*ff-1),:) = -1*active(ff,:); % pos branch
                                deltamat(2*ff,:) = 1*active(ff,:); % neg branch
                            end
                        end
                    end
                else % Ultimate branch:look one gate back and adjust weights
                    if gates(jxnfirst) == 1 || gates(jxnlast) == 1 && (current(jxnfirst)-current(jxnlast) > 0) && isInhibited(junctures,n) ~= 1
                        gates(junctures) = 1;
                        outtr(n,ind) = outtr(n,ind) + 1;
                        deltamat(sum,:) = 1;
                    else
                        deltamat(sum,:) = -1;
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
                    %if outtr(n,ind) == 0 && assigntable(c,n)== 1
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
                            len = size(X3);
                            if len > 0
                                 for i = 1:len
                                    conductancemat(f,X3(i),n) = condevolvein(conductancemat(f,X3(i),n), vminus,vt1,vt2,writelev,gmax,gmin,1);
                                end
                            end
                        end
                    elseif outtr(n,ind) == 1 && assigntable(c,n) == 1 % REWARD
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
                            len = size(X3);
                            if len > 0 
                                for i = 1:len
                                    conductancemat(f,X3(i),n) = condevolvein(conductancemat(f,X3(i),n), vplus,vt1,vt2,writelev,gmax,gmin,2);
                                end
                            end
                        end
                        
                    end
                end
            end
            %coinc
            %current
            %gates
            %spikes
            %deltamat
 end