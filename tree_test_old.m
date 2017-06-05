function outte  = tree_test_old(n,ind,sum,junctures,levels,levelsback,andnumb,isInhibited,inputspikeste,dim,mat,syn,conductancemat,outte)
		    vst = 0.5;
             % vectors to clear each time
		    coinc = zeros(sum,1);
            spikes = zeros(sum,1);
            gates = zeros(junctures,1);
            current = zeros(sum,1);
            currentdiff = zeros(junctures,1);
            for k = 1:levels
                if k == 1 % synapses affected directly by retina
                    first = 2^(levels);
                    jxnret = first/2;
                   %potentiate synapses and produce output spikes
                    for f = 1:first
                        branchcons = zeros(syn,1);
                        branchcons(:) = mat(n,f,:);
                        for dd = 1:dim
                            if inputspikeste(ind,dd) == 1
                                connected = find(branchcons,dd);
                                length = size(connected);
                                for ii = 1:length
                                    if connected(ii) == dd
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
                elseif k < levels % forward prop
                    if k == 2
                        jxnfirst = jxnret + 1;
                        jxnlast = jxnret + 2^(levels-k+1)/2;
                    else
                        jxnfirst = jxnlast + 1;
                        jxnlast = jxnfirst + 2^(levels-k+1)/2-1;
                    end
                     first = jxnfirst*2 - 1;
                    last = jxnlast*2 + 1;
                    for ff = first:last
                        %potentiate synapses and produce output spikes
                        for dd = 1:dim
                              if inputspikeste(ind,dd) == 1
                                 X = mat(n,ff,:);
                                 connected = find(X,dd);
                                 length = size(connected);
                                 for ii = 1:length
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
                else
                    if gates(jxnfirst) == 1 || gates(jxnlast) == 1 && (current(jxnfirst)-current(jxnlast) > 0) && isInhibited(junctures,n) ~= 1
                        gates(junctures) = 1;
                        outte(n,ind) = outte(n,ind) + 1;
                    end
                end
            end