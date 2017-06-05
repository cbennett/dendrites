function  synmat = assigntrees(dim,tot,sp,nneuron,ntrees,nsyn,negperc) 
          synper = (sp*tot)/dim; % how many redundant per pixel
          synmat = zeros(nneuron,ntrees,nsyn); %to fill
          for i = 1:dim
              for j = 1:synper
                  neur = randi(nneuron);  % pick random neuron
                  tree = randi(ntrees);     %pick random tree
                  syn = randi(nsyn);
                  if synmat(neur,tree,syn) ~= 0; %check redundt
                     randn = rand;
                    if randn < (0.5 - (1-negperc)/2) || randn > (0.5 + (1-negperc)/2)
                        synmat(neur,tree,syn) = -i; %create connectivity
                    else 
                       synmat(neur,tree,syn) = i; %create connectivity
                    end
                  else  %again
                    neur = randi(nneuron);  % pick random neuron
                    tree = randi(ntrees);     %pick random tree
                    syn = randi(nsyn);
                    randn = rand;
                    if randn < (0.5 - (1-negperc)/2) || randn > (0.5 + (1-negperc)/2)
                        synmat(neur,tree,syn) = -i; %create connectivity
                    else 
                       synmat(neur,tree,syn) = i; %create connectivity
                    end 
                  end
              end
         end
end

% OLD
          %out=~rem(dim,nb)*dim/nb;
          %if out == 0
          %   syn1 = dim/nb;
          %   syn = floor(syn1);
          %else
          %   syn = dim/nb;
          %end
          %p = randperm(dim)';
          %counter = 0;
                    %for i = 1:nb
          %    counter = counter +1;
          %    first = (i-1)*syn +1;
           %   last = first + syn -1;
          %    synmat(:,i) = p(first:last);
         % end