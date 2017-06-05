clear all
%setup db - USPS
load('usps_all.mat');
[trainind,valind,testind] = dividerand(1100,0.7,0.05,0.25);
K = 10; %classes
dim = 256;
training = length(trainind);
trimages = zeros(dim,training,K);
%trtotal = training;
for k = 1:K
    for j = 1:training
        trimages(:,j,k) = data(:,j,k);
    end
end
tests = length(testind);
teimages = zeros(dim,tests,K);
for k = 1:K
    for j = 1:tests
        teimages(:,j,k) = data(:,j,k);
    end
end
%set up neurons and trees
levels = 4;
levelsback = 4; % number of levels  trained in supervised fashion by by the backprop spike
junctures = 0;
sum = 0;
for i  = 1:levels
    sum = sum + 2^(levels-i+1);
    junctures = junctures + (2^(levels-i+1))/2;
end
nbtrees = sum; %number branches each tree
syn = 32; % number of synapses/tree
nbneurons = 10*K; % number of neurons in pop; should be some multiple of K
tot = nbtrees*nbneurons*syn;
spars = 0.2; % probability of any given syn being connected to an input
mat = assigntrees(dim,tot,spars,nbneurons,nbtrees,syn);
[ne,bra,syn] = size(mat); % m = neurons ; n = number branches per tree; l = syn per branch
conductancemat = 0.0008*rand(bra,syn,ne);
deltamat = zeros(bra,syn);
andnumb = (syn/10)-1;
trtotal = 2*K;
inputspikes = zeros(trtotal*K,dim);
outtr = zeros(nbneurons,K*trtotal);
expected  = zeros(K,K*trtotal);
assigntable = zeros(K,nbneurons);
for i = 1:nbneurons
    n = randi(10,1);
    assigntable(n,i) =  assigntable(n,i) + 1;
end
inhvec = 0:0.05:1;
results = zeros(length(inhvec),1);
for ii = 1:length(inhvec)
    inhibProb = inhvec(ii);
    inhibProb
    isInhibited = zeros(junctures,nbneurons);
    for n = 1:nbneurons
          for ff = 1:junctures
              randn = rand;
              if randn < (0.5 - (1-inhibProb)/2) || randn > (0.5 + (1-inhibProb)/2)
                    isInhibited(ff,n) = 1;
              end
          end 
    end
    for tr = 1:trtotal %epochs
        tr
        for c = 1:K
            image = trimages(:,tr,c);
            ind = (tr -1)*K + c;
            expvec = zeros(K,1);
            expvec(c) = 1;
            expected(:,ind) = expvec;
            for k = 1:dim % binarized into spike
                if image(k) > 0
                    inputspikes(ind,k) = 1;
                end
            end
            for n = 1:nbneurons
               [conductancemat,outtr] = tree_prop(n,c,ind,sum,junctures,levels,levelsback,isInhibited,inputspikes,dim,mat,syn,conductancemat,outtr,assigntable);
              %[conductancemat,outtr] = tree_prop_old(n,c,ind,sum,junctures,levels,levelsback,andnumb,isInhibited,inputspikes,dim,mat,syn,conductancemat,outtr,assigntable); 
            end
        end
    end
    Wout =  expected*outtr'*inv(outtr*outtr'+eye(nbneurons)); % currents bias
    %test patterns
    testsgiven = round(tests/20);
    prediction = zeros(K,testsgiven);
    inputspikeste = zeros(testsgiven*K,dim);
    outte = zeros(nbneurons,K*testsgiven);
    actual = zeros(K,testsgiven*K);
    found = 0;
    for te = 1:testsgiven
        te
        for c = 1:K
            image = teimages(:,te,c);
            ind = (te -1)*K + c;
            for k = 1:dim % binarized into spike
                if image(k) > 0
                    inputspikeste(ind,k) = 1;
                end
            end
            labels = zeros(K,1);
            labels(c) = 1;
            actual(:,ind) = labels;
            for n = 1:nbneurons
                outte  = tree_test(n,ind,sum,junctures,levels,levelsback,isInhibited,inputspikeste,dim,mat,syn,conductancemat,outte);
                %outte  = tree_test_old(n,ind,sum,junctures,levels,levelsback,andnumb,isInhibited,inputspikeste,dim,mat,syn,conductancemat,outte);
            end
            outvec = outte(:,ind);
            output = outvec'*Wout.';
            [Max,I] = max(output);
            if I(1) == c
                found = found + 1;
            end
        end
    end
    percentage = (found/(testsgiven*K))*100;
    finalmsg = ['Final test percentage is:',num2str(percentage) '%'];
    disp(finalmsg)
    results(ii) = percentage;
    plot(inhvec,results)
end