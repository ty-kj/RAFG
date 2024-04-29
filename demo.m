addpath(genpath(pwd));

%RAFG-UFS
clear
clc
load('./datasets/jaffe.mat');
% NewFeaNum=[50 100 150 200 250 300];
NewFeaNum=[20 30 40 50 60 70 80 90 100];
nc=length(unique(gnd));
alpha=1e-1;
beta = 1e2;
range=[1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2,1e3,1e4,1e5];
prange=[0.01,0.1,0.4,0.7,1];
model=0;
for i = 1:length(range)
    alpha = range(i);
    for j = 11:length(range)
        beta = range(j);
        for k = 1:length(prange)
            p = prange(k);
            model = model+1
            [idx] = RAFG(fea',alpha,beta,nc,p);

            for m=1:size(NewFeaNum,2)
                Newfea=fea(:,idx(1:NewFeaNum(m)));
                label=litekmeans(Newfea,nc,'Replicates',20);
                [Acc(m),Nmi(m),~]=ClusteringMeasure(gnd,label);
            end
            meanacc(model)=mean(Acc);
            meannmi(model)=mean(Nmi);
            stdacc(model)=std(Acc);
            stdnmi(model)=std(Nmi);  
            palpha(model)=alpha;
            pbeta(model)=beta;
        end
    end
end
% clearvars -except meanacc meannmi stdacc stdnmi
