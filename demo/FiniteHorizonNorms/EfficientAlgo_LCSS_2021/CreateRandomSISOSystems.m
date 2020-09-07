% We use the following state dimentions for this study.
NxAll = [1 10:10:200];
m = length(NxAll);
nPlant = 5;

% Set seed
rng(0);

%  We create 5 random SISO plants using rss for a given state dimention and
% average the computational time to represent a sample.
k = 0;
for i = 1:m
    for j = 1:nPlant
        k = k + 1;
        while true
            Gk = rss(NxAll(i),1,1);
            ginf = hinfnorm(Gk);
            if ginf >= 0.1 && ginf <= 100
                fact = 1/sqrt(ginf);
                Gall(:,:,k) = fact*Gk*fact; %#ok<SAGROW>
                break;
            end
        end
    end
end

% Save MAT file
filename = sprintf('LTIDataSet_nPlant%d',nPlant);
save(filename,'Gall','k','nPlant','NxAll');