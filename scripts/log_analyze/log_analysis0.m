clear;
clc;
close all;


matTmp=load('enfwi_noise.txt');
numSample=size(matTmp,2);
indAnalysisStep=15;
indBaseLine=(numSample+3)*(indAnalysisStep-1);

vecL1=matTmp(indBaseLine+1,:)';
vecL2=matTmp(indBaseLine+2,:)';
vecResidual=matTmp(indBaseLine+3,:)';
matUpd=matTmp(indBaseLine+4:indBaseLine+3+numSample,:);


%L1
vecCoeSumL1=sum(matUpd,1);
matTmpCov=cov(vecCoeSumL1,vecL1);
tmpSumCovL1=matTmpCov(2,1)/sqrt(matTmpCov(1,1)*matTmpCov(2,2));
fprintf('Sum L1 cov: %f\n',tmpSumCovL1);
vecCovL1=zeros(numSample,1);
for indSample=1:numSample
    matTmpCov=cov(matUpd(indSample,:),vecL1);
    vecCovL1(indSample)=matTmpCov(2,1)/sqrt(matTmpCov(1,1)*matTmpCov(2,2));
end
figure;
subplot(3,1,1);
plot(vecL1);
colorbar();
title('L1');
subplot(3,1,2);
imagesc(matUpd);
colorbar();
title('matUpd');
subplot(3,1,3);
plot(vecCovL1);
colorbar();
title('Cov of the coefficients and L1 for each sample');

%L2
vecCoeSumL2=sum(matUpd,1);
matTmpCov=cov(vecCoeSumL2,vecL2);
tmpSumCovL2=matTmpCov(2,1)/sqrt(matTmpCov(1,1)*matTmpCov(2,2));
fprintf('Sum L2 cov: %f\n',tmpSumCovL2);
vecCovL2=zeros(numSample,1);
for indSample=1:numSample
    matTmpCov=cov(matUpd(indSample,:),vecL2);
    vecCovL2(indSample)=matTmpCov(2,1)/sqrt(matTmpCov(1,1)*matTmpCov(2,2));
end
figure;
subplot(3,1,1);
plot(vecL2);
colorbar();
title('L2');
subplot(3,1,2);
imagesc(matUpd);
colorbar();
title('matUpd');
subplot(3,1,3);
plot(vecCovL2);
colorbar();
title('Cov of the coefficients and L2 for each sample');

%Residual
vecCoeSumResidual=sum(matUpd,1);
matTmpCov=cov(vecCoeSumResidual,vecResidual);
tmpSumCovResidual=matTmpCov(2,1)/sqrt(matTmpCov(1,1)*matTmpCov(2,2));
fprintf('Sum Residual cov: %f\n',tmpSumCovResidual);
vecCovResidual=zeros(numSample,1);
for indSample=1:numSample
    matTmpCov=cov(matUpd(indSample,:),vecResidual);
    vecCovResidual(indSample)=matTmpCov(2,1)/sqrt(matTmpCov(1,1)*matTmpCov(2,2));
end
figure;
subplot(3,1,1);
plot(vecResidual);
colorbar();
title('Residual');
subplot(3,1,2);
imagesc(matUpd);
colorbar();
title('matUpd');
subplot(3,1,3);
plot(vecCovResidual);
colorbar();
title('Cov of the coefficients and residual for each sample');