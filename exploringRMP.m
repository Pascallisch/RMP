%Replace nans in tags (we checked this is legit - these people have the
%tags they should have)
%TAGS(1501553,1) = 0; 
%NAMES{295559,1} = 0;
figureLabels = {'TG','GF','R','RR','PM','SC','HW','I','PQ','A','MP','CG','H','TH','FT','AL','C','EC','GP','LH'};

temp = find(isnan(TAGS)==1);
TAGS(temp) = 0;
bs = 1; %Do the bootstrap
bsCI = 99; %Get the 99% confidence interval
%% 3 ways of doing this:
%1) Just use the raw number of tags. Problem: People will have more tags
%just because they have more ratings.
%2) Divide the number of tags by number of ratings. Problem: The tags were
%introduced midway through (can you find out when exactly?)
%If we only look at profiles after TAGS were introduced this is the way to
%go
%3) Normalize so that all TAGS add up to 1.
%Problem: If someone has a single tag - just 1 - it will look like they are only that. Distorting it. Actually, that's probably the only way to get numbers above 0.3
%4) Only take tags after January 2015
%5) Do 3) - to normalize. But make sure that there are at least 10 tags.

%Normalize the data, but with more than 15 tags
temp = sum(TAGS'); %Determine the sum of the tags
%Pick those above the cutoff
tagCutoff = 15; 
indicesforTagAnalysis = find (temp>tagCutoff);
DATAtags = DATA(indicesforTagAnalysis,:);
NAMEStags = NAMES(indicesforTagAnalysis,:);
TAGStags = TAGS(indicesforTagAnalysis,:); 

%Find male and females indices
maleIndices = find(DATAtags(:,7)==1);
femaleIndices = find(DATAtags(:,8)==1); 

minimin = -2;
maximax = 2;
textOffset = 0.25;

nTags = nan(length(TAGStags),20); %Preallocate
for ii = 1:length(TAGStags)
nTags(ii,:) = TAGStags(ii,:)./DATAtags(ii,3);%sum(TAGStags(200,:)); %Normalize
if mod(ii,1e4) == 0
    ii
end
end
% Regression models
%a) Overall
[betas,betaCI,resid,RINT,STATS] = regress(DATAtags(:,1),[ones(length(DATAtags),1) nTags]);
orderedBetas = [[1:20]' betas(2:end) betaCI(2:end,:)];
orderedBetas = sortrows(orderedBetas,2);

eBars = abs(orderedBetas(:,2)-orderedBetas(:,3));

[betas2,betaCI2,resid2,RINT2,STATS2] = regress(DATAtags(:,2),[ones(length(DATAtags),1) nTags]);
orderedBetas2 = [[1:20]' betas2(2:end) betaCI2(2:end,:)];
orderedBetas2 = sortrows(orderedBetas2,2);

eBars2 = abs(orderedBetas2(:,2)-orderedBetas2(:,3));

figure
subplot(2,1,1)
x = 1:length(orderedBetas);
miniMin = min(x)-1;
maxiMax = max(x)+1;
h1 = bar(x,orderedBetas(:,2));
for ii = x
    if orderedBetas(ii,2) < 0
text(ii,orderedBetas(ii,2)-textOffset,figureLabels{orderedBetas(ii,1)},'Horizontalalignment','center');
    else
text(ii,orderedBetas(ii,2)+textOffset,figureLabels{orderedBetas(ii,1)},'Horizontalalignment','center');        
    end
end
h1.FaceColor = [0.9 0.9 0.9];
hold on
h2 = errorbar(x,orderedBetas(:,2),eBars(:,1));
h2.Color = 'r';
h2.LineWidth = 3;
h2.LineStyle = 'none';
xlim([miniMin maxiMax])
%h3 = line([miniMin maxiMax],[0 0], 'linestyle', '--');
ylabel('Tag beta on quality')
ylim([minimin maximax])
set(gca,'XTick',[])
movshonize(26,1)
axis normal
title('A')
subplot(2,1,2)
h2 = bar(x,orderedBetas2(:,2));
for ii = x
    if orderedBetas2(ii,2) < 0
text(ii,orderedBetas2(ii,2)-textOffset,figureLabels{orderedBetas2(ii,1)},'Horizontalalignment','center');
    else
text(ii,orderedBetas2(ii,2)+textOffset,figureLabels{orderedBetas2(ii,1)},'Horizontalalignment','center');        
    end
end
h2.FaceColor = [0.9 0.9 0.9];
hold on
h2 = errorbar(x,orderedBetas2(:,2),eBars2(:,1));
h2.Color = 'r';
h2.LineWidth = 3;
h2.LineStyle = 'none';
xlim([miniMin maxiMax])
%h3 = line([miniMin maxiMax],[0 0], 'linestyle', '--');
ylabel('Tag beta on difficulty')
ylim([minimin maximax])
set(gca,'XTick',[])
movshonize(26,1)
axis normal
title('B')
makeWhite
print('fig3','-dpng')

%%
figure
%Quality
subplot(2,1,1) 

% b) Male only
[betas,betaCI,resid,RINT,STATS] = regress(DATAtags(maleIndices,1),[ones(length(maleIndices),1) nTags(maleIndices,:)]);
betaMale = betas(2:end);
eBarsMale = abs(betas(2:end)-betaCI(2:end,:));

x = 1:length(betas)-1;
miniMin = min(x)-1;
maxiMax = max(x)+1;
h1 = bar(x,betaMale(orderedBetas(:,1)));
for ii = x
    if orderedBetas(ii,2) < 0
text(ii,orderedBetas(ii,2)-textOffset,figureLabels{orderedBetas(ii,1)},'Horizontalalignment','center');
    else
text(ii,orderedBetas(ii,2)+textOffset,figureLabels{orderedBetas(ii,1)},'Horizontalalignment','center');        
    end
end
set(h1,'facecolor','b')
set(h1,'facealpha',0.5)

hodl
h2 = errorbar(x,betaMale(orderedBetas(:,1)),eBarsMale(orderedBetas(:,1),1),eBarsMale(orderedBetas(:,1),2));

h2.Color = 'k';
h2.LineWidth = 2;
h2.LineStyle = 'none';
xlim([miniMin maxiMax])
%h3 = line([miniMin maxiMax],[0 0], 'linestyle', '--');


% c) Female only
[betas,betaCI,resid,RINT,STATS] = regress(DATAtags(femaleIndices,1),[ones(length(femaleIndices),1) nTags(femaleIndices,:)]);
betaFemale = betas(2:end);
eBarsFemale = abs(betas(2:end)-betaCI(2:end,:));
x = 1:length(betas)-1;
miniMin = min(x)-1;
maxiMax = max(x)+1;
%h2 = bar(x,betas(2:end));
h2 = bar(x,betaFemale(orderedBetas(:,1)));
set(h2,'facecolor',[1 0 1]) %'r'
set(h2,'facealpha',0.5)
%h3 = errorbar(x,betas(2:end),eBarsFemale(:,1),eBarsFemale(:,2));
h3 = errorbar(x,betaFemale(orderedBetas(:,1)),eBarsFemale(orderedBetas(:,1),1),eBarsFemale(orderedBetas(:,1),2));
h3.Color = 'k';
h3.LineWidth = 2;
h3.LineStyle = 'none';
xlim([miniMin maxiMax])
%h3 = line([miniMin maxiMax],[0 0], 'linestyle', '--');
ylabel('Tag beta on quality')
set(gca,'XTick',[])
movshonize(26,1)
axis normal


%Difficulty
subplot(2,1,2) 

% b) Male only
[betas,betaCI,resid,RINT,STATS] = regress(DATAtags(maleIndices,2),[ones(length(maleIndices),1) nTags(maleIndices,:)]);
betaMale = betas(2:end);
eBarsMale = abs(betas(2:end)-betaCI(2:end,:));

x = 1:length(betas)-1;
miniMin = min(x)-1;
maxiMax = max(x)+1;
h1 = bar(x,betaMale(orderedBetas2(:,1)));
for ii = x
    if orderedBetas2(ii,2) < 0
text(ii,orderedBetas2(ii,2)-textOffset,figureLabels{orderedBetas2(ii,1)},'Horizontalalignment','center');
    else
text(ii,orderedBetas2(ii,2)+textOffset,figureLabels{orderedBetas2(ii,1)},'Horizontalalignment','center');        
    end
end
set(h1,'facecolor','b')
set(h1,'facealpha',0.5)

hodl
h2 = errorbar(x,betaMale(orderedBetas2(:,1)),eBarsMale(orderedBetas2(:,1),1),eBarsMale(orderedBetas2(:,1),2));

h2.Color = 'k';
h2.LineWidth = 2;
h2.LineStyle = 'none';
xlim([miniMin maxiMax])
%h3 = line([miniMin maxiMax],[0 0], 'linestyle', '--');


% c) Female only
[betas,betaCI,resid,RINT,STATS] = regress(DATAtags(femaleIndices,2),[ones(length(femaleIndices),1) nTags(femaleIndices,:)]);
betaFemale = betas(2:end);
eBarsFemale = abs(betas(2:end)-betaCI(2:end,:));
x = 1:length(betas)-1;
miniMin = min(x)-1;
maxiMax = max(x)+1;
%h2 = bar(x,betas(2:end));
h2 = bar(x,betaFemale(orderedBetas2(:,1)));
set(h2,'facecolor',[1 0 1]) %'r'
set(h2,'facealpha',0.5)
%h3 = errorbar(x,betas(2:end),eBarsFemale(:,1),eBarsFemale(:,2));
h3 = errorbar(x,betaFemale(orderedBetas2(:,1)),eBarsFemale(orderedBetas2(:,1),1),eBarsFemale(orderedBetas2(:,1),2));
h3.Color = 'k';
h3.LineWidth = 2;
h3.LineStyle = 'none';
xlim([miniMin maxiMax])
%h3 = line([miniMin maxiMax],[0 0], 'linestyle', '--');
ylabel('Tag beta on difficulty')
set(gca,'XTick',[])
movshonize(26,1)
axis normal


makeWhite

print('fig4','-dpng')

% %%
% figure
% diffWeights = [[1:20]' abs(betaMale)-abs(betaFemale)];
% diffWeights2 = sortrows(diffWeights,-2);
% 
% figureLabels = {'TG','GF','R','RR','PM','SC','HW','I','PQ','A','MP','CG','H','TH','FT','AL','C','EC','GP','LH'};
% bar(1:20,diffWeights2(:,2))
% for ii = 1:20
%     if diffWeights2(ii,2) > 0
% text(ii,diffWeights2(ii,2)+0.02,figureLabels(diffWeights2(ii,1)),'HorizontalAlignment','center');
%     else
% text(ii,diffWeights2(ii,2)-0.02,figureLabels(diffWeights2(ii,1)),'HorizontalAlignment','center');        
%     end
% end
% title('Absolute difference in tag weight, male vs. female')
% ylim([-0.25 0.15])
% movshonize(26,1)
% axis normal
% makeWhite
%% Intercorrelation analysis of tags 
figure
CCF1 = corrcoef(nTags);
subplot(2,2,1)
colormap(jet)
CCF1(find(CCF1==1)) = 0;
imagesc(CCF1); colorbar; 
axis square
title('Overall')

CCF2 = corrcoef(nTags(maleIndices,:));
subplot(2,2,3)

CCF2(find(CCF2==1)) = 0;
imagesc(CCF2); colorbar; 
axis square
title('Male')
colormap(jet)


CCF3 = corrcoef(nTags(femaleIndices,:));
subplot(2,2,4)
CCF3(find(CCF3==1)) = 0;
imagesc(CCF3); colorbar; 
axis square
title('Female')
colormap(jet)

subplot(2,2,2)
CCF4 = CCF2-CCF3;
imagesc(CCF4); colorbar; 
axis square
title('Difference - male minus female')
shg
colormap(jet)

sum(abs(CCF4))

%% Tag differences
for ii = 1:20
MD(ii,1) = mean(nTags(:,ii));
MD(ii,2) = mean(nTags(maleIndices,ii));
MD(ii,3) = mean(nTags(femaleIndices,ii));
MD(ii,4) = MD(ii,2)-MD(ii,3);
MD(ii,5) = ii;
end

figure
figureLabels = {'TG','GF','R','RR','PM','SC','HW','I','PQ','A','MP','CG','H','TH','FT','AL','C','EC','GP','LH'};
sortedEffects = sortrows(MD,-4);
hold on

for ii = 1:20
    if sortedEffects(ii,4) > 0
text(ii,sortedEffects(ii,4)+0.005,figureLabels(sortedEffects(ii,5)),'HorizontalAlignment','center');
bar(ii,sortedEffects(ii,4),'facecolor','b','facealpha',0.5)
    else
text(ii,sortedEffects(ii,4)-0.005,figureLabels(sortedEffects(ii,5)),'HorizontalAlignment','center');        
bar(ii,sortedEffects(ii,4),'facecolor',[1 0 1],'facealpha',0.5)
    end
end
title('Difference in mean tag prevalence, male vs. female')
ylim([-0.045 0.075])
movshonize(26,1)
axis normal
set(gca,'XTick',[])

%%
if bs == 1
% Adding the bootstrapped confidence intervals to this
nReps = 1e4;
nIndices = length(nTags);
nMale = length(maleIndices);
nFemale = length(femaleIndices);

bMD = nan(nReps,20,5); %Preallocate
for nn = 1:nReps
tempOrder = randperm(nIndices); %Shuffle all
bMaleIndices = tempOrder(1:nMale);
bFemaleIndices = tempOrder(nMale+1:nMale+nFemale); 

for ii = 1:20
bMD(nn,ii,1) = mean(nTags(:,ii));
bMD(nn,ii,2) = mean(nTags(bMaleIndices,ii));
bMD(nn,ii,3) = mean(nTags(bFemaleIndices,ii));
bMD(nn,ii,4) = bMD(nn,ii,2)-bMD(nn,ii,3);
bMD(nn,ii,5) = ii;
end    
nn
end
end
%% Hardcoding in the interest of time
lowerBoundIndex = 50;
higherBoundIndex = 9950; 
for ii = 1:20
    tempSort = sortrows(bMD(:,ii,4));
    bounds(ii,1) = tempSort(lowerBoundIndex);
    bounds(ii,2) = tempSort(higherBoundIndex);
end

%%
hq = errorbar(1:20,zeros(1,20),bounds(sortedEffects(:,5),1),bounds(sortedEffects(:,5),2))
hq.Color = 'k';
hq.LineStyle = 'none';
hq.LineWidth = 3;

makeWhite
print('fig5','-dpng')

%% Logistic regression
[B,q,stats] = mnrfit(nTags,DATA(:,4)+1);
phat = mnrval(B,nTags)