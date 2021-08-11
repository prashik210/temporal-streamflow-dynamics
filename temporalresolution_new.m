function[thresholdoldtycru,cctycru,PY,networkdensity,FNNtycru,slopefinaltycru,critthresholdfinaltycru,critavgcctycru,globalefficiencytycru,harmonicmeantycru]=temporalresolution_new(ty)
% initialization of variables
thresholdoldtycru = zeros(size(ty,2),30);
slopefinaltycru = zeros(size(ty,2),29);
critthresholdfinaltycru = zeros(size(ty,2),1);
critavgcctycru = zeros(size(ty,2),1);
globalefficiencytycru = zeros(1,size(ty,2));
FNNtycru = zeros(size(ty,2),1);

 

%Normalization of data

 

for o =1:size(ty,2)
x = ty(:,o);
minimum = min(x);
maximum = max(x);
y2 = zeros(size(x,1),1);
for j = 1:size(x,1)
        y2(j,1) = (x(j,1)- minimum(1,1))/(maximum(1,1)- minimum(1,1));
end
x=y2;

 

% FNN determination
mmax = 20; %mmax : maximum embedding dimension
rtol = 15; %rtol: distance tolerance
atol = 2; %atol: lonliness tolerance
N=length(x); %N: length of data
Ra=std(x,1); %Ra: tandard deviation of data
xlabel('Lag');
ylabel('Autocorrelation');
title('Auto Correlation Function');
tau = 1;
for m=1:mmax
    M=N-m*tau;
    Y=zeros(M,m);
for i=1:m
    Y(:,i)=x((1:M)+(i-1)*tau)';
end
    FNN(m,1)=0;
    for n=1:M
        y0=ones(M,1)*Y(n,:);
        distance=sqrt(sum((Y-y0).^2,2));
        [neardis ,nearpos]=sort(distance);
        
        D=abs(x(n+m*tau)-x(nearpos(2)+m*tau));
        R=sqrt(D.^2+neardis(2).^2);
        if D/neardis(2) > rtol || R/Ra > atol
             FNN(m,1)=FNN(m,1)+1;
        end
    end
end
FNN=(FNN./FNN(1,1))*100;
%%
c = min(FNN);
[row,column] = find(FNN==c);
ed = 2;
FNNtycru(o,1)= ed;   % FNN value for each station
% Phase Space Reconstruction

 

PM = size(x,1)-((1+(ed-1)*tau)-1);
PY = zeros(PM,ed);
for i = 1:ed
    PY(:,i)=x((1:PM)+(i-1)*tau)';           % Reconstructed Vector
    
end

 

% Determination of distance between the reconstructed vectors

 

distance = zeros(size(PY,1),size(PY,1));
for i =1:size(PY,1)
R =[];
R(1,:)= PY(i,:);
for j = 1:size(PY,1)
distance (j,i)= sum(abs((R(1,:)-PY(j,:))));
end
end
s = size(distance);
index = 1:s(1)+1:s(1)*s(2);
distance(index) = nan;

 

% Threshold (consider 30 distance threshold between min and maximum value
% of distance

 

mindis = min(min(distance));
maxdis = max(max(distance));
thresholdtycru = (mindis:(maxdis-mindis)/30:maxdis);     % 30 threshold values

 

% Determination of critical threshold based on slope (for that for each
% threshold value we have to determine the network density

 

thresholdtycru = thresholdtycru(1,2:31);
thresholdoldtycru(o,:) = thresholdtycru';
links = zeros(size(thresholdtycru,2),size(distance,1));
actualconnection = zeros(size(thresholdtycru,2),size(distance,1));
r = zeros(size(distance,1),size(distance,1));
for kj = [1:size(thresholdtycru,2);thresholdtycru(1,1:size(thresholdtycru,2)) ]  %tricky tycrue of for!
  k = kj(1);
  j = kj(2);
for i = 1:size(distance,2)
r(:,i)= distance(:,i)< j;
links(k,i)= sum(r(:,i));
actualconnection(k,i) = sum(r(i:size(distance,1),i));
end
end
potentialconnections = (size(distance,1)*(size(distance,1)-1))/2;
networkdensity = zeros(size(thresholdtycru,2),1);
for i =1:size(thresholdtycru,2)
networkdensity(i,1) = sum(actualconnection(i,:))/potentialconnections;    % equation of network density
end

 

% determination of slope based on network density

 

slopetycru1 = zeros(size(thresholdtycru,2)-1,1);
for i = 1:size(thresholdtycru,2)-1
slopetycru1(i,1) = (networkdensity(i+1,1)-networkdensity(i,1))/(thresholdtycru(1,i+1)-thresholdtycru(1,i));
end
x = thresholdtycru(1,2:size(thresholdtycru,2));
crittycru = zeros(1,1);
maximum = max(slopetycru1(:,1));
[row]= find(slopetycru1(:,1)==maximum);
row = row(1,1);
crittycru(1,1)= x(1,row);
slopefinaltycru(o,:) = slopetycru1';
critthresholdfinaltycru(o,1) = crittycru;  % threshold corresponding to maximum slope is identified as the critical threshold
thresholdold = crittycru;

 

% clustering coefficient calculation based on critical threshold

 

column = thresholdtycru == crittycru;
links = links(column,:);
actualcon = zeros(size(thresholdold,2),size(distance,1));
r = zeros(size(distance,1),size(distance,1));
for ip = [1:size(thresholdold,2);thresholdold(1,1:size(thresholdold,2))]
     i = ip(1);
     j = ip(2);
row = zeros(max(links(i,:)),size(links,2));
for a = 1:size(distance,1)
    r(:,a)= distance(:,a)< j;
if length(find(r(:,a)==1))< max(links(i,:))
    row(1:length(find(r(:,a)==1)),a)= find(r(:,a)==1);
    row(length(find(r(:,a)==1))+1:max(links(i,:)),a)=0;
else
   [row(:,a)]= find(r(:,a)==1);
end
end
s = zeros(size(row,2),size(row,1));
for b = 1:size(row,2)
    for k =1:links(i,b)-1
      s(b,k)= sum(r(row((k+1):links(i,b),b),row(k,b)));
    end
end
for n =1:size(s,1)
    actualcon(i,n) = sum(s(n,:));
end
end
cctycru = zeros(size(thresholdold,2),size(distance,1));
for j = 1:size(thresholdold,2)
for i =1:size(distance,1)
 if links(j,i)==1 || links(j,i)==0
  cctycru(j,i)=0;
 else
     cctycru(j,i)= 2*actualcon(j,i)/(links(j,i)*(links(j,i)-1));   % clustering coefficient for each node in the network
 end
end
end
avgcctycru = zeros(size(thresholdold,2),1);
for i = 1:size(thresholdold,2)
avgcctycru(i,1) = mean(cctycru(i,:)); 
end
critavgcctycru(o,1)=avgcctycru; % avgclustering coefficient for all stations
plot(thresholdtycru(1,2:30),slopetycru1);

 

% steps for shortest pathlength

 

linksall = zeros(size(distance,1),size(distance,1));
for i = 1:size(distance,2)
linksall(:,i)= distance(:,i)< crittycru;
end
kl = zeros(1,size(linksall,2));
for i =1:size(linksall,2)
 kl(1,i) = sum(linksall(i:size(linksall,2),i));
end
s = zeros(1,sum(kl));
p = cumsum(kl);
h = p+1;
d = [1,h(1,1:size(h,2)-1)];
for ijk =[1:1:size(linksall,2);d;p]
    i = ijk(1);
    j =ijk(2);
    k = ijk(3);
s(1,j:k) = repmat(i,kl(1,i),1);
end
t = zeros(1,sum(kl));
for ijk  =[1:1:size(linksall,2);d;p]
    i= ijk(1);
    j =ijk(2);
    k =ijk(3);
    row = [];
    [row,~] = find(linksall(i:size(linksall,2),i)==1);
    roww =[];
    roww = row+i-1;
    t(1,j:k)= roww;
end

 

% graph formation

 

G = graph(s,t);
plot(G);
dk = distances(G);
hj = zeros(size(dk,1),size(dk,1));
for i = 1:size(dk,1)
    for j =1:size(dk,1)
        hj(j,i) = 1/dk(j,i);
    end
end
hj(~isfinite(hj))=0;
%dk(~isfinite(dk))=0;
shortestlength = zeros(size(dk,2),1);
for  m = 1:size(dk,2)
shortestlength(m,1) = sum((hj(m:size(dk,2),m)));
end
globalefficiencytycru(1,o) = sum(shortestlength)/(size(dk,2)*(size(dk,2)-1));
end
harmonicmeantycru = zeros(1,size(x,2));
for i =1:size(globalefficiencytycru,2)
   harmonicmeantycru(1,i) = 1/globalefficiencytycru(1,i);   % shortest path length for each station
end
save('critthresholdfinaltycru','critthresholdfinaltycru');
save('critavgcctycru','critavgcctycru');
save('slopefinaltycru','slopefinaltycru');
save('thresholdoldtycru','thresholdoldtycru');
save('globalefficiencytycru','globalefficiencytycru');
save('harmonicmeantycru','harmonicmeantycru');
save('FNNtycru','FNNtycru');
%save('degdistributiontycru','degdistributiontycru');

% to get the result for m =2 change the line 61 to ed =2;
% load ('onestation');
% ty = table2array(onestation);
% [thresholdoldtycru,cctycru,FNNtycru,slopefinaltycru,critthresholdfinaltycru,critavgcctycru,globalefficiencytycru,harmonicmeantycru]=temporalresolution_new(ty);
% plot(PY(:,1),PY(:,2));
%plot(cctycru,'o');