numFiles = 331;
gridSize = 256;
outputFreq = 20;

% LAYER 1
alldata1 = zeros(gridSize,gridSize,numFiles);

for k = 0:numFiles-1
  csvFileName = sprintf('out1_%08d.csv', outputFreq*k);
  csvdata = csvread(csvFileName);
  alldata1(:,:,k+1) = csvdata(:,1:gridSize);
end

k = 1;
while k <= numFiles
  h = pcolor(alldata1(:,:,k));
  set(h, 'EdgeColor', 'none');
  pause(.1);
  k = k + 1;
  if k == numFiles
    k = 1;
  end
end

% LAYER 2
alldata2 = zeros(gridSize,gridSize,numFiles);

for k = 0:numFiles-1
  csvFileName=sprintf('out2_%08d.csv', outputFreq*k);
  csvdata = csvread(csvFileName);
  alldata2(:,:,k+1) = csvdata(:,1:gridSize);
end

k = 1;
while k <= numFiles
  h = pcolor(alldata2(:,:,k));
  set(h, 'EdgeColor', 'none');
  pause(.1);
  k = k + 1;
  if k == numFiles
    k = 1;
  end
end