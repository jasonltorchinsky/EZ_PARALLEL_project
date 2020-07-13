function ANIMATE_DATA(numPts, numSteps, outputFreq, dt)
% Creates an animation of the output data.
%   ARGUMENTS: - numPts: Number of lattice points in one direction.
%   - numSteps: Number of time-steps executed.
%   - outputFreq: The output frequency of the simulation.
%   - dt: Time-step size

% We will store all of the output in a single array, and call them for when
% we plot them.
fileCount = floor(numSteps/outputFreq)+1;
all_data = zeros(numPts, numPts, fileCount);
% Fill in all_data.
for k = 0:fileCount-1
    csvFileName=sprintf('out_%08d.csv', outputFreq*k);
    csvdata = csvread(csvFileName);
    all_data(:,:,k+1) = csvdata(:,1:numPts);
end

% Animate the plot by rewriting it repeatedly.
k = 1;
while k <= fileCount
    h = pcolor(all_data(:,:,k));
    set(h, 'EdgeColor', 'none');
    %colorbar;
    caxis([65.99 66.01]);
    colormap(gray);
    label = ["Column-Integrated Water Vapor: ", 
        num2str(outputFreq*(k-1)*dt,'%03.3f')];
    title(join(label));
    pause(.1);
    saveas(h,sprintf('out_%08d.png', outputFreq*(k-1)));
    k = k + 1;
    %if k == fileCount
    %    k = 1;
    %end
end

end

