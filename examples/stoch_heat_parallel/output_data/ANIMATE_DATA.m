function ANIMATE_DATA(dim1_len, dim2_len, file_count, output_freq)
% Creates an animation of the output data.
%   ARGUMENTS: - dim1_len, dim2_len: The dimensions of the grid in the
%   Fortran code. Note that MATLAB is row-major and not column major, so
%   dim1_len is the length along the second index and vice-versa in MATLAB.
%   - num_files: The number of output files.
%   - output_freq: The output frequency of the simulation.

% We will store all of the output in a single array, and call them for when
% we plot them.
all_data = zeros(dim2_len, dim1_len, file_count);
% Fill in all_data.
for file_num = 0:file_count-1
    file_name = sprintf('out_%08d.csv', output_freq*file_num);
    file_data = readmatrix(file_name);
    all_data(:,:,file_num+1) = file_data(:,:);
end

% Animate the plot by rewriting it repeatedly.
file_num = 1;
while file_num <= file_count
    h = pcolor(all_data(:,:,file_num));
    set(h, 'EdgeColor', 'none');
    colorbar;
    caxis([60 70]);
    title(int2str(file_num-1));
    pause(.1);
    file_num = file_num + 1;
    if file_num == file_count
        file_num = 1;
    end
end

end

