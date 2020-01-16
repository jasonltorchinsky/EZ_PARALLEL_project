function cfft2df_unit_test_MATLAB
% Recreation of the outputs of paralel_fft_unit_test in EZ_PARALLEL
% library.

dim2_len = 10; % dim1_len in Fortran.
dim1_len = 8; % dim2_len in Fortran.
dim2_ref = 0.0; % dim1_ref in Fortran.
dim1_ref = 0.0; % dim2_ref in Fortran.
dim2_spc = 0.25; % dim1_spc in Fortran.
dim1_spc = 0.25; % dim2_spc in Fortran.

grid = ones(dim1_len, dim2_len);

for j = 1:dim1_len
    for i = 1:dim2_len
        grid(j,i) = sin((dim1_ref + (i-1) * dim1_spc) * pi);
    end
end

disp(grid);

grid = fft2(grid);

disp(grid(:,:));

end

