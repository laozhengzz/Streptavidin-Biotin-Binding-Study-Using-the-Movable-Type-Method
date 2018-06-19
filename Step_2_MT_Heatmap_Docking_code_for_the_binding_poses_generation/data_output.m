%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  data_output(data_mat,output_dir,code_dir)

cd(output_dir);

save data.txt data_mat -ascii

cd(code_dir);


