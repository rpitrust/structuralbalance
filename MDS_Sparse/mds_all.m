function [] = mds_all(data_name, dim)

data=data_name(1:size(data_name,2)-4);
%path = strcat(pwd,'/', data, '_expand');
path = strcat(pwd,'/', data, 'D');
mkdir(path)
for i=0:9,
    %file_name=strcat('cut_expand_re_train', int2str(i), '_', data_name);
    file_name=strcat('re_train', int2str(i), '_', data_name);
    mds_sparse(file_name, dim, path);
    %mds_weight(file_name, dim, path);
end


