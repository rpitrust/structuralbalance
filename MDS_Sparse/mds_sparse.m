function [X,hist,total_time] = mds_sparse(file_name, dim, path)

%read the file to a sparse matrix Adj
graph_data=dlmread(file_name, '\t');
graph_data(:,1:2)=graph_data(:,1:2)+1;
graph_data(:,3)=graph_data(:,3)+2;
graph_data=graph_data(:,1:3);
graph_data=[graph_data;[graph_data(:,2),graph_data(:,1),graph_data(:,3)]];
Adj= spconvert(graph_data);


fprintf('Sparse Adj size:\n');
fprintf(1,'%4d   \n', size(Adj,1)); 


wp=1;
wn=3;
wu=0.01;
dp=0.1;
dn=1.3;
du=0.6;

% initialization
X0= randn(size(Adj,1),dim);
X0(1,:)=zeros(1,dim);
% set defauls
rtol = 0.001; 
atol=0;
iter = 1000; 
verbose='iter';
tol=0.00000001;
xhistory = 'off';



% start optimization
[X,hist,total_time] = smacof_sparse(Adj,X0,iter,verbose,xhistory,rtol,atol, tol, wp, wn, wu, dp, dn, du);


name = strcat('loc_', file_name);
path = strcat(path, '/');
dlmwrite([path name], X, 'delimiter', '\t', ...  
    'precision', 6)





