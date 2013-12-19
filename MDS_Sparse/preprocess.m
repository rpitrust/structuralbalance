function preprocess(file_name)

graph_data=dlmread(file_name, '\t');
graph_data(:,1:2)=graph_data(:,1:2)+1;
graph_data(:,3)=graph_data(:,3)+2;
graph_data=graph_data(:,1:3);
graph_data=[graph_data;[graph_data(:,2),graph_data(:,1),graph_data(:,3)]];
name = strcat('MatlabReady_', file_name);
dlmwrite(name, graph_data, 'delimiter', '\t');