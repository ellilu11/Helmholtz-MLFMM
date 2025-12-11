clear; clc;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";
% srcs = readmatrix(dir+"config\n120\vertices.txt");
nodes = readmatrix(dir+"out\nodes.txt");

nodePos = nodes(:,1:3);
nodeLengs = nodes(:,4);
nodeVec = [nodes(:,1:4),nodeLengs,nodeLengs];
nodeVec(:,1:3) = nodeVec(:,1:3) - nodeVec(:,4:6)/2;

%% Plot nodes
rootLeng = 11.0;
lim = [-rootLeng/2 rootLeng/2];

figure(1)
scatter3(nodes(:,1),nodes(:,2),nodes(:,3),stat2rgb(5),'filled');

%%
function rgb = stat2rgb(stat)
    switch stat
        case 0
            rgb = "none";
        case 1
            rgb = "black"; % self
        case 2
            rgb = "green"; % list 1
        case 3 
            rgb = "yellow";  % list 2
        case 4 
            rgb = "cyan"; % list 3 (leaf)
        case 5 
            rgb = "blue";  % list 3 (stem)
        case 6 
            rgb = "magenta"; % list 4
        case 7 
            rgb = "red"; 
        case 8 
            rgb = "yellow"; 
    end
end