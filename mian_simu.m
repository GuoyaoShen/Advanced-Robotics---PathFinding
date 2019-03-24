% This script is used as main script to start the simulation.

% map = load_map(filename, res_xy, res_z, margin)
map = load_map('map1.txt', 0.1, 0.8, 0.5);
% function [path, num_expanded] = dijkstra(map, start, goal, astar)
start = [5.0, -5.0, 4.0];
goal = [8.0 ,20.0 ,6.0];
% [path_simu, num_expanded_simu] = dijkstra_new(map, start, goal);
% [path_simu, num_expanded_simu] = dijkstra(map, start, goal);
[path_simu, num_expanded_simu] = dijkstra_array(map, start, goal, true);
% function plot_path(map, path)
plot_path(map, path_simu);