function [I_distances, L_ratios] = sort_quality_wenyi(file_name)
% calculates Lratio and Idist from the neuralynx output.
%
% file_name is the path to the spike (.ntt) file, and should be a string
% e.g., 'C:\Users\ampm1\Desktop\TT2_sorted.ntt'

%
% Calculations are based on Schmitzer-Torbert, N., Jackson, J., Henze, D., 
% Harris, K., and Redish, A.D.(2005). Quantitative measures of cluster 
% quality for use in extracellular recordings. Neuroscience 131, 1–11.
%


%load waveform features
[TimestampsTT, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(file_name, [1 1 1 1 1], 1, 1, [] );

%preallocate
I_distances = nan(length(unique(CellNumbers(CellNumbers>0))), 1);
L_ratios = nan(size(I_distances));

%calculate sort quality measures for each neuron
cell_ct = 0;
for ic = 1:unique(CellNumbers(CellNumbers>0))

    %new cell
    cell_ct = cell_ct+1;
    
    %isolate waveforms of interest
    cluster_locs = Features(CellNumbers==clusts(ic),:)';
    noise_locs = Features(CellNumbers~=clusts(ic),:)';

    %quality measures
    L_rat = Lratio(cluster_locs, noise_locs);
    I_dist = idist(cluster_locs, noise_locs);

    %load output
    L_ratios(cell_ct) = L_rat;
    I_distances(cell_ct) = I_dist;

end

%I distance function
function id_out = idist(cluster_locs, noise_locs)
%calculates isolation distance, a measure of how far sorted spikes are from
%noise (and other spikes)
%
%Isolation distance is the radius of the smallest ellipsoid form the
%cluster center that contains all of the cluster spikes and an equal number
%of noise spikes. Larger distances indicate greater segregation.
%
%
% From Schmitzer-Torbert, Jackson, Henze, Harris & Redish, 2005:
% "Isolation Distance was first introduced by Harris et al. (2001)
% and applied to hippocampal data sets. If a cluster contains nC
% cluster spikes, the Isolation Distance of the cluster is the D2
% value of the nCth closest noise spike (see Fig. 1). Isolation
% Distance is therefore the radius of the smallest ellipsoid from
% the cluster center containing all of the cluster spikes and an
% equal number of noise spikes. As such, Isolation Distance
% estimates how distant the cluster spikes are from the other
% spikes recorded on the same electrode. Isolation Distance is
% not defined for cases in which the number of cluster spikes is
% greater than the number of noise spikes."
%

%if fewer noise spikes than cluster spikes
if size(noise_locs, 1) < size(cluster_locs, 1)
    %downsample cluster locs
    dwn_samp_idx = randperm( size(cluster_locs, 1));
    dwn_samp_idx = dwn_samp_idx(1:size(noise_locs, 1));
    cluster_locs = cluster_locs(dwn_samp_idx, :);
end

%find number of cluster spikes
num_clust_spks = size(cluster_locs,1);

%find distances from each noise point to the cluster loc cloud
noise_dists = mahal(noise_locs, cluster_locs);

%sort those distances
noise_dists = sort(noise_dists);

%find the nth smallest distance, where n is the number of cluster spikes.
id_out = noise_dists(num_clust_spks);

%L ratio function
function lratio_out = Lratio(cluster_locs, noise_locs)
%calculates L ratio, a measure of how close noise spikes are to the center 
% of the sorted cluster compared to sorted spikes.
%

%calculate all distances from each cluster_loc to cluster loc cloud
cluster_distances = mahal(cluster_locs, cluster_locs);

%and for each noise loc to the cluster loc cloud
noise_distances = mahal(noise_locs, cluster_locs);

%calculate the the proportion of cluster distance smaller than each noise
%distance
L_vals = nan(length(noise_distances),1);
for ind = 1:length(noise_distances)
    L_vals(ind) = 1 - sum(cluster_distances<noise_distances(ind))...
        /length(cluster_distances);
end

%divide by the number of spikes in the cluster
lratio_out = sum(L_vals)/length(cluster_locs);