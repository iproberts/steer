function [nbr,nbr_az,nbr_el] = construct_neighborhood(Delta_az,Delta_el,delta_az,delta_el)

% azimuthal neighborhood
val = floor(Delta_az / delta_az);
nbr_az = (-val:1:val).' .* delta_az;

% elevational neighborhood
val = floor(Delta_el / delta_el);
nbr_el = (-val:1:val).' .* delta_el;

% full az-el neighborhood
nbr = [repmat(nbr_az,length(nbr_el),1) repelem(nbr_el,length(nbr_az),1)];

% sort in ascending distance
[~,idx_sort] = sort(sum(nbr.^2,2),'ascend');
nbr = nbr(idx_sort,:);

end