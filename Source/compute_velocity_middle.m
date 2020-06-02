function [um vm] = compute_velocity_middle(umac,vmac)
um = .5*(umac(1:end-1,:) + umac(2:end,:));
vm = .5*(vmac(:,1:end-1) + vmac(:,2:end,:));


