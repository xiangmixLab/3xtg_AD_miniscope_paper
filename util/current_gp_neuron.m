function neuron=current_gp_neuron(neuron_ori,gp,gp_idx)

neuron=neuron_ori;
neuron.C=neuron.C(gp==gp_idx,:);
neuron.S=neuron.S(gp==gp_idx,:);
neuron.C_raw=neuron.C_raw(gp==gp_idx,:);
neuron.A=neuron.A(:,gp==gp_idx);
neuron.centroid=neuron.centroid(gp==gp_idx);


