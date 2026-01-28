[z1, z2] = meshgrid(-1 : 0.01 : 1, -1 : 0.01 : 1);
h = hamiltonian(z1, z2);
surf(z1, z2, h);
