To run this code, open the script `poisson.m` in MATLAB and press the RUN button. Change the `hfun` parameter to control the mesh resolution. The mesh processing can take some time due to the neighbor-finding step, and unfortunately the current implementation does not scale well. Expect to wait a while if you set `hfun` less than 0.025 :)


This unstructured mesh generation in this code is heavily dependent on the [MESH2d library](https://github.com/dengwirda/mesh2d) by Darren Engwirda.
