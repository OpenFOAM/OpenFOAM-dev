constant/
    cavity blockMesh
    point (0.05 0.05 0.01) moved to (0.05 0.05 0.001)

0.005/
    collapseEdges with

        collapseEdgesCoeffs.minimumEdgeLength   0.0011;

    so it collapses the one edge


processor*/0.005/
    decomposePar


mpirun -np 2 Test-patchRegion -parallel
