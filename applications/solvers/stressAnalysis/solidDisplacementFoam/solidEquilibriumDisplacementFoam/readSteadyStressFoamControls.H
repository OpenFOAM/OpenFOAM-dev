const scalar accFac
(
    mesh.solution().dict().subDict("stressAnalysis")
   .lookup<scalar>("accelerationFactor")
);
