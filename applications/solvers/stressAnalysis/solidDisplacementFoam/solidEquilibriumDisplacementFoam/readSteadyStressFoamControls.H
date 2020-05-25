const scalar accFac
(
    mesh.solutionDict().subDict("stressAnalysis")
   .lookup<scalar>("accelerationFactor")
);
