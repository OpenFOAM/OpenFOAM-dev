const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

const label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));

const label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));

const bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
const bool alphaApplyPrevCorr
(
    alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false)
);
