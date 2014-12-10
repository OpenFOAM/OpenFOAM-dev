const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
bool alphaApplyPrevCorr
(
    alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false)
);

// Isotropic compression coefficient
scalar icAlpha
(
    alphaControls.lookupOrDefault<scalar>("icAlpha", 0)
);
