/* RAS turbulence models */
RASModel/RASModel.C
laminar/laminar.C
kEpsilon/kEpsilon.C
RNGkEpsilon/RNGkEpsilon.C
LaunderSharmaKE/LaunderSharmaKE.C
LRR/LRR.C
LaunderGibsonRSTM/LaunderGibsonRSTM.C
realizableKE/realizableKE.C
SpalartAllmaras/SpalartAllmaras.C
kOmegaSST/kOmegaSST.C
v2f/v2f.C

/* Wall functions */
wallFunctions = derivedFvPatchFields/wallFunctions

alphatWallFunctions = $(wallFunctions)/alphatWallFunctions
$(alphatWallFunctions)/alphatWallFunction/alphatWallFunctionFvPatchScalarField.C
$(alphatWallFunctions)/alphatJayatillekeWallFunction/alphatJayatillekeWallFunctionFvPatchScalarField.C

mutWallFunctions = $(wallFunctions)/mutWallFunctions
$(mutWallFunctions)/mutWallFunction/mutWallFunctionFvPatchScalarField.C

$(mutWallFunctions)/mutkWallFunction/mutkWallFunctionFvPatchScalarField.C
$(mutWallFunctions)/mutkRoughWallFunction/mutkRoughWallFunctionFvPatchScalarField.C

$(mutWallFunctions)/mutUWallFunction/mutUWallFunctionFvPatchScalarField.C
$(mutWallFunctions)/mutURoughWallFunction/mutURoughWallFunctionFvPatchScalarField.C
$(mutWallFunctions)/mutUSpaldingWallFunction/mutUSpaldingWallFunctionFvPatchScalarField.C
$(mutWallFunctions)/mutLowReWallFunction/mutLowReWallFunctionFvPatchScalarField.C

epsilonWallFunctions = $(wallFunctions)/epsilonWallFunctions
$(epsilonWallFunctions)/epsilonWallFunction/epsilonWallFunctionFvPatchScalarField.C
$(epsilonWallFunctions)/epsilonLowReWallFunction/epsilonLowReWallFunctionFvPatchScalarField.C

fWallFunctions = $(wallFunctions)/fWallFunctions
$(fWallFunctions)/fWallFunction/fWallFunctionFvPatchScalarField.C

omegaWallFunctions = $(wallFunctions)/omegaWallFunctions
$(omegaWallFunctions)/omegaWallFunction/omegaWallFunctionFvPatchScalarField.C

kqRWallFunctions = $(wallFunctions)/kqRWallFunctions
$(kqRWallFunctions)/kqRWallFunction/kqRWallFunctionFvPatchFields.C
$(kqRWallFunctions)/kLowReWallFunction/kLowReWallFunctionFvPatchScalarField.C

v2WallFunctions = $(wallFunctions)/v2WallFunctions
$(v2WallFunctions)/v2WallFunction/v2WallFunctionFvPatchScalarField.C

/* Patch fields */
derivedFvPatchFields/convectiveHeatTransfer/convectiveHeatTransferFvPatchScalarField.C
derivedFvPatchFields/turbulentMixingLengthDissipationRateInlet/turbulentMixingLengthDissipationRateInletFvPatchScalarField.C
derivedFvPatchFields/turbulentMixingLengthFrequencyInlet/turbulentMixingLengthFrequencyInletFvPatchScalarField.C
backwardsCompatibility/wallFunctions/backwardsCompatibilityWallFunctions.C

LIB = $(FOAM_LIBBIN)/libcompressibleRASModels
