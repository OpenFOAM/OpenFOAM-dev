RASModel/RASModel.C

laminar/laminar.C
kEpsilon/kEpsilon.C
RNGkEpsilon/RNGkEpsilon.C
realizableKE/realizableKE.C
kOmega/kOmega.C
kOmegaSST/kOmegaSST.C
SpalartAllmaras/SpalartAllmaras.C
LRR/LRR.C
LaunderGibsonRSTM/LaunderGibsonRSTM.C
LaunderSharmaKE/LaunderSharmaKE.C
qZeta/qZeta.C
LienCubicKE/LienCubicKE.C
LienCubicKELowRe/LienCubicKELowRe.C
NonlinearKEShih/NonlinearKEShih.C
LienLeschzinerLowRe/LienLeschzinerLowRe.C
LamBremhorstKE/LamBremhorstKE.C
kkLOmega/kkLOmega.C
v2f/v2f.C

/* Wall functions */
wallFunctions = derivedFvPatchFields/wallFunctions

nutWallFunctions = $(wallFunctions)/nutWallFunctions
$(nutWallFunctions)/nutWallFunction/nutWallFunctionFvPatchScalarField.C

$(nutWallFunctions)/nutkWallFunction/nutkWallFunctionFvPatchScalarField.C
$(nutWallFunctions)/nutkRoughWallFunction/nutkRoughWallFunctionFvPatchScalarField.C
$(nutWallFunctions)/nutkAtmRoughWallFunction/nutkAtmRoughWallFunctionFvPatchScalarField.C

$(nutWallFunctions)/nutUWallFunction/nutUWallFunctionFvPatchScalarField.C
$(nutWallFunctions)/nutUSpaldingWallFunction/nutUSpaldingWallFunctionFvPatchScalarField.C
$(nutWallFunctions)/nutUTabulatedWallFunction/nutUTabulatedWallFunctionFvPatchScalarField.C
$(nutWallFunctions)/nutLowReWallFunction/nutLowReWallFunctionFvPatchScalarField.C
$(nutWallFunctions)/nutURoughWallFunction/nutURoughWallFunctionFvPatchScalarField.C

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

alphatWallFunctions = $(wallFunctions)/alphatWallFunctions
$(alphatWallFunctions)/alphatJayatillekeWallFunction/alphatJayatillekeWallFunctionFvPatchScalarField.C

v2WallFunctions = $(wallFunctions)/v2WallFunctions
$(v2WallFunctions)/v2WallFunction/v2WallFunctionFvPatchScalarField.C

/* Patch fields */
derivedFvPatchFields/turbulentMixingLengthDissipationRateInlet/turbulentMixingLengthDissipationRateInletFvPatchScalarField.C
derivedFvPatchFields/turbulentMixingLengthFrequencyInlet/turbulentMixingLengthFrequencyInletFvPatchScalarField.C
derivedFvPatchFields/atmBoundaryLayerInletEpsilon/atmBoundaryLayerInletEpsilonFvPatchScalarField.C

backwardsCompatibility/wallFunctions/backwardsCompatibilityWallFunctions.C

LIB = $(FOAM_LIBBIN)/libincompressibleRASModels
