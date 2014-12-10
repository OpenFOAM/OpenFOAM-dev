turbulenceModel.C

LESdelta = LES/LESdeltas

$(LESdelta)/LESdelta/LESdelta.C
$(LESdelta)/cubeRootVolDelta/cubeRootVolDelta.C
$(LESdelta)/PrandtlDelta/PrandtlDelta.C
$(LESdelta)/smoothDelta/smoothDelta.C
$(LESdelta)/maxDeltaxyz/maxDeltaxyz.C

LESfilters = LES/LESfilters

$(LESfilters)/LESfilter/LESfilter.C
$(LESfilters)/simpleFilter/simpleFilter.C
$(LESfilters)/laplaceFilter/laplaceFilter.C
$(LESfilters)/anisotropicFilter/anisotropicFilter.C


derivedFvPatchFields/fixedShearStress/fixedShearStressFvPatchVectorField.C
derivedFvPatchFields/porousBafflePressure/porousBafflePressureFvPatchField.C

/* Wall functions */
wallFunctions = RAS/derivedFvPatchFields/wallFunctions

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

omegaWallFunctions = $(wallFunctions)/omegaWallFunctions
$(omegaWallFunctions)/omegaWallFunction/omegaWallFunctionFvPatchScalarField.C

kqRWallFunctions = $(wallFunctions)/kqRWallFunctions
$(kqRWallFunctions)/kqRWallFunction/kqRWallFunctionFvPatchFields.C
$(kqRWallFunctions)/kLowReWallFunction/kLowReWallFunctionFvPatchScalarField.C

/*
v2WallFunctions = $(wallFunctions)/v2WallFunctions
$(v2WallFunctions)/v2WallFunction/v2WallFunctionFvPatchScalarField.C

fWallFunctions = $(wallFunctions)/fWallFunctions
$(fWallFunctions)/fWallFunction/fWallFunctionFvPatchScalarField.C
*/

/* Patch fields */
RAS/derivedFvPatchFields/turbulentMixingLengthDissipationRateInlet/turbulentMixingLengthDissipationRateInletFvPatchScalarField.C
RAS/derivedFvPatchFields/turbulentMixingLengthFrequencyInlet/turbulentMixingLengthFrequencyInletFvPatchScalarField.C
RAS/derivedFvPatchFields/atmBoundaryLayerInletEpsilon/atmBoundaryLayerInletEpsilonFvPatchScalarField.C
RAS/derivedFvPatchFields/atmBoundaryLayerInletVelocity/atmBoundaryLayerInletVelocityFvPatchVectorField.C

LIB = $(FOAM_LIBBIN)/libturbulenceModels
