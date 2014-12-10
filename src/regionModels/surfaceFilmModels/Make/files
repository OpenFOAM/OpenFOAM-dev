/* Surface film models */
surfaceFilmModel/surfaceFilmModel.C
surfaceFilmModel/surfaceFilmModelNew.C
noFilm/noFilm.C
kinematicSingleLayer/kinematicSingleLayer.C
thermoSingleLayer/thermoSingleLayer.C


/* Sub-models */
submodels/filmSubModelBase.C

KINEMATICMODELS=submodels/kinematic
$(KINEMATICMODELS)/force/force/force.C
$(KINEMATICMODELS)/force/force/forceNew.C
$(KINEMATICMODELS)/force/forceList/forceList.C
$(KINEMATICMODELS)/force/contactAngleForce/contactAngleForce.C
$(KINEMATICMODELS)/force/thermocapillaryForce/thermocapillaryForce.C

$(KINEMATICMODELS)/injectionModel/injectionModel/injectionModel.C
$(KINEMATICMODELS)/injectionModel/injectionModel/injectionModelNew.C
$(KINEMATICMODELS)/injectionModel/injectionModelList/injectionModelList.C
$(KINEMATICMODELS)/injectionModel/drippingInjection/drippingInjection.C
$(KINEMATICMODELS)/injectionModel/removeInjection/removeInjection.C
$(KINEMATICMODELS)/injectionModel/curvatureSeparation/curvatureSeparation.C

$(KINEMATICMODELS)/filmThermoModel/filmThermoModel/filmThermoModel.C
$(KINEMATICMODELS)/filmThermoModel/filmThermoModel/filmThermoModelNew.C
$(KINEMATICMODELS)/filmThermoModel/constantFilmThermo/constantFilmThermo.C
$(KINEMATICMODELS)/filmThermoModel/liquidFilmThermo/liquidFilmThermo.C

$(KINEMATICMODELS)/filmTurbulenceModel/filmTurbulenceModel/filmTurbulenceModel.C
$(KINEMATICMODELS)/filmTurbulenceModel/filmTurbulenceModel/filmTurbulenceModelNew.C
$(KINEMATICMODELS)/filmTurbulenceModel/laminar/laminar.C

THERMOMODELS=submodels/thermo
$(THERMOMODELS)/phaseChangeModel/phaseChangeModel/phaseChangeModel.C
$(THERMOMODELS)/phaseChangeModel/phaseChangeModel/phaseChangeModelNew.C
$(THERMOMODELS)/phaseChangeModel/noPhaseChange/noPhaseChange.C
$(THERMOMODELS)/phaseChangeModel/standardPhaseChange/standardPhaseChange.C
$(THERMOMODELS)/phaseChangeModel/solidification/solidification.C

$(THERMOMODELS)/heatTransferModel/heatTransferModel/heatTransferModel.C
$(THERMOMODELS)/heatTransferModel/heatTransferModel/heatTransferModelNew.C
$(THERMOMODELS)/heatTransferModel/constantHeatTransfer/constantHeatTransfer.C
$(THERMOMODELS)/heatTransferModel/mappedConvectiveHeatTransfer/mappedConvectiveHeatTransfer.C

$(THERMOMODELS)/filmRadiationModel/filmRadiationModel/filmRadiationModel.C
$(THERMOMODELS)/filmRadiationModel/filmRadiationModel/filmRadiationModelNew.C
$(THERMOMODELS)/filmRadiationModel/noRadiation/noRadiation.C
$(THERMOMODELS)/filmRadiationModel/constantRadiation/constantRadiation.C
$(THERMOMODELS)/filmRadiationModel/primaryRadiation/primaryRadiation.C
$(THERMOMODELS)/filmRadiationModel/standardRadiation/standardRadiation.C

$(THERMOMODELS)/filmViscosityModel/filmViscosityModel/filmViscosityModel.C
$(THERMOMODELS)/filmViscosityModel/filmViscosityModel/filmViscosityModelNew.C
$(THERMOMODELS)/filmViscosityModel/constantViscosity/constantViscosity.C
$(THERMOMODELS)/filmViscosityModel/liquidViscosity/liquidViscosity.C
$(THERMOMODELS)/filmViscosityModel/thixotropicViscosity/thixotropicViscosity.C


/* Boundary conditions */
PATCHFIELDS=derivedFvPatchFields
$(PATCHFIELDS)/filmHeightInletVelocity/filmHeightInletVelocityFvPatchVectorField.C
$(PATCHFIELDS)/inclinedFilmNusseltHeight/inclinedFilmNusseltHeightFvPatchScalarField.C
$(PATCHFIELDS)/inclinedFilmNusseltInletVelocity/inclinedFilmNusseltInletVelocityFvPatchVectorField.C


LIB = $(FOAM_LIBBIN)/libsurfaceFilmModels
