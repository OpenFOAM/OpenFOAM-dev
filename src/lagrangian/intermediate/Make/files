PARCELS=parcels
BASEPARCELS=$(PARCELS)/baseClasses
DERIVEDPARCELS=$(PARCELS)/derived

CLOUDS=clouds
BASECLOUDS=$(CLOUDS)/baseClasses
DERIVEDCLOUDS=$(CLOUDS)/derived


/* Cloud base classes */
$(BASECLOUDS)/kinematicCloud/kinematicCloud.C
$(BASECLOUDS)/thermoCloud/thermoCloud.C
$(BASECLOUDS)/reactingCloud/reactingCloud.C
$(BASECLOUDS)/reactingMultiphaseCloud/reactingMultiphaseCloud.C


/* kinematic parcel sub-models */
KINEMATICPARCEL=$(DERIVEDPARCELS)/basicKinematicParcel
$(KINEMATICPARCEL)/defineBasicKinematicParcel.C
$(KINEMATICPARCEL)/makeBasicKinematicParcelSubmodels.C


/* kinematic colliding parcel sub-models */
KINEMATICCOLLIDINGPARCEL=$(DERIVEDPARCELS)/basicKinematicCollidingParcel
$(KINEMATICCOLLIDINGPARCEL)/defineBasicKinematicCollidingParcel.C
$(KINEMATICCOLLIDINGPARCEL)/makeBasicKinematicCollidingParcelSubmodels.C


/* thermo parcel sub-models */
THERMOPARCEL=$(DERIVEDPARCELS)/basicThermoParcel
$(THERMOPARCEL)/defineBasicThermoParcel.C
$(THERMOPARCEL)/makeBasicThermoParcelSubmodels.C


/* reacting parcel sub-models */
REACTINGPARCEL=$(DERIVEDPARCELS)/basicReactingParcel
$(REACTINGPARCEL)/defineBasicReactingParcel.C
$(REACTINGPARCEL)/makeBasicReactingParcelSubmodels.C


/* reacting multiphase parcel sub-models */
REACTINGMPPARCEL=$(DERIVEDPARCELS)/basicReactingMultiphaseParcel
$(REACTINGMPPARCEL)/defineBasicReactingMultiphaseParcel.C
$(REACTINGMPPARCEL)/makeBasicReactingMultiphaseParcelSubmodels.C


/* kinematic MPPIC parcel sub-models */
KINEMATICMPPICPARCEL=$(DERIVEDPARCELS)/basicKinematicMPPICParcel
$(KINEMATICMPPICPARCEL)/defineBasicKinematicMPPICParcel.C
$(KINEMATICMPPICPARCEL)/makeBasicKinematicMPPICParcelSubmodels.C


/* bolt-on models */
RADIATION=submodels/addOns/radiation
$(RADIATION)/absorptionEmission/cloudAbsorptionEmission/cloudAbsorptionEmission.C
$(RADIATION)/scatter/cloudScatter/cloudScatter.C

submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionData.C
submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionDataList.C

KINEMATICINJECTION=submodels/Kinematic/InjectionModel
$(KINEMATICINJECTION)/KinematicLookupTableInjection/kinematicParcelInjectionData.C
$(KINEMATICINJECTION)/KinematicLookupTableInjection/kinematicParcelInjectionDataIO.C
$(KINEMATICINJECTION)/KinematicLookupTableInjection/kinematicParcelInjectionDataIOList.C
$(KINEMATICINJECTION)/PatchInjection/patchInjectionBase.C

THERMOINJECTION=submodels/Thermodynamic/InjectionModel
$(THERMOINJECTION)/ThermoLookupTableInjection/thermoParcelInjectionData.C
$(THERMOINJECTION)/ThermoLookupTableInjection/thermoParcelInjectionDataIO.C
$(THERMOINJECTION)/ThermoLookupTableInjection/thermoParcelInjectionDataIOList.C

REACTINGINJECTION=submodels/Reacting/InjectionModel
$(REACTINGINJECTION)/ReactingLookupTableInjection/reactingParcelInjectionData.C
$(REACTINGINJECTION)/ReactingLookupTableInjection/reactingParcelInjectionDataIO.C
$(REACTINGINJECTION)/ReactingLookupTableInjection/reactingParcelInjectionDataIOList.C

REACTINGMPINJECTION=submodels/ReactingMultiphase/InjectionModel
$(REACTINGMPINJECTION)/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionData.C
$(REACTINGMPINJECTION)/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIO.C
$(REACTINGMPINJECTION)/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIOList.C

MPPICPARTICLESTRESS=submodels/MPPIC/ParticleStressModels
$(MPPICPARTICLESTRESS)/ParticleStressModel/ParticleStressModel.C
$(MPPICPARTICLESTRESS)/HarrisCrighton/HarrisCrighton.C
$(MPPICPARTICLESTRESS)/Lun/Lun.C
$(MPPICPARTICLESTRESS)/exponential/exponential.C

MPPICCORRECTIONLIMITING=submodels/MPPIC/CorrectionLimitingMethods
$(MPPICCORRECTIONLIMITING)/CorrectionLimitingMethod/CorrectionLimitingMethod.C
$(MPPICCORRECTIONLIMITING)/noCorrectionLimiting/noCorrectionLimiting.C
$(MPPICCORRECTIONLIMITING)/absolute/absolute.C
$(MPPICCORRECTIONLIMITING)/relative/relative.C

MPPICTIMESCALE=submodels/MPPIC/TimeScaleModels
$(MPPICTIMESCALE)/TimeScaleModel/TimeScaleModel.C
$(MPPICTIMESCALE)/equilibrium/equilibrium.C
$(MPPICTIMESCALE)/nonEquilibrium/nonEquilibrium.C
$(MPPICTIMESCALE)/isotropic/isotropic.C


/* integration schemes */
IntegrationScheme/makeIntegrationSchemes.C


/* phase properties */
phaseProperties/phaseProperties/phaseProperties.C
phaseProperties/phaseProperties/phasePropertiesIO.C
phaseProperties/phasePropertiesList/phasePropertiesList.C


/* additional helper classes */
clouds/Templates/KinematicCloud/cloudSolution/cloudSolution.C


/* averaging methods */
submodels/MPPIC/AveragingMethods/makeAveragingMethods.C


LIB = $(FOAM_LIBBIN)/liblagrangianIntermediate
