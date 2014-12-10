PARCELS=parcels
DERIVEDPARCELS=$(PARCELS)/derived

KINEMATICPARCEL=$(DERIVEDPARCELS)/basicKinematicParcel
$(KINEMATICPARCEL)/makeBasicKinematicParcelSubmodels.C

KINEMATICCOLLIDINGPARCEL=$(DERIVEDPARCELS)/basicKinematicCollidingParcel
$(KINEMATICCOLLIDINGPARCEL)/makeBasicKinematicCollidingParcelSubmodels.C

THERMOPARCEL=$(DERIVEDPARCELS)/basicThermoParcel
$(THERMOPARCEL)/makeBasicThermoParcelSubmodels.C

REACTINGPARCEL=$(DERIVEDPARCELS)/basicReactingParcel
$(REACTINGPARCEL)/makeBasicReactingParcelSubmodels.C

REACTINGMPPARCEL=$(DERIVEDPARCELS)/basicReactingMultiphaseParcel
$(REACTINGMPPARCEL)/makeBasicReactingMultiphaseParcelSubmodels.C

LIB = $(FOAM_LIBBIN)/liblagrangianTurbulence
