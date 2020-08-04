argList::addBoolOption
(
    "listSwitches",
    "List all available debug, info and optimisation switches"
);

#ifdef fvPatchField_H
argList::addBoolOption
(
    "listScalarBCs",
    "List scalar field boundary conditions (fvPatchField<scalar>)"
);
argList::addBoolOption
(
    "listVectorBCs",
    "List vector field boundary conditions (fvPatchField<vector>)"
);
#endif

#ifdef functionObject_H
argList::addBoolOption
(
    "listFunctionObjects",
    "List functionObjects"
);
#endif

#ifdef fvOption_H
argList::addBoolOption
(
    "listFvOptions",
    "List fvOptions"
);
#endif

#if defined(kinematicMomentumTransportModel_H) || defined(fluidThermoMomentumTransportModel_H)
argList::addBoolOption
(
    "listMomentumTransportModels",
    "List momentumTransportModels"
);
#endif
