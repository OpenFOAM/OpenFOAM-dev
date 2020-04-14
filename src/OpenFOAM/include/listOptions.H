argList::addBoolOption
(
    "listSwitches",
    "List switches declared in libraries but not set in etc/controlDict"
);
argList::addBoolOption
(
    "listRegisteredSwitches",
    "List switches registered for run-time modification"
);
argList::addBoolOption
(
    "listUnsetSwitches",
    "List switches declared in libraries but not set in etc/controlDict"
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
