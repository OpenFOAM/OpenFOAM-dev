bool listOptions = false;

if (args.optionFound("listSwitches"))
{
    debug::listSwitches();
    listOptions = true;
}

#ifdef fvPatchField_H
if (args.optionFound("listScalarBCs"))
{
    Info<< "scalarBCs"
        << fvPatchField<scalar>::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}

if (args.optionFound("listVectorBCs"))
{
    Info<< "vectorBCs"
        << fvPatchField<vector>::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef functionObject_H
if (args.optionFound("listFunctionObjects"))
{
    Info<< "functionObjects"
        << functionObject::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef fvOption_H
if (args.optionFound("listFvOptions"))
{
    Info<< "fvOptions"
        << fv::option::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef kinematicMomentumTransportModel_H
if (args.optionFound("listMomentumTransportModels"))
{
    Info<< "Turbulence models"
        << incompressible::momentumTransportModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "RAS models"
        << incompressible::RASModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "LES models"
        << incompressible::LESModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#elif defined(fluidThermoMomentumTransportModel_H)
if (args.optionFound("listMomentumTransportModels"))
{
    Info<< "Turbulence models"
        << compressible::momentumTransportModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "RAS models"
        << compressible::RASModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "LES models"
        << compressible::LESModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

if (listOptions)
{
    exit(0);
}
