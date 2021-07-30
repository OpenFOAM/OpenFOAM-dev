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

#ifdef fvModel_H
if (args.optionFound("listFvModels"))
{
    Info<< "fvModels"
        << Foam::fvModel::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef fvConstraint_H
if (args.optionFound("listFvConstraints"))
{
    Info<< "fvConstraints"
        << Foam::fvConstraint::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef incompressibleMomentumTransportModels_H
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
#elif defined(compressibleMomentumTransportModels_H)
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
