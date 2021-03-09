/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "massSource.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(massSource, 0);
    addToRunTimeSelectionTable(fvModel, massSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::massSource::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    rhoName_ =
        coeffs().lookupOrDefault<word>
        (
            "rho",
            IOobject::groupName("rho", phaseName_)
        );

    if
    (
        mesh().foundObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phaseName_)
        )
    )
    {
        const basicThermo& thermo =
            mesh().lookupObject<basicThermo>
            (
                IOobject::groupName(basicThermo::dictName, phaseName_)
            );
        heName_ = thermo.he().name();
        TName_ = thermo.T().name();
    }

    fieldValues_.clear();
    const dictionary& fieldCoeffs = coeffs().subDict("fieldValues");
    forAllConstIter(dictionary, fieldCoeffs, iter)
    {
        fieldValues_.set(iter().keyword(), nullptr);
    }

    #define callReadFieldValues(Type, nullArg) readFieldValues<Type>();
    FOR_ALL_FIELD_TYPES(callReadFieldValues);
    #undef callReadFieldValues

    massFlowRate_.reset(Function1<scalar>::New("massFlowRate", coeffs()).ptr());
}


template<class Type>
void Foam::fv::massSource::readFieldValues()
{
    const dictionary& fieldCoeffs = coeffs().subDict("fieldValues");

    forAllConstIter(dictionary, fieldCoeffs, iter)
    {
        const word& fieldName = iter().keyword();

        if (mesh().foundObject<VolField<Type>>(fieldName))
        {
            fieldValues_.set
            (
                fieldName,
                objectFunction1::New<VolField>
                (
                    fieldName,
                    fieldCoeffs,
                    fieldName,
                    mesh()
                ).ptr()
            );
        }
    }
}


template<class Type>
void Foam::fv::massSource::addGeneralSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const scalar t = mesh().time().value();
    const scalar massFlowRate = massFlowRate_->value(t);
    const Type value = fieldValues_[fieldName]->value<Type>(t);

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        eqn.source()[cells[i]] -=
            mesh().V()[cells[i]]/set_.V()*massFlowRate*value;
    }
}


template<class Type>
void Foam::fv::massSource::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addGeneralSupType(eqn, fieldName);
}


void Foam::fv::massSource::addSupType
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const labelList& cells = set_.cells();

    if (fieldName == rhoName_)
    {
        const scalar t = mesh().time().value();
        const scalar massFlowRate = massFlowRate_->value(t);

        forAll(cells, i)
        {
            eqn.source()[cells[i]] -=
                mesh().V()[cells[i]]/set_.V()*massFlowRate;
        }
    }
    else if (fieldName == heName_ && fieldValues_.found(TName_))
    {
        if (fieldValues_.found(heName_))
        {
            WarningInFunction
                << "Source " << name() << " defined for both field " << heName_
                << " and " << TName_ << ". Only one of these should be present."
                << endl;
        }

        const scalar t = mesh().time().value();
        const scalar massFlowRate = massFlowRate_->value(t);
        const scalar T = fieldValues_[TName_]->value<scalar>(t);
        const basicThermo& thermo =
            mesh().lookupObject<basicThermo>
            (
                IOobject::groupName(basicThermo::dictName, phaseName_)
            );
        const scalarField hs
        (
            thermo.hs(scalarField(cells.size(), T), cells)
        );

        forAll(cells, i)
        {
            eqn.source()[cells[i]] -=
                mesh().V()[cells[i]]/set_.V()*massFlowRate*hs[i];
        }
    }
    else
    {
        addGeneralSupType(eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::massSource::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addSupType(eqn, fieldName);
}


template<class Type>
void Foam::fv::massSource::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addSupType(eqn, fieldName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::massSource::massSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    phaseName_(),
    rhoName_(),
    heName_(),
    TName_(),
    massFlowRate_(),
    fieldValues_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::massSource::addsSupToField(const word& fieldName) const
{
    const bool isThisPhase = IOobject::group(fieldName) == phaseName_;

    if
    (
        isThisPhase
     && !(fieldName == rhoName_)
     && !(fieldName == heName_ && fieldValues_.found(TName_))
     && !fieldValues_.found(fieldName)
    )
    {
        WarningInFunction
            << "No value supplied for field " << fieldName << " in "
            << type() << " fvModel " << name() << endl;

        return false;
    }

    return isThisPhase;
}


Foam::wordList Foam::fv::massSource::addSupFields() const
{
    wordList fieldNames = fieldValues_.toc();

    if (fieldValues_.found(TName_))
    {
        fieldNames[findIndex(fieldNames, TName_)] = heName_;
    }

    return fieldNames;
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_SUP, fv::massSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_SUP, fv::massSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_SUP, fv::massSource);


void Foam::fv::massSource::updateMesh(const mapPolyMesh& mpm)
{
    set_.updateMesh(mpm);
}


bool Foam::fv::massSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
