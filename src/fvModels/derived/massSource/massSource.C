/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    defineTypeNameAndDebug(massSourceBase, 0);
    defineTypeNameAndDebug(massSource, 0);
    addToRunTimeSelectionTable(fvModel, massSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::massSourceBase::readCoeffs()
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
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        )
    )
    {
        const basicThermo& thermo =
            mesh().lookupObject<basicThermo>
            (
                IOobject::groupName(physicalProperties::typeName, phaseName_)
            );
        heName_ = thermo.he().name();
        TName_ = thermo.T().name();
    }
}


template<class Type>
void Foam::fv::massSourceBase::addGeneralSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const labelUList cells = set_.cells();

    const scalar massFlowRate = this->massFlowRate();

    if (massFlowRate > 0)
    {
        const Type value =
            fieldValues_[fieldName]->value<Type>(mesh().time().userTimeValue());

        forAll(cells, i)
        {
            eqn.source()[cells[i]] -=
                mesh().V()[cells[i]]/set_.V()*massFlowRate*value;
        }
    }
    else
    {
        forAll(cells, i)
        {
            eqn.diag()[cells[i]] +=
                mesh().V()[cells[i]]/set_.V()*massFlowRate;
        }
    }
}


template<class Type>
void Foam::fv::massSourceBase::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addGeneralSupType(eqn, fieldName);
}


void Foam::fv::massSourceBase::addSupType
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const labelUList cells = set_.cells();

    if (fieldName == rhoName_)
    {
        const scalar massFlowRate = this->massFlowRate();

        forAll(cells, i)
        {
            eqn.source()[cells[i]] -=
                mesh().V()[cells[i]]/set_.V()*massFlowRate;
        }
    }
    else if (fieldName == heName_ && fieldValues_.found(TName_))
    {
        const scalar massFlowRate = this->massFlowRate();

        if (massFlowRate > 0)
        {
            if (fieldValues_.found(heName_))
            {
                WarningInFunction
                    << "Source " << name() << " defined for both field "
                    << heName_ << " and " << TName_
                    << ". Only one of these should be present." << endl;
            }

            const basicThermo& thermo =
                mesh().lookupObject<basicThermo>
                (
                    IOobject::groupName
                    (
                        physicalProperties::typeName,
                        phaseName_
                    )
                );

            const scalar T =
                fieldValues_[TName_]->value<scalar>
                (
                    mesh().time().userTimeValue()
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
            forAll(cells, i)
            {
                eqn.diag()[cells[i]] +=
                    mesh().V()[cells[i]]/set_.V()*massFlowRate;
            }
        }
    }
    else
    {
        addGeneralSupType(eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::massSourceBase::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addSupType(eqn, fieldName);
}


template<class Type>
void Foam::fv::massSourceBase::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    addSupType(eqn, fieldName);
}


void Foam::fv::massSource::readCoeffs()
{
    readSet();

    readFieldValues();

    massFlowRate_.reset
    (
        Function1<scalar>::New("massFlowRate", coeffs()).ptr()
    );
}


Foam::scalar Foam::fv::massSource::massFlowRate() const
{
    return massFlowRate_->value(mesh().time().userTimeValue());
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::massSourceBase::readSet()
{
    set_.read(coeffs());
}


void Foam::fv::massSourceBase::readFieldValues()
{
    fieldValues_.clear();
    const dictionary& fieldCoeffs = coeffs().subDict("fieldValues");
    forAllConstIter(dictionary, fieldCoeffs, iter)
    {
        fieldValues_.set
        (
            iter().keyword(),
            new unknownTypeFunction1(iter().keyword(), fieldCoeffs)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::massSourceBase::massSourceBase
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    phaseName_(),
    rhoName_(),
    heName_(),
    TName_(),
    set_(fvCellSet(mesh)),
    fieldValues_()
{
    readCoeffs();
}


Foam::fv::massSource::massSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    massSourceBase(name, modelType, mesh, dict),
    massFlowRate_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::massSourceBase::addsSupToField(const word& fieldName) const
{
    const bool isThisPhase = IOobject::group(fieldName) == phaseName_;

    if
    (
        isThisPhase
     && massFlowRate() > 0
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


Foam::wordList Foam::fv::massSourceBase::addSupFields() const
{
    wordList fieldNames = fieldValues_.toc();

    if (fieldValues_.found(TName_))
    {
        fieldNames[findIndex(fieldNames, TName_)] = heName_;
    }

    return fieldNames;
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_SUP, fv::massSourceBase);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_SUP, fv::massSourceBase);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_SUP, fv::massSourceBase);


bool Foam::fv::massSourceBase::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::massSourceBase::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::massSourceBase::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::massSourceBase::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::massSourceBase::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::fv::massSource::read(const dictionary& dict)
{
    if (massSourceBase::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
