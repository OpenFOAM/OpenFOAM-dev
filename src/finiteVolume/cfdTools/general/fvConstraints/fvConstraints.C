/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvConstraints.H"
#include "fvModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvConstraints, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fvConstraints::createIOobject
(
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        typeName,
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        Info<< "Creating fvConstraints from "
            << io.instance()/io.name() << nl
            << endl;

        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        // For backward-compatibility
        // check if the fvOptions file is in system
        io.rename("fvOptions");

        if (io.headerOk())
        {
            Warning
                << "Creating fvConstraints from "
                << io.instance()/io.name() << nl
                << endl;

            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            return io;
        }
        else
        {
            // For backward-compatibility
            // check if the fvOptions file is in constant
            io.instance() = mesh.time().constant();

            if (io.headerOk())
            {
                Warning
                    << "Creating fvConstraints from "
                    << io.instance()/io.name()
                    << " rather than system/fvConstraints"
                    << endl;

                io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
                return io;
            }
            else
            {
                io.rename(typeName);
                io.readOpt() = IOobject::NO_READ;
                return io;
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvConstraints::checkApplied() const
{
    if (mesh().time().timeIndex() > checkTimeIndex_)
    {
        const PtrListDictionary<fvConstraint>& constraintList(*this);

        forAll(constraintList, i)
        {
            const fvConstraint& constraint = constraintList[i];

            wordHashSet notConstrainedFields(constraint.constrainedFields());
            notConstrainedFields -= constrainedFields_[i];

            forAllConstIter(wordHashSet, notConstrainedFields, iter)
            {
                WarningInFunction
                    << "Constraint " << constraint.name()
                    << " defined for field " << iter.key()
                    << " but never used" << endl;
            }
        }

        checkTimeIndex_ = mesh().time().timeIndex();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvConstraints::fvConstraints
(
    const fvMesh& mesh
)
:
    DemandDrivenMeshObject<fvMesh, Foam::UpdateableMeshObject, fvConstraints>
    (
        mesh,
        createIOobject(mesh)
    ),
    PtrListDictionary<fvConstraint>(0),
    checkTimeIndex_(mesh.time().timeIndex() + 1),
    constrainedFields_()
{
    readHeaderOk(IOstream::ASCII, typeName);

    const bool readFromFvConstraints(IOobject::name() == typeName);

    const dictionary& dict(*this);

    // Count number of active fvConstraints
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    PtrListDictionary<fvConstraint>::setSize(count);

    constrainedFields_.setSize(count);

    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& constraintDict = iter().dict();

            const word constraintType(constraintDict.lookup("type"));

            if
            (
                readFromFvConstraints
             || !fvModel::dictionaryConstructorTablePtr_->found
                (
                    constraintType
                )
            )
            {
                PtrListDictionary<fvConstraint>::set
                (
                    i,
                    name,
                    fvConstraint::New(name, constraintDict, mesh).ptr()
                );

                constrainedFields_.set(i, new wordHashSet());

                i++;
            }
        }
    }

    PtrListDictionary<fvConstraint>::setSize(i);
    constrainedFields_.setSize(i);

    if (readFromFvConstraints)
    {
        // Add file watch on the fvConstraints dictionary for
        // MUST_READ_IF_MODIFIED
        addWatch();
    }
    else
    {
        // If the fvOptions file was read re-name to fvConstraints
        rename(typeName);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvConstraints::constrainsField(const word& fieldName) const
{
    const PtrListDictionary<fvConstraint>& constraintList(*this);

    forAll(constraintList, i)
    {
        if (constraintList[i].constrainsField(fieldName))
        {
            return true;
        }
    }

    return false;
}


bool Foam::fvConstraints::movePoints()
{
    bool allOk = true;

    PtrListDictionary<fvConstraint>& constraintList(*this);

    forAll(constraintList, i)
    {
        allOk = allOk && constraintList[i].movePoints();
    }

    return allOk;
}


void Foam::fvConstraints::topoChange(const polyTopoChangeMap& map)
{
    PtrListDictionary<fvConstraint>& constraintList(*this);

    forAll(constraintList, i)
    {
        constraintList[i].topoChange(map);
    }
}


void Foam::fvConstraints::mapMesh(const polyMeshMap& map)
{
    PtrListDictionary<fvConstraint>& constraintList(*this);

    forAll(constraintList, i)
    {
        constraintList[i].mapMesh(map);
    }
}


void Foam::fvConstraints::distribute(const polyDistributionMap& map)
{
    PtrListDictionary<fvConstraint>& constraintList(*this);

    forAll(constraintList, i)
    {
        constraintList[i].distribute(map);
    }
}


bool Foam::fvConstraints::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::fvConstraints::writeData(Ostream& os) const
{
    dictionary::write(os, false);
    return os.good();
}


bool Foam::fvConstraints::read()
{
    if (IOdictionary::regIOobject::read())
    {
        checkTimeIndex_ = mesh().time().timeIndex() + 1;

        bool allOk = true;

        PtrListDictionary<fvConstraint>& constraintList(*this);

        forAll(constraintList, i)
        {
            fvConstraint& constraint = constraintList[i];
            bool ok = constraint.read(subDict(constraint.name()));
            allOk = (allOk && ok);
        }
        return allOk;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const fvConstraints& constraints
)
{
    constraints.writeData(os);
    return os;
}


// ************************************************************************* //
