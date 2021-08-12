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

#include "fvModels.H"
#include "fvConstraint.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvModels, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fvModels::createIOobject
(
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        Info<< "Creating fvModels from "
            << io.instance()/io.name() << nl
            << endl;

        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        // For backward-compatibility
        // check if the fvOptions file is in constant
        io.rename("fvOptions");

        if (io.headerOk())
        {
            Warning
                << "Creating fvModels from "
                << io.instance()/io.name() << nl
                << endl;

            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            return io;
        }
        else
        {
            // For backward-compatibility
            // check if the fvOptions file is in system
            io.instance() = mesh.time().system();

            if (io.headerOk())
            {
                Warning
                    << "Creating fvModels from "
                    << io.instance()/io.name()
                    << " rather than constant/fvModels"
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

void Foam::fvModels::checkApplied() const
{
    if (mesh().time().timeIndex() > checkTimeIndex_)
    {
        const PtrListDictionary<fvModel>& modelList(*this);

        forAll(modelList, i)
        {
            const fvModel& model = modelList[i];

            wordHashSet notAddSupFields(model.addSupFields());
            notAddSupFields -= addSupFields_[i];

            forAllConstIter(wordHashSet, notAddSupFields, iter)
            {
                WarningInFunction
                    << "Model " << model.name()
                    << " defined for field " << iter.key()
                    << " but never used" << endl;
            }
        }

        checkTimeIndex_ = mesh().time().timeIndex();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvModels::fvModels
(
    const fvMesh& mesh
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, fvModels>
    (
        mesh,
        createIOobject(mesh)
    ),
    PtrListDictionary<fvModel>(0),
    checkTimeIndex_(mesh.time().timeIndex() + 1),
    addSupFields_()
{
    readHeaderOk(IOstream::ASCII, typeName);

    const bool readFromFvModels(IOobject::name() == typeName);

    const dictionary& dict(*this);

    // Count number of active fvModels
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    PtrListDictionary<fvModel>::setSize(count);

    addSupFields_.setSize(count);

    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            const word modelType(modelDict.lookup("type"));

            if
            (
                readFromFvModels
             || !fvConstraint::dictionaryConstructorTablePtr_->found
                (
                    modelType
                )
            )
            {
                PtrListDictionary<fvModel>::set
                (
                    i,
                    name,
                    fvModel::New(name, modelDict, mesh).ptr()
                );

                addSupFields_.set(i, new wordHashSet());

                i++;
            }
        }
    }

    PtrListDictionary<fvModel>::setSize(i);
    addSupFields_.setSize(i);

    if (readFromFvModels)
    {
        // Add file watch on the fvModels dictionary for
        // MUST_READ_IF_MODIFIED
        addWatch();
    }
    else
    {
        // If the fvOptions file was read re-name to fvModels
        rename(typeName);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvModels::addsSupToField(const word& fieldName) const
{
    const PtrListDictionary<fvModel>& modelList(*this);

    forAll(modelList, i)
    {
        if (modelList[i].addsSupToField(fieldName))
        {
            return true;
        }
    }

    return false;
}


void Foam::fvModels::preUpdateMesh()
{
    PtrListDictionary<fvModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].preUpdateMesh();
    }
}


void Foam::fvModels::updateMesh(const mapPolyMesh& mpm)
{
    PtrListDictionary<fvModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].updateMesh(mpm);
    }
}


bool Foam::fvModels::movePoints()
{
    bool allOk = true;

    PtrListDictionary<fvModel>& modelList(*this);

    forAll(modelList, i)
    {
        allOk = allOk && modelList[i].movePoints();
    }

    return allOk;
}


bool Foam::fvModels::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::fvModels::writeData(Ostream& os) const
{
    dictionary::write(os, false);
    return os.good();
}


bool Foam::fvModels::read()
{
    if (IOdictionary::regIOobject::read())
    {
        checkTimeIndex_ = mesh().time().timeIndex() + 1;

        bool allOk = true;

        PtrListDictionary<fvModel>& modelList(*this);

        forAll(modelList, i)
        {
            fvModel& model = modelList[i];
            bool ok = model.read(subDict(model.name()));
            allOk = (allOk && ok);
        }
        return allOk;
    }
    else
    {
        return false;
    }
}


void Foam::fvModels::correct()
{
    PtrListDictionary<fvModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].correct();
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const fvModels& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
