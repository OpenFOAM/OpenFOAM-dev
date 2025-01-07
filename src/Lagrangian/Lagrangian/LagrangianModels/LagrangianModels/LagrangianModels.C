/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianModels.H"
#include "LagrangianMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LagrangianModels, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::LagrangianModels::io(const LagrangianMesh& mesh) const
{
    typeIOobject<IOdictionary> result
    (
        typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!result.headerOk())
    {
        result.readOpt() = IOobject::NO_READ;
    }

    return result;
}


void Foam::LagrangianModels::checkApplied() const
{
    const label timeIndex =
        mesh().time().subCycling()
      ? mesh().time().prevTimeState().timeIndex()
      : mesh().time().timeIndex();

    if (timeIndex > checkTimeIndex_)
    {
        const PtrListDictionary<LagrangianModel>& modelList(*this);

        forAll(modelList, i)
        {
            const LagrangianModel& model = modelList[i];

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

        checkTimeIndex_ = timeIndex;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianModels::LagrangianModels(const LagrangianMesh& mesh)
:
    DemandDrivenMeshObject
    <
        LagrangianMesh,
        TopoChangeableMeshObject,
        LagrangianModels
    >
    (
        io(mesh),
        mesh
    ),
    dictionary(),
    PtrListDictionary<LagrangianModel>(0),
    checkTimeIndex_(mesh.time().timeIndex() + 1),
    addSupFields_()
{
    readHeaderOk(IOstream::ASCII, typeName);

    const dictionary& dict(*this);

    // Iterate through the dictionary to determine the number of models
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        i += iter().isDict();
    }

    // Size the storage
    PtrListDictionary<LagrangianModel>::setSize(i);
    addSupFields_.setSize(i);

    // Iterate through the dictionary to construct the models
    i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict()) continue;

        const word& name = iter().keyword();
        const dictionary& modelDict = iter().dict();

        PtrListDictionary<LagrangianModel>::set
        (
            i,
            name,
            LagrangianModel::New(name, mesh, modelDict).ptr()
        );

        addSupFields_.set(i, new wordHashSet());

        i ++;
    }

    // Do post-construction
    {
        PtrListDictionary<LagrangianModel>& modelList(*this);

        forAll(modelList, i)
        {
            modelList[i].postConstruct();
        }
    }

    // Enable re-reading
    addWatch();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianModels::~LagrangianModels()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::LagrangianModels::addsSupToField(const word& fieldName) const
{
    const PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        const LagrangianModel& model = modelList[i];

        if (model.addsSupToField(fieldName))
        {
            return true;
        }
    }

    return false;
}


void Foam::LagrangianModels::correct()
{
    PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].correct();
    }
}


Foam::LagrangianSubMesh Foam::LagrangianModels::preModify
(
    LagrangianMesh& mesh
) const
{
    const PtrListDictionary<LagrangianModel>& modelList(*this);

    // Get the models to produce a list of which elements are to be modified
    // and which are to be removed
    DynamicList<LagrangianModel::elementModification> elementModifications;
    forAll(modelList, i)
    {
        modelList[i].preModify(mesh, elementModifications);
    }

    // Partition the mesh into modified and removed blocks of elements (and
    // unmodified elements, at the start)
    const labelList offsets = mesh.partition(2, elementModifications);

    // Return the sub-mesh associated with the modified block
    return
        LagrangianSubMesh
        (
            mesh,
            LagrangianGroup::none,
            offsets[1] - offsets[0],
            offsets[0]
        );

    // Note: We don't remove the removed elements here. The outer code might
    // need to do something with them (e.g., remove them from averages, or
    // transfer their mass to the carrier) before they get deleted. They will
    // actually get removed in LagrangianModels::modify.
}


Foam::LagrangianSubMesh Foam::LagrangianModels::modify
(
    LagrangianMesh& mesh,
    const LagrangianSubMesh& modifiedMesh
) const
{
    const PtrListDictionary<LagrangianModel>& modelList(*this);

    // Remove the removed block from the mesh
    mesh.remove(mesh.size() - modifiedMesh.end());

    // Get the models to modify and create elements
    LagrangianSubMesh modifiedAndCreatedMesh(mesh.subNone());
    modifiedAndCreatedMesh += modifiedMesh;
    forAll(modelList, i)
    {
        modifiedAndCreatedMesh += modelList[i].modify(mesh, modifiedMesh);
    }

    // Return a sub-mesh identifying all modified and created elements
    return modifiedAndCreatedMesh;
}


void Foam::LagrangianModels::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].calculate(deltaT, final);
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::LagrangianModels::source(const LagrangianSubScalarField& deltaT) const
{
    checkApplied();

    tmp<LagrangianSubScalarField> tS =
        LagrangianSubScalarField::New
        (
            word::null,
            deltaT.mesh(),
            dimensionedScalar(dimRate, 0)
        );
    LagrangianSubScalarField& S = tS.ref();

    const PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        const LagrangianModel& model = modelList[i];

        if (model.addsSupToField(word::null))
        {
            addSupFields_[i].insert(word::null);

            model.addSup(deltaT, S);
        }
    }

    return tS;
}


bool Foam::LagrangianModels::movePoints()
{
    return true;
}


void Foam::LagrangianModels::topoChange(const polyTopoChangeMap& map)
{
    PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].topoChange(map);
    }
}


void Foam::LagrangianModels::mapMesh(const polyMeshMap& map)
{
    PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].mapMesh(map);
    }
}


void Foam::LagrangianModels::distribute(const polyDistributionMap& map)
{
    PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        modelList[i].distribute(map);
    }
}


bool Foam::LagrangianModels::read()
{
    if (regIOobject::read())
    {
        checkTimeIndex_ = mesh().time().timeIndex() + 1;

        bool allOk = true;

        PtrListDictionary<LagrangianModel>& modelList(*this);

        forAll(modelList, i)
        {
            const dictionary& modelDict = subDict(modelList[i].name());

            const bool ok =
                modelList[i].read
                (
                    modelDict.optionalSubDict(modelList[i].type() + "Coeffs")
                );
            allOk = allOk && ok;
        }

        return allOk;
    }
    else
    {
        return false;
    }
}


bool Foam::LagrangianModels::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::LagrangianModels::writeData(Ostream& os) const
{
    dictionary::write(os, false);
    return os.good();
}


bool Foam::LagrangianModels::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    bool allOk = true;

    const PtrListDictionary<LagrangianModel>& modelList(*this);

    forAll(modelList, i)
    {
        const bool ok = modelList[i].write(write);
        allOk = allOk && ok;

        const bool okState = modelList[i].writeState(write);
        allOk = allOk && okState;
    }

    return allOk;
}


// ************************************************************************* //
