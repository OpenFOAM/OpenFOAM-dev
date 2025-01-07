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

#include "cloud.H"
#include "cloudFunctionObjectUList.H"
#include "cloudVelocityLagrangianPatch.H"
#include "processorLagrangianPatch.H"
#include "LagrangianFields.H"
#include "LagrangianmDdt.H"
#include "LagrangianSubFields.H"
#include "dimensionedTypes.H"
#include "pimpleNoLoopControl.H"
#include "Time.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloud, 0);
    defineRunTimeSelectionTable(cloud, polyMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::LagrangianMesh& Foam::cloud::mesh
(
    const polyMesh& pMesh,
    const word& name,
    const IOobject::readOption readOption,
    const IOobject::writeOption writeOption
)
{
    if (!pMesh.foundObject<LagrangianMesh>(name))
    {
        wordList wantedPatchTypes(pMesh.boundaryMesh().size());

        forAll(pMesh.boundaryMesh(), patchi)
        {
            const polyPatch& patch = pMesh.boundaryMesh()[patchi];

            wantedPatchTypes[patchi] =
                polyPatch::constraintType(patch.type())
              ? patch.type()
              : cloudVelocityLagrangianPatch::typeName;
        }

        LagrangianMesh* lMesh =
            new LagrangianMesh
            (
                pMesh,
                name,
                wantedPatchTypes,
                readOption,
                writeOption
            );

        lMesh->store();
    }

    return pMesh.lookupObjectRef<LagrangianMesh>(name);
}


#define ACCESS_STATE_FIELDS(Type, nullArg)                                     \
namespace Foam                                                                 \
{                                                                              \
    template<>                                                                 \
    PtrList<Foam::CloudStateField<Type>>& cloud::stateFields() const           \
    {                                                                          \
        return CAT3(state, CAPITALIZE(Type), Fields_);                         \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_STATE_FIELDS)
#undef ACCESS_STATE_FIELDS


void Foam::cloud::clearStateFields()
{
    #define CLEAR_TYPE_STATE_FIELDS(Type, nullArg)                             \
        forAll(stateFields<Type>(), i)                                         \
        {                                                                      \
            stateFields<Type>()[i].clear();                                    \
        }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_STATE_FIELDS);
    #undef CLEAR_TYPE_STATE_FIELDS
}


#define ACCESS_DERIVED_FIELDS(Type, nullArg)                                   \
namespace Foam                                                                 \
{                                                                              \
    template<>                                                                 \
    PtrList<Foam::CloudDerivedField<Type>>& cloud::derivedFields() const       \
    {                                                                          \
        return CAT3(derived, CAPITALIZE(Type), Fields_);                       \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_DERIVED_FIELDS)
#undef ACCESS_DERIVED_FIELDS


void Foam::cloud::clearDerivedFields(const bool final)
{
    #define CLEAR_TYPE_DERIVED_FIELDS(Type, nullArg)                           \
        forAll(derivedFields<Type>(), i)                                       \
        {                                                                      \
            derivedFields<Type>()[i].clear(final);                             \
        }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_DERIVED_FIELDS);
    #undef CLEAR_TYPE_DERIVED_FIELDS
}


#define ACCESS_AVERAGE_FIELDS(Type, nullArg)                                   \
namespace Foam                                                                 \
{                                                                              \
    template<>                                                                 \
    PtrList<Foam::CloudAverageField<Type>>& cloud::averageFields() const       \
    {                                                                          \
        return CAT3(average, CAPITALIZE(Type), Fields_);                       \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_AVERAGE_FIELDS)
#undef ACCESS_AVERAGE_FIELDS


void Foam::cloud::removeFromAverageFields(const LagrangianSubMesh& subMesh)
{
    #define REMOVE_FROM_TYPE_AVERAGE_FIELDS(Type, nullArg)                     \
        forAll(averageFields<Type>(), i)                                       \
        {                                                                      \
            averageFields<Type>()[i].remove(subMesh);                          \
        }
    FOR_ALL_FIELD_TYPES(REMOVE_FROM_TYPE_AVERAGE_FIELDS);
    #undef REMOVE_FROM_TYPE_AVERAGE_FIELDS
}


void Foam::cloud::addToAverageFields
(
    const LagrangianSubMesh& subMesh,
    const bool final
)
{
    #define ADD_TO_TYPE_AVERAGE_FIELDS(Type, nullArg)                          \
        forAll(averageFields<Type>(), i)                                       \
        {                                                                      \
            averageFields<Type>()[i].add(subMesh, !final);                     \
        }
    FOR_ALL_FIELD_TYPES(ADD_TO_TYPE_AVERAGE_FIELDS);
    #undef ADD_TO_TYPE_AVERAGE_FIELDS
}


void Foam::cloud::correctAverageFields
(
    const LagrangianSubMesh& subMesh,
    const bool final
)
{
    #define CORRECT_TYPE_AVERAGE_FIELDS(Type, nullArg)                         \
        forAll(averageFields<Type>(), i)                                       \
        {                                                                      \
            averageFields<Type>()[i].correct(subMesh, !final);                 \
        }
    FOR_ALL_FIELD_TYPES(CORRECT_TYPE_AVERAGE_FIELDS);
    #undef CORRECT_TYPE_AVERAGE_FIELDS
}


void Foam::cloud::clearAverageFields()
{
    #define CLEAR_TYPE_AVERAGE_FIELDS(Type, nullArg)                           \
        forAll(averageFields<Type>(), i)                                       \
        {                                                                      \
            averageFields<Type>()[i].clear(true);                              \
        }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_AVERAGE_FIELDS);
    #undef CLEAR_TYPE_AVERAGE_FIELDS
}


void Foam::cloud::resetAverageFields()
{
    #define RESET_TYPE_AVERAGE_FIELDS(Type, nullArg)                           \
        forAll(averageFields<Type>(), i)                                       \
        {                                                                      \
            averageFields<Type>()[i].reset();                                  \
        }
    FOR_ALL_FIELD_TYPES(RESET_TYPE_AVERAGE_FIELDS);
    #undef RESET_TYPE_AVERAGE_FIELDS
}


Foam::IOobject Foam::cloud::stateIo(const IOobject::readOption r) const
{
    return
        IOobject
        (
            "state",
            time().name(),
            mesh(),
            r,
            IOobject::AUTO_WRITE
        );
}


Foam::autoPtr<Foam::LagrangianLabelInternalField>
Foam::cloud::readStates() const
{
    typeIOobject<LagrangianLabelInternalField> stateIo
    (
        this->stateIo(IOobject::MUST_READ)
    );

    return
        autoPtr<Foam::LagrangianLabelInternalField>
        (
            stateIo.headerOk()
          ? new LagrangianLabelInternalField(stateIo, mesh_)
          : nullptr
        );
}


Foam::autoPtr<Foam::List<Foam::LagrangianState>>
Foam::cloud::initialStates() const
{
    // Return an empty pointer of we have no states stored
    if (!statePtr_.valid()) return autoPtr<List<LagrangianState>>();

    // Allocate a list of states
    autoPtr<List<LagrangianState>> resultPtr
    (
        new List<LagrangianState>(mesh_.size(), LagrangianState::inCell)
    );
    List<LagrangianState>& result = resultPtr();

    // Convert from the state labels to state enumerations
    forAll(mesh_, i)
    {
        const LagrangianState s =
            static_cast<LagrangianState>(statePtr_()[i]);

        if (s != LagrangianState::none)
        {
            result[i] = s;
        }
    }

    return resultPtr;
}


bool Foam::cloud::storeStates()
{
    // Determine whether states are needed
    bool needStates = false;

    forAll(mesh_.boundary(), patchi)
    {
        const LagrangianSubMesh& patchMesh =
            mesh_.boundary()[patchi].mesh();

        if (!patchMesh.empty()) needStates = true;
    }

    reduce(needStates, orOp<bool>());

    // Create or destroy the state labels field as appropriate
    if (needStates && !statePtr_.valid())
    {
        statePtr_.set
        (
            new LagrangianLabelInternalField
            (
                stateIo(IOobject::NO_READ),
                mesh_,
                dimensioned<label>
                (
                    dimless,
                    static_cast<label>(LagrangianState::none)
                )
            )
        );
    }
    if (!needStates && statePtr_.valid())
    {
        statePtr_.clear();
    }

    return needStates;
}


void Foam::cloud::storeStates(const LagrangianSubMesh& subMesh)
{
    const LagrangianState state = groupToState(subMesh.group());

    forAll(subMesh, subi)
    {
        const label i = subMesh.start() + subi;

        // If this particle is complete, then store the sub-mesh state within
        // the retained state labels. If it is not complete, then reset its
        // state so that it continues to be associated with this sub-mesh.
        if (mesh_.states()[i] == LagrangianState::complete)
        {
            statePtr_()[i] = static_cast<label>(state);
        }
        else
        {
            mesh_.states()[i] = state;
        }
    }
}


Foam::tmp<Foam::LagrangianSubScalarField> Foam::cloud::cellLengthScale
(
    const LagrangianSubMesh& subMesh
) const
{
    return
        LagrangianSubScalarField::New
        (
            "cellLengthScale",
            subMesh,
            dimLength,
            scalarField(cellLengthScaleVf_, subMesh.sub(mesh_.celli()))
        );
}


void Foam::cloud::track
(
    LagrangianSubScalarSubField& fraction,
    const scalar maxTimeStepFraction,
    const scalar maxCellLengthScaleFraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    // Initialise tracking to completion over the remaining fraction
    List<LagrangianState> endState(subMesh.size(), LagrangianState::complete);
    LagrangianSubScalarField deltaFraction(1 - fraction);

    // Evaluate the displacement over the entire time-step
    const LagrangianSubVectorField displacement
    (
        time().deltaT()*U(subMesh)
    );

    // Limit the track to the maximum allowed fraction of the time-step
    if (maxTimeStepFraction < 1)
    {
        forAll(deltaFraction, subi)
        {
            if (deltaFraction[subi] > maxTimeStepFraction)
            {
                endState[subi] = LagrangianState::inCell;
                deltaFraction[subi] = maxTimeStepFraction;
            }
        }
    }

    // Limit the track to the maximum allowed fraction of the cell length scale
    if (maxCellLengthScaleFraction < rootGreat)
    {
        const LagrangianSubScalarField maxDFraction
        (
            maxCellLengthScaleFraction
           *cellLengthScale(subMesh)
           /max
            (
                mag(displacement),
                dimensionedScalar(dimLength, rootVSmall)
            )
        );

        forAll(deltaFraction, subi)
        {
            if (deltaFraction[subi] > maxDFraction[subi])
            {
                endState[subi] = LagrangianState::inCell;
                deltaFraction[subi] = maxDFraction[subi];
            }
        }
    }

    // Track the particles
    switch (tracking)
    {
        case trackingType::linear:
            mesh_.track
            (
                endState,
                LagrangianMesh::linearDisplacement
                (
                    deltaFraction*displacement
                ),
                deltaFraction,
                fraction
            );
            break;

        case trackingType::parabolic:
        {
            mesh_.track
            (
                endState,
                LagrangianMesh::parabolicDisplacement
                (
                    deltaFraction*displacement,
                    sqr(deltaFraction*time().deltaT())/2*dUdt(subMesh)
                ),
                deltaFraction,
                fraction
            );
            break;
        }
    }
}


bool Foam::cloud::writeData(Ostream&) const
{
    NotImplemented;
    return false;
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::cloud::initialise(const bool predict)
{
    clearStateFields();
    clearDerivedFields(true);
    clearAverageFields();
}


void Foam::cloud::partition()
{
    clearStateFields();
    clearDerivedFields(true);
    clearAverageFields();
}


void Foam::cloud::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<enum cloud::trackingType, 2>::names[] =
        {"linear", "parabolic"};

    const NamedEnum<enum cloud::trackingType, 2> cloudTrackingNames;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloud::cloud
(
    const polyMesh& pMesh,
    const word& name,
    const contextType context,
    const IOobject::readOption readOption,
    const IOobject::writeOption writeOption
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            pMesh.time().name(),
            mesh(pMesh, name, readOption, writeOption)
        )
    ),
    mesh_(mesh(pMesh, name, readOption, writeOption)),
    LagrangianModelsPtr_(nullptr),
    statePtr_(readStates()),
    cellLengthScaleVf_(mag(cbrt(mesh_.mesh().cellVolumes()))),
    context(context),
    tracking
    (
        cloudTrackingNames
        [
            mesh().schemes().schemesDict().lookup<word>("tracking")
        ]
    ),
    U
    (
        stateField<vector>
        (
            IOobject
            (
                "U",
                time().name(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cloud> Foam::cloud::New
(
    const polyMesh& pMesh,
    const word& name,
    const word& type
)
{
    Info<< "Selecting " << typeName
        << " with name " << name
        << " of type " << type << endl;

    if (!polyMeshConstructorTablePtr_)
    {
        FatalErrorInFunction
            << typeName << "s table is empty"
            << exit(FatalError);
    }

    polyMeshConstructorTable::iterator cstrIter;

    cstrIter = polyMeshConstructorTablePtr_->find(type);

    if (cstrIter == polyMeshConstructorTablePtr_->end())
    {
        libs.open("lib" + type + typeName.capitalise() + ".so");
    }

    cstrIter = polyMeshConstructorTablePtr_->find(type);

    if (cstrIter == polyMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << type << nl << nl
            << "Valid " << typeName << "s are :" << endl
            << polyMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<cloud> cloudPtr(cstrIter()(pMesh, name, contextType::unknown));

    // Ensure LagrangianModels are constructed before time is incremented
    cloudPtr->LagrangianModels();

    return cloudPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cloud::~cloud()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::cloud& Foam::cloud::lookup(const LagrangianMesh& mesh)
{
    return mesh.lookupObject<cloud>(cloud::typeName);
}


Foam::LagrangianModels& Foam::cloud::LagrangianModels() const
{
    if (!LagrangianModelsPtr_)
    {
        LagrangianModelsPtr_ = &Foam::LagrangianModels::New(mesh_);
    }

    return *LagrangianModelsPtr_;
}


void Foam::cloud::solve()
{
    // Create the functions list
    cloudFunctionObjectUList functions(*this);

    // Handle outer correctors
    bool predict = false;
    if (context == contextType::fvModel)
    {
        if
        (
            mesh_.mesh().foundObject<pimpleNoLoopControl>
            (
                solutionControl::typeName
            )
        )
        {
            const pimpleNoLoopControl& pimple =
                mesh_.mesh().lookupObject<pimpleNoLoopControl>
                (
                    solutionControl::typeName
                );

            if (pimple.nCorr() > 1)
            {
                if (mesh_.solution().lookup<bool>("outerCorrectors"))
                {
                    mesh_.reset(pimple.firstIter(), pimple.finalIter());
                    resetAverageFields();
                }
                else if (!pimple.firstIter())
                {
                    return;
                }
            }

            predict = pimple.firstIter();
        }
        else
        {
            predict = true;
        }
    }

    // Initial reset of cached objects
    initialise(predict);

    Info<< "Solving cloud " << mesh_.name() << ':' << endl << incrIndent;

    // Get solution controls
    const scalar maxTimeStepFraction =
        mesh_.solution().lookup<scalar>("maxTimeStepFraction");
    const scalar maxCellLengthScaleFraction =
        mesh_.solution().lookup<scalar>("maxCellLengthScaleFraction");
    const label nCorrectors = mesh_.solution().lookup<label>("nCorrectors");

    // Correct the models
    LagrangianModels().correct();

    // Initialise the tracked fraction to zero, representing all particles
    // being at the start of the time-step
    LagrangianScalarInternalDynamicField fraction
    (
        IOobject
        (
            LagrangianMesh::fractionName,
            time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0)
    );

    // Let the models do any instantaneous modifications, removals and
    // injections/creations of existing particles
    LagrangianSubMesh preModifiedMesh =
        LagrangianModels().preModify(mesh_);
    removeFromAverageFields(preModifiedMesh);
    LagrangianSubMesh modifiedMesh =
        LagrangianModels().modify(mesh_, preModifiedMesh);
    addToAverageFields(modifiedMesh, true);

    // Do a calculation on modified particles (if necessary)
    if (reCalculateModified())
    {
        removeFromAverageFields(modifiedMesh);

        const LagrangianSubScalarField zeroDeltaT
        (
            IOobject("zeroDeltaT", mesh().time().name(), mesh()),
            modifiedMesh,
            dimensionedScalar(dimTime, scalar(0))
        );

        addToAverageFields(modifiedMesh, false);

        LagrangianModels().calculate(zeroDeltaT, true);
        calculate(zeroDeltaT, true);
        functions.calculate(zeroDeltaT, true);

        correctAverageFields(modifiedMesh, true);
    }

    // Construct the object defining the scope of the Lagrangian mesh changes
    autoPtr<List<LagrangianState>> statesPtr = initialStates();
    LagrangianMesh::changer changer =
        statesPtr.valid()
      ? LagrangianMesh::changer(mesh_, statesPtr())
      : LagrangianMesh::changer(mesh_, LagrangianState::inCell);
    statesPtr.clear();

    // If we have initial states then we need to evaluate the
    // derived/non-constraint boundary conditions so that any affected
    // properties are up date. This evaluation of the velocity boundary fields
    // is what ensures that "stuck" particles move consistently with their
    // respective patches. In theory, though, it could apply to any field.
    #define EVAL_TYPE_DERIVED_PATCH_FIELDS(Type, GeoField)                     \
    {                                                                          \
        HashTable<GeoField<Type>*> fields                                      \
        (                                                                      \
            mesh_.lookupClass<GeoField<Type>>()                                \
        );                                                                     \
                                                                               \
        forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)          \
        {                                                                      \
            forAll(mesh_.boundary(), patchi)                                   \
            {                                                                  \
                const LagrangianPatch& patch = mesh_.boundary()[patchi];       \
                                                                               \
                if                                                             \
                (                                                              \
                    patch.mesh().size()                                        \
                 && !polyPatch::constraintType(patch.type())                   \
                )                                                              \
                {                                                              \
                    iter()->boundaryFieldRef()[patchi].evaluate                \
                    (                                                          \
                        NullObjectNonConstRef<PstreamBuffers>(),               \
                        fraction                                               \
                    );                                                         \
                }                                                              \
            }                                                                  \
        }                                                                      \
    }
    FOR_ALL_FIELD_TYPES(EVAL_TYPE_DERIVED_PATCH_FIELDS, LagrangianField);
    FOR_ALL_FIELD_TYPES(EVAL_TYPE_DERIVED_PATCH_FIELDS, LagrangianDynamicField);
    #undef EVAL_TYPE_DERIVED_PATCH_FIELDS

    // Iterate whilst there are incomplete particles
    while
    (
        returnReduce
        (
            mesh_.sub(LagrangianGroup::complete).size() != mesh_.size(),
            orOp<bool>()
        )
    )
    {
        // Internal tracking and calculation
        {
            LagrangianSubMesh internalMesh
            (
                mesh_.sub(LagrangianGroup::inInternalMesh)
            );

            LagrangianSubScalarSubField internalFraction
            (
                internalMesh.sub(fraction)
            );

            removeFromAverageFields(internalMesh);

            track
            (
                internalFraction,
                maxTimeStepFraction,
                maxCellLengthScaleFraction
            );

            const LagrangianSubScalarField deltaT
            (
                (internalFraction - internalFraction.oldTime())
               *mesh_.time().deltaT()
            );

            addToAverageFields(internalMesh, false);

            for (label i = 0; i <= nCorrectors; ++ i)
            {
                const bool final = i == nCorrectors;

                LagrangianModels().calculate(deltaT, final);
                calculate(deltaT, final);
                functions.calculate(deltaT, final);

                clearDerivedFields(final);
                correctAverageFields(internalMesh, final);
            }

            clearStateFields();
        }

        // Boundary tracking and calculation (if necessary)
        if (storeStates())
        {
            const labelList subMeshGlobalSizes = mesh_.subMeshGlobalSizes();

            forAll(mesh_.boundary(), patchi)
            {
                static const label onPatchZeroi =
                    static_cast<label>(LagrangianGroup::onPatchZero);

                if (subMeshGlobalSizes[onPatchZeroi + patchi] <= 0) continue;

                LagrangianSubMesh patchMesh
                (
                    mesh_.boundary()[patchi].mesh()
                );

                LagrangianSubScalarSubField patchFraction
                (
                    patchMesh.sub(fraction)
                );

                removeFromAverageFields(patchMesh);

                track
                (
                    patchFraction,
                    maxTimeStepFraction,
                    maxCellLengthScaleFraction
                );

                storeStates(patchMesh);

                const LagrangianSubScalarField deltaT
                (
                    (patchFraction - patchFraction.oldTime())
                   *mesh_.time().deltaT()
                );

                addToAverageFields(patchMesh, false);

                for (label i = 0; i <= nCorrectors; ++ i)
                {
                    const bool final = i == nCorrectors;

                    LagrangianModels().calculate(deltaT, final);
                    calculate(deltaT, final);
                    functions.calculate(deltaT, final);

                    clearDerivedFields(final);
                    correctAverageFields(patchMesh, final);
                }

                clearStateFields();
            }
        }

        // Intermediate partitioning
        mesh_.partition();
        partition();

        removeFromAverageFields(mesh_.subIncomplete());

        // Cross the faces
        functions.preCrossFaces(fraction);
        mesh_.crossFaces(fraction);
        functions.postCrossFaces(fraction);

        addToAverageFields(mesh_.subIncomplete(), true);

        // Final partitioning
        mesh_.partition();
        partition();
    };

    Info<< decrIndent;
}


void Foam::cloud::storePosition()
{
    mesh_.storePosition();
}


void Foam::cloud::movePoints(const polyMesh&)
{
    cellLengthScaleVf_ = mag(cbrt(mesh_.mesh().cellVolumes()));
}


void Foam::cloud::topoChange(const polyTopoChangeMap& map)
{
    mesh_.topoChange(map);

    cellLengthScaleVf_ = mag(cbrt(mesh_.mesh().cellVolumes()));
}


void Foam::cloud::mapMesh(const polyMeshMap& map)
{
    mesh_.mapMesh(map);

    cellLengthScaleVf_ = mag(cbrt(mesh_.mesh().cellVolumes()));
}


void Foam::cloud::distribute(const polyDistributionMap& map)
{
    mesh_.distribute(map);

    cellLengthScaleVf_ = mag(cbrt(mesh_.mesh().cellVolumes()));
}


// ************************************************************************* //
