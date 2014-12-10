/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "DsmcCloud.H"
#include "BinaryCollisionModel.H"
#include "WallInteractionModel.H"
#include "InflowBoundaryModel.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"
#include "polyMeshTetDecomposition.H"

using namespace Foam::constant;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] =
        typename ParcelType::constantProperties::constantProperties(molDict);
    }
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::initialise
(
    const IOdictionary& dsmcInitialiseDict
)
{
    Info<< nl << "Initialising particles" << endl;

    const scalar temperature
    (
        readScalar(dsmcInitialiseDict.lookup("temperature"))
    );

    const vector velocity(dsmcInitialiseDict.lookup("velocity"));

    const dictionary& numberDensitiesDict
    (
        dsmcInitialiseDict.subDict("numberDensities")
    );

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar> numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );
    }

    numberDensities /= nParticle_;

    forAll(mesh_.cells(), cellI)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh_,
            cellI
        );

        forAll(cellTets, tetI)
        {
            const tetIndices& cellTetIs = cellTets[tetI];

            tetPointRef tet = cellTetIs.tet(mesh_);

            scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName(molecules[i]);

                label typeId(findIndex(typeIdList_, moleculeName));

                if (typeId == -1)
                {
                    FatalErrorIn("Foam::DsmcCloud<ParcelType>::initialise")
                        << "typeId " << moleculeName << "not defined." << nl
                        << abort(FatalError);
                }

                const typename ParcelType::constantProperties& cP =
                constProps(typeId);

                scalar numberDensity = numberDensities[i];

                // Calculate the number of particles required
                scalar particlesRequired = numberDensity*tetVolume;

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.scalar01()
                )
                {
                    nParticlesToInsert++;
                }

                for (label pI = 0; pI < nParticlesToInsert; pI++)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U = equipartitionLinearVelocity
                    (
                        temperature,
                        cP.mass()
                    );

                    scalar Ei = equipartitionInternalEnergy
                    (
                        temperature,
                        cP.internalDegreesOfFreedom()
                    );

                    U += velocity;

                    addNewParcel
                    (
                        p,
                        U,
                        Ei,
                        cellI,
                        cellTetIs.face(),
                        cellTetIs.tetPt(),
                        typeId
                    );
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const typename ParcelType::constantProperties& cP = constProps
    (
        mostAbundantType
    );

    sigmaTcRMax_.internalField() = cP.sigmaT()*maxwellianMostProbableSpeed
    (
        temperature,
        cP.mass()
    );

    sigmaTcRMax_.correctBoundaryConditions();
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::collisions()
{
    if (!binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    scalar deltaT = mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    forAll(cellOccupancy_, cellI)
    {
        const DynamicList<ParcelType*>& cellParcels(cellOccupancy_[cellI]);

        label nC(cellParcels.size());

        if (nC > 1)
        {

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            const point& cC = mesh_.cellCentres()[cellI];

            forAll(cellParcels, i)
            {
                const ParcelType& p = *cellParcels[i];

                vector relPos = p.position() - cC;

                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                subCells[subCell].append(i);

                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = sigmaTcRMax_[cellI];

            scalar selectedPairs =
                collisionSelectionRemainder_[cellI]
              + 0.5*nC*(nC - 1)*nParticle_*sigmaTcRMax*deltaT
               /mesh_.cellVolumes()[cellI];

            label nCandidates(selectedPairs);

            collisionSelectionRemainder_[cellI] = selectedPairs - nCandidates;

            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen_.integer(0, nC - 1);

                // Declare the second collision candidate
                label candidateQ = -1;

                List<label> subCellPs = subCells[whichSubCell[candidateP]];

                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        candidateQ = subCellPs[rndGen_.integer(0, nSC - 1)];

                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        candidateQ = rndGen_.integer(0, nC - 1);

                    } while (candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = rndGen_.integer(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = rndGen_.integer(0, nC-1);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = rndGen_.integer(0, nC-1);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                ParcelType& parcelP = *cellParcels[candidateP];
                ParcelType& parcelQ = *cellParcels[candidateQ];

                scalar sigmaTcR = binaryCollision().sigmaTcR
                (
                    parcelP,
                    parcelQ
                );

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this

                if (sigmaTcR > sigmaTcRMax_[cellI])
                {
                    sigmaTcRMax_[cellI] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > rndGen_.scalar01())
                {
                    binaryCollision().collide
                    (
                        parcelP,
                        parcelQ
                    );

                    collisions++;
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    sigmaTcRMax_.correctBoundaryConditions();

    if (collisionCandidates)
    {
        Info<< "    Collisions                      = "
            << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl
            << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::resetFields()
{
    q_ = dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0);

    fD_ = dimensionedVector
    (
        "zero",
        dimensionSet(1, -1, -2, 0, 0),
        vector::zero
    );

    rhoN_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL);

    rhoM_ =  dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), VSMALL);

    dsmcRhoN_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0);

    linearKE_ = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    internalE_ = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    iDof_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL);

    momentum_ = dimensionedVector
    (
        "zero",
        dimensionSet(1, -2, -1, 0, 0),
        vector::zero
    );
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::calculateFields()
{
    scalarField& rhoN = rhoN_.internalField();

    scalarField& rhoM = rhoM_.internalField();

    scalarField& dsmcRhoN = dsmcRhoN_.internalField();

    scalarField& linearKE = linearKE_.internalField();

    scalarField& internalE = internalE_.internalField();

    scalarField& iDof = iDof_.internalField();

    vectorField& momentum = momentum_.internalField();

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        rhoN[cellI]++;

        rhoM[cellI] += constProps(p.typeId()).mass();

        dsmcRhoN[cellI]++;

        linearKE[cellI] += 0.5*constProps(p.typeId()).mass()*(p.U() & p.U());

        internalE[cellI] += p.Ei();

        iDof[cellI] += constProps(p.typeId()).internalDegreesOfFreedom();

        momentum[cellI] += constProps(p.typeId()).mass()*p.U();
    }

    rhoN *= nParticle_/mesh().cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoM *= nParticle_/mesh().cellVolumes();
    rhoM_.correctBoundaryConditions();

    dsmcRhoN_.correctBoundaryConditions();

    linearKE *= nParticle_/mesh().cellVolumes();
    linearKE_.correctBoundaryConditions();

    internalE *= nParticle_/mesh().cellVolumes();
    internalE_.correctBoundaryConditions();

    iDof *= nParticle_/mesh().cellVolumes();
    iDof_.correctBoundaryConditions();

    momentum *= nParticle_/mesh().cellVolumes();
    momentum_.correctBoundaryConditions();
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const vector& U,
    const scalar Ei,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId
)
{
    ParcelType* pPtr = new ParcelType
    (
        mesh_,
        position,
        U,
        Ei,
        cellI,
        tetFaceI,
        tetPtI,
        typeId
    );

    this->addParticle(pPtr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DsmcCloud<ParcelType>::DsmcCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    Cloud<ParcelType>(mesh, cloudName, false),
    DsmcBaseCloud(),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(mesh_.nCells()),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    collisionSelectionRemainder_(mesh_.nCells(), 0),
    q_
    (
        IOobject
        (
            "q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    fD_
    (
        IOobject
        (
            "fD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    dsmcRhoN_
    (
        IOobject
        (
            "dsmcRhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    linearKE_
    (
        IOobject
        (
            "linearKE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    internalE_
    (
        IOobject
        (
            "internalE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    iDof_
    (
        IOobject
        (
            "iDof",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    momentum_
    (
        IOobject
        (
            "momentum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    constProps_(),
    rndGen_(label(149382906) + 7183*Pstream::myProcNo()),
    boundaryT_
    (
        volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    boundaryU_
    (
        volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    binaryCollisionModel_
    (
        BinaryCollisionModel<DsmcCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    wallInteractionModel_
    (
        WallInteractionModel<DsmcCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    inflowBoundaryModel_
    (
        InflowBoundaryModel<DsmcCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    )
{
    buildConstProps();

    buildCellOccupancy();

    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, i)
    {
        collisionSelectionRemainder_[i] = rndGen_.scalar01();
    }

    if (readFields)
    {
        ParcelType::readFields(*this);
    }
}


template<class ParcelType>
Foam::DsmcCloud<ParcelType>::DsmcCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& dsmcInitialiseDict
)
    :
    Cloud<ParcelType>(mesh, cloudName, false),
    DsmcBaseCloud(),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 3, -1, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    collisionSelectionRemainder_(),
    q_
    (
        IOobject
        (
            this->name() + "q_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    fD_
    (
        IOobject
        (
            this->name() + "fD_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    rhoN_
    (
        IOobject
        (
            this->name() + "rhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL)
    ),
    rhoM_
    (
        IOobject
        (
            this->name() + "rhoM_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), VSMALL)
    ),
    dsmcRhoN_
    (
        IOobject
        (
            this->name() + "dsmcRhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0)
    ),
    linearKE_
    (
        IOobject
        (
            this->name() + "linearKE_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    internalE_
    (
        IOobject
        (
            this->name() + "internalE_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    iDof_
    (
        IOobject
        (
            this->name() + "iDof_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL)
    ),
    momentum_
    (
        IOobject
        (
            this->name() + "momentum_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -2, -1, 0, 0),
            vector::zero
        )
    ),
    constProps_(),
    rndGen_(label(971501) + 1526*Pstream::myProcNo()),
    boundaryT_
    (
        volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",  dimensionSet(0, 0, 0, 1, 0), 0.0)
        )
    ),
    boundaryU_
    (
        volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                dimensionSet(0, 1, -1, 0, 0),
                vector::zero
            )
        )
    ),
    binaryCollisionModel_(),
    wallInteractionModel_(),
    inflowBoundaryModel_()
{
    clear();

    buildConstProps();

    initialise(dsmcInitialiseDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DsmcCloud<ParcelType>::~DsmcCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::evolve()
{
    typename ParcelType::trackingData td(*this);

    // Reset the data collection fields
    resetFields();

    if (debug)
    {
        this->dumpParticlePositions();
    }

    // Insert new particles from the inflow boundary
    this->inflowBoundary().inflow();

    // Move the particles ballistically with their current velocities
    Cloud<ParcelType>::move(td, mesh_.time().deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    // Calculate new velocities via stochastic collisions
    collisions();

    // Calculate the volume field data
    calculateFields();
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::info() const
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    scalar nMol = nDsmcParticles*nParticle_;

    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar internalEnergy = internalEnergyOfSystem();
    reduce(internalEnergy, sumOp<scalar>());

    Info<< "Cloud name: " << this->name() << nl
        << "    Number of dsmc particles        = "
        << nDsmcParticles
        << endl;

    if (nDsmcParticles)
    {
        Info<< "    Number of molecules             = "
            << nMol << nl
            << "    Mass in system                  = "
            << returnReduce(massInSystem(), sumOp<scalar>()) << nl
            << "    Average linear momentum         = "
            << linearMomentum/nMol << nl
            << "   |Average linear momentum|        = "
            << mag(linearMomentum)/nMol << nl
            << "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average internal energy         = "
            << internalEnergy/nMol << nl
            << "    Average total energy            = "
            << (internalEnergy + linearKineticEnergy)/nMol
            << endl;
    }
}


template<class ParcelType>
Foam::vector Foam::DsmcCloud<ParcelType>::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *vector
        (
            rndGen_.GaussNormal(),
            rndGen_.GaussNormal(),
            rndGen_.GaussNormal()
        );
}


template<class ParcelType>
Foam::scalar Foam::DsmcCloud<ParcelType>::equipartitionInternalEnergy
(
    scalar temperature,
    scalar iDof
)
{
    scalar Ei = 0.0;

    if (iDof < SMALL)
    {
        return Ei;
    }
    else if (iDof < 2.0 + SMALL && iDof > 2.0 - SMALL)
    {
        // Special case for iDof = 2, i.e. diatomics;
        Ei = -log(rndGen_.scalar01())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*iDof - 1;

        scalar energyRatio;

        scalar P = -1;

        do
        {
            energyRatio = 10*rndGen_.scalar01();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < rndGen_.scalar01());

        Ei = energyRatio*physicoChemical::k.value()*temperature;
    }

    return Ei;
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}


// ************************************************************************* //
