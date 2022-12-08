/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "viewFactor.H"
#include "surfaceFields.H"
#include "constants.H"
#include "greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(viewFactor, 0);
    addToRadiationRunTimeSelectionTables(viewFactor);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiationModels::viewFactor::initialise()
{
    const polyBoundaryMesh& coarsePatches = coarseMesh_.boundaryMesh();
    const volScalarField::Boundary& qrp = qr_.boundaryField();

    label count = 0;
    forAll(qrp, patchi)
    {
        // const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const fvPatchScalarField& qrPatchi = qrp[patchi];

        if ((isA<fixedValueFvPatchScalarField>(qrPatchi)))
        {
            selectedPatches_[count] = qrPatchi.patch().index();
            nLocalCoarseFaces_ += coarsePatches[patchi].size();
            count++;
        }
    }

    selectedPatches_.resize(count--);

    if (debug)
    {
        Pout<< "radiationModels::viewFactor::initialise() "
            << "Selected patches:" << selectedPatches_ << endl;
        Pout<< "radiationModels::viewFactor::initialise() "
            << "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = nLocalCoarseFaces_;
    reduce(totalNCoarseFaces_, sumOp<label>());

    if (debug && Pstream::master())
    {
        InfoInFunction
            << "Total number of clusters : " << totalNCoarseFaces_ << endl;
    }

    labelListIOList subMap
    (
        IOobject
        (
            "subMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelListIOList constructMap
    (
        IOobject
        (
            "constructMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    IOList<label> consMapDim
    (
        IOobject
        (
            "constructMapDim",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    map_.reset
    (
        new distributionMap
        (
            consMapDim[0],
            move(subMap),
            move(constructMap)
        )
    );

    scalarListIOList FmyProc
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelListIOList globalFaceFaces
    (
        IOobject
        (
            "globalFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    List<labelListList> globalFaceFacesProc(Pstream::nProcs());
    globalFaceFacesProc[Pstream::myProcNo()] = globalFaceFaces;
    Pstream::gatherList(globalFaceFacesProc);

    List<scalarListList> F(Pstream::nProcs());
    F[Pstream::myProcNo()] = FmyProc;
    Pstream::gatherList(F);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    if (Pstream::master())
    {
        Fmatrix_.reset
        (
            new scalarSquareMatrix(totalNCoarseFaces_, 0.0)
        );

        if (debug)
        {
            InfoInFunction
                << "Insert elements in the matrix..." << endl;
        }

        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            insertMatrixElements
            (
                globalNumbering,
                proci,
                globalFaceFacesProc[proci],
                F[proci],
                Fmatrix_()
            );
        }


        bool smoothing = readBool(coeffs_.lookup("smoothing"));
        if (smoothing)
        {
            if (debug)
            {
                InfoInFunction
                    << "Smoothing the matrix..." << endl;
            }

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                scalar sumF = 0.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    sumF += Fmatrix_()(i, j);
                }

                const scalar delta = sumF - 1.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    Fmatrix_()(i, j) *= (1.0 - delta/(sumF + 0.001));
                }
            }
        }

        constEmissivity_ = readBool(coeffs_.lookup("constantEmissivity"));
        if (constEmissivity_)
        {
            CLU_.reset
            (
                new scalarSquareMatrix(totalNCoarseFaces_, 0.0)
            );

            pivotIndices_.setSize(CLU_().m());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::viewFactor::viewFactor(const volScalarField& T)
:
    radiationModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


Foam::radiationModels::viewFactor::viewFactor
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::viewFactor::~viewFactor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::viewFactor::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiationModels::viewFactor::insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label proci,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& Fmatrix
)
{
    forAll(viewFactors, facei)
    {
        const scalarList& vf = viewFactors[facei];
        const labelList& globalFaces = globalFaceFaces[facei];

        label globalI = globalNumbering.toGlobal(proci, facei);
        forAll(globalFaces, i)
        {
            Fmatrix[globalI][globalFaces[i]] = vf[i];
        }
    }
}


void Foam::radiationModels::viewFactor::calculate()
{
    // Store previous iteration
    qr_.storePrevIter();

    scalarField compactCoarseT4(map_->constructSize(), 0.0);
    scalarField compactCoarseE(map_->constructSize(), 0.0);
    scalarField compactCoarseHo(map_->constructSize(), 0.0);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Fill local averaged(T), emissivity(E) and external heatFlux(Ho)
    DynamicList<scalar> localCoarseT4ave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseEave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    forAll(selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& Tp = T_.boundaryField()[patchID];
        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        fvPatchScalarField& qrPatch = qrBf[patchID];

        greyDiffusiveViewFactorFixedValueFvPatchScalarField& qrp =
            refCast
            <
                greyDiffusiveViewFactorFixedValueFvPatchScalarField
            >(qrPatch);

        const scalarList eb = qrp.emissivity();

        const scalarList& Hoi = qrp.qro();

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID];
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID];

        scalarList T4ave(pp.size(), 0.0);
        scalarList Eave(pp.size(), 0.0);
        scalarList Hoiave(pp.size(), 0.0);

        if (pp.size() > 0)
        {
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom) + 1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            forAll(coarseToFine, coarseI)
            {
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                UIndirectList<scalar> fineSf
                (
                    sf,
                    fineFaces
                );

                const scalar area = sum(fineSf());

                // Temperature, emissivity and external flux area weighting
                forAll(fineFaces, j)
                {
                    label facei = fineFaces[j];
                    T4ave[coarseI] += (pow4(Tp[facei])*sf[facei])/area;
                    Eave[coarseI] += (eb[facei]*sf[facei])/area;
                    Hoiave[coarseI] += (Hoi[facei]*sf[facei])/area;
                }
            }
        }

        localCoarseT4ave.append(T4ave);
        localCoarseEave.append(Eave);
        localCoarseHoave.append(Hoiave);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseT4, nLocalCoarseFaces_) = localCoarseT4ave;
    SubList<scalar>(compactCoarseE, nLocalCoarseFaces_) = localCoarseEave;
    SubList<scalar>(compactCoarseHo, nLocalCoarseFaces_) = localCoarseHoave;

    // Distribute data
    map_->distribute(compactCoarseT4);
    map_->distribute(compactCoarseE);
    map_->distribute(compactCoarseHo);

    // Distribute local global ID
    labelList compactGlobalIds(map_->constructSize(), 0.0);

    labelList localGlobalIds(nLocalCoarseFaces_);

    for(label k = 0; k < nLocalCoarseFaces_; k++)
    {
        localGlobalIds[k] = globalNumbering.toGlobal(Pstream::myProcNo(), k);
    }

    SubList<label>
    (
        compactGlobalIds,
        nLocalCoarseFaces_
    ) = localGlobalIds;

    map_->distribute(compactGlobalIds);

    // Create global size vectors
    scalarField T4(totalNCoarseFaces_, 0.0);
    scalarField E(totalNCoarseFaces_, 0.0);
    scalarField qrExt(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll(compactCoarseT4, i)
    {
        T4[compactGlobalIds[i]] = compactCoarseT4[i];
        E[compactGlobalIds[i]] = compactCoarseE[i];
        qrExt[compactGlobalIds[i]] = compactCoarseHo[i];
    }

    Pstream::listCombineGather(T4, maxEqOp<scalar>());
    Pstream::listCombineGather(E, maxEqOp<scalar>());
    Pstream::listCombineGather(qrExt, maxEqOp<scalar>());

    Pstream::listCombineScatter(T4);
    Pstream::listCombineScatter(E);
    Pstream::listCombineScatter(qrExt);

    // Net radiation
    scalarField q(totalNCoarseFaces_, 0.0);

    if (Pstream::master())
    {
        // Variable emissivity
        if (!constEmissivity_)
        {
            scalarSquareMatrix C(totalNCoarseFaces_, 0.0);

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    const scalar invEj = 1.0/E[j];
                    const scalar sigmaT4 = physicoChemical::sigma.value()*T4[j];

                    if (i==j)
                    {
                        C(i, j) = invEj - (invEj - 1.0)*Fmatrix_()(i, j);
                        q[i] += (Fmatrix_()(i, j) - 1.0)*sigmaT4 - qrExt[j];
                    }
                    else
                    {
                        C(i, j) = (1.0 - invEj)*Fmatrix_()(i, j);
                        q[i] += Fmatrix_()(i, j)*sigmaT4;
                    }

                }
            }

            Info<< "\nSolving view factor equations..." << endl;

            // Negative coming into the fluid
            LUsolve(C, q);
        }
        else // Constant emissivity
        {
            // Initial iter calculates CLU and chaches it
            if (iterCounter_ == 0)
            {
                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    {
                        const scalar invEj = 1.0/E[j];
                        if (i==j)
                        {
                            CLU_()(i, j) = invEj-(invEj-1.0)*Fmatrix_()(i, j);
                        }
                        else
                        {
                            CLU_()(i, j) = (1.0 - invEj)*Fmatrix_()(i, j);
                        }
                    }
                }

                if (debug)
                {
                    InfoInFunction
                        << "\nDecomposing C matrix..." << endl;
                }

                LUDecompose(CLU_(), pivotIndices_);
            }

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    const scalar sigmaT4 =
                        constant::physicoChemical::sigma.value()*T4[j];

                    if (i==j)
                    {
                        q[i] += (Fmatrix_()(i, j) - 1.0)*sigmaT4  - qrExt[j];
                    }
                    else
                    {
                        q[i] += Fmatrix_()(i, j)*sigmaT4;
                    }
                }
            }

            if (debug)
            {
                InfoInFunction
                    << "\nLU Back substitute C matrix.." << endl;
            }

            LUBacksubstitute(CLU_(), pivotIndices_, q);
            iterCounter_ ++;
        }
    }

    // Scatter q and fill qr
    Pstream::listCombineScatter(q);
    Pstream::listCombineGather(q, maxEqOp<scalar>());

    label globCoarseId = 0;
    forAll(selectedPatches_, i)
    {
        const label patchID = selectedPatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        if (pp.size() > 0)
        {
            scalarField& qrp = qrBf[patchID];
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom)+1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            const labelList& coarsePatchFace =
                coarseMesh_.patchFaceMap()[patchID];

            forAll(coarseToFine, coarseI)
            {
                label globalCoarse =
                    globalNumbering.toGlobal(Pstream::myProcNo(), globCoarseId);
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                forAll(fineFaces, k)
                {
                    label facei = fineFaces[k];

                    qrp[facei] = q[globalCoarse];
                }
                globCoarseId ++;
            }
        }
    }

    if (debug)
    {
        forAll(qrBf, patchID)
        {
            const scalarField& qrp = qrBf[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            const scalar heatFlux = gSum(qrp*magSf);

            InfoInFunction
                << "Total heat transfer rate at patch: "
                << patchID << " "
                << heatFlux << endl;
        }
    }

    // Relax qr if necessary
    qr_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::viewFactor::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        mesh_,
        dimensionedScalar
        (
            dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
            0
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::viewFactor::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}

// ************************************************************************* //
