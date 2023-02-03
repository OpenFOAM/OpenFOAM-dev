/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "BlendedInterfacialModel.H"
#include "phaseSystem.H"
#include "dispersedDisplacedPhaseInterface.H"
#include "segregatedDisplacedPhaseInterface.H"
#include "fixedValueFvsPatchFields.H"
#include "surfaceInterpolate.H"
#include "zeroDimensionalFvMesh.H"
#include "mathematicalConstants.H"
#include "writeFile.H"
#include "triFace.H"
#include "noSetWriter.H"
#include "noSurfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendedInterfacialModel
{

template<class GeoField>
inline tmp<GeoField> interpolate(tmp<volScalarField> f);

template<>
inline tmp<Foam::volScalarField> interpolate(tmp<volScalarField> f)
{
    return f;
}

template<>
inline tmp<Foam::surfaceScalarField> interpolate(tmp<volScalarField> f)
{
    return fvc::interpolate(f);
}

template<class ModelType>
inline bool valid(const PtrList<ModelType>& l)
{
    forAll(l, i)
    {
        if (l.set(i))
        {
            return true;
        }
    }
    return false;
}

} // End namespace blendedInterfacialModel
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ModelType>
void Foam::BlendedInterfacialModel<ModelType>::check() const
{
    // Only generate warnings once per timestep
    const label timeIndex = interface_.mesh().time().timeIndex();
    if (checkTimeIndex_ == timeIndex) return;
    checkTimeIndex_ = timeIndex;

    const phaseModel& phase1 = interface_.phase1();
    const phaseModel& phase2 = interface_.phase2();

    const bool can1In2 = blending_->canBeContinuous(1);
    const bool can2In1 = blending_->canBeContinuous(0);
    const bool canS = blending_->canSegregate();

    // Warnings associated with redundant model specification

    if
    (
       !can1In2
     && (
            model1DispersedIn2_.valid()
         || blendedInterfacialModel::valid(models1DispersedIn2Displaced_)
        )
    )
    {
        WarningInFunction
            << "A " << ModelType::typeName << " was provided for "
            << dispersedPhaseInterface(phase1, phase2).name()
            << " but the associated blending does not permit " << phase2.name()
            << " to be continuous so this model will not be used" << endl;
    }

    if
    (
       !can2In1
     && (
            model2DispersedIn1_.valid()
         || blendedInterfacialModel::valid(models2DispersedIn1Displaced_)
        )
    )
    {
        WarningInFunction
            << "A " << ModelType::typeName << " was provided for "
            << dispersedPhaseInterface(phase2, phase1).name()
            << " but the associated blending does not permit " << phase1.name()
            << " to be continuous so this model will not be used" << endl;
    }

    if
    (
       !canS
     && (
            model1SegregatedWith2_.valid()
         || blendedInterfacialModel::valid(models1SegregatedWith2Displaced_)
        )
    )
    {
        WarningInFunction
            << "A " << ModelType::typeName << " was provided for "
            << segregatedPhaseInterface(phase1, phase2).name()
            << " but the associated blending does not permit segregation"
            << " so this model will not be used" << endl;
    }

    if
    (
        modelGeneral_.valid()
     && (can1In2 || can2In1 || canS)
     && (!can1In2 || model1DispersedIn2_.valid())
     && (!can2In1 || model2DispersedIn1_.valid())
     && (!canS || model1SegregatedWith2_.valid())
    )
    {
        WarningInFunction
            << "A " << ModelType::typeName << " was provided for "
            << phaseInterface(phase1, phase2).name()
            << " but other displaced and/or segregated models apply"
            << " across the entire phase fraction space so this model"
            << " will not be used" << endl;
    }

    forAll(interface_.fluid().phases(), phasei)
    {
        const phaseModel& phaseD = interface_.fluid().phases()[phasei];

        if
        (
            modelsGeneralDisplaced_.set(phasei)
         && (can1In2 || can2In1 || canS)
         && (!can1In2 || models1DispersedIn2Displaced_.set(phasei))
         && (!can2In1 || models2DispersedIn1Displaced_.set(phasei))
         && (!canS || models1SegregatedWith2Displaced_.set(phasei))
        )
        {
            WarningInFunction
                << "A " << ModelType::typeName << " was provided for "
                << displacedPhaseInterface(phase1, phase2, phaseD).name()
                << " but other displaced and/or segregated models apply"
                << " across the entire phase fraction space so this model"
                << " will not be used" << endl;
        }
    }

    // Warnings associated with gaps in the blending space

    if (can1In2 && !modelGeneral_.valid() && !model1DispersedIn2_.valid())
    {
        WarningInFunction
            << "Blending for " << ModelType::typeName << "s permits "
            << phase2.name() << " to become continuous, but no model was "
            << "provided for " << dispersedPhaseInterface(phase1, phase2).name()
            << ". Consider adding a model for this configuration (or for "
            << phaseInterface(phase1, phase2).name() << "), or if no model is "
            << "needed then add a \"none\" model to suppress this warning."
            << endl;
    }

    if (can2In1 && !modelGeneral_.valid() && !model2DispersedIn1_.valid())
    {
        WarningInFunction
            << "Blending for " << ModelType::typeName << "s permits "
            << phase1.name() << " to become continuous, but no model was "
            << "provided for " << dispersedPhaseInterface(phase2, phase1).name()
            << ". Consider adding a model for this configuration (or for "
            << phaseInterface(phase1, phase2).name() << "), or if no model is "
            << "needed then add a \"none\" model to suppress this warning."
            << endl;
    }

    if (canS && !modelGeneral_.valid() && !model1SegregatedWith2_.valid())
    {
        WarningInFunction
            << "Blending for " << ModelType::typeName << "s permits "
            << "segregation but no model was provided for "
            << segregatedPhaseInterface(phase2, phase1).name()
            << ". Consider adding a model for this configuration (or for "
            << phaseInterface(phase1, phase2).name() << "), or if no model is "
            << "needed then add a \"none\" model to suppress this warning."
            << endl;
    }
}


template<class ModelType>
template<template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::calculateBlendingCoeffs
(
    const UPtrList<const volScalarField>& alphas,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& fG,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f1D2,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f2D1,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& fS,
    PtrList<GeometricField<scalar, PatchField, GeoMesh>>& fGD,
    PtrList<GeometricField<scalar, PatchField, GeoMesh>>& f1D2D,
    PtrList<GeometricField<scalar, PatchField, GeoMesh>>& f2D1D,
    PtrList<GeometricField<scalar, PatchField, GeoMesh>>& fSD,
    const bool subtract
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;

    // Create a constant field
    auto constant = [&](const scalar k)
    {
        return
            scalarGeoField::New
            (
                Foam::name(k),
                alphas.first().mesh(),
                dimensionedScalar(dimless, k)
            );
    };

    // Get the dispersed blending functions
    tmp<scalarGeoField> F1D2, F2D1;
    if
    (
        model1DispersedIn2_.valid()
     || model1SegregatedWith2_.valid()
     || blendedInterfacialModel::valid(models1DispersedIn2Displaced_)
     || blendedInterfacialModel::valid(models1SegregatedWith2Displaced_)
    )
    {
        F1D2 =
            blendedInterfacialModel::interpolate<scalarGeoField>
            (
                blending_->f1DispersedIn2(alphas)
            );
    }
    if
    (
        model2DispersedIn1_.valid()
     || model1SegregatedWith2_.valid()
     || blendedInterfacialModel::valid(models2DispersedIn1Displaced_)
     || blendedInterfacialModel::valid(models1SegregatedWith2Displaced_)
    )
    {
        F2D1 =
            blendedInterfacialModel::interpolate<scalarGeoField>
            (
                blending_->f2DispersedIn1(alphas)
            );
    }

    // Construct non-displaced coefficients
    {
        if (model1DispersedIn2_.valid())
        {
            f1D2 = F1D2().clone();
        }

        if (model2DispersedIn1_.valid())
        {
            f2D1 = F2D1().clone();
        }

        if (model1SegregatedWith2_.valid())
        {
            fS = constant(1);
            if (f1D2.valid()) { fS.ref() -= f1D2(); }
            if (f2D1.valid()) { fS.ref() -= f2D1(); }
        }

        if (modelGeneral_.valid())
        {
            fG = constant(1);
            if (f1D2.valid()) { fG.ref() -= f1D2(); }
            if (f2D1.valid()) { fG.ref() -= f2D1(); }
            if (fS.valid()) { fG.ref() -= fS(); }
        }
    }

    // Get the displaced blending function
    tmp<scalarGeoField> FD;
    if
    (
        blendedInterfacialModel::valid(modelsGeneralDisplaced_)
     || blendedInterfacialModel::valid(models1DispersedIn2Displaced_)
     || blendedInterfacialModel::valid(models2DispersedIn1Displaced_)
     || blendedInterfacialModel::valid(models1SegregatedWith2Displaced_)
    )
    {
        FD =
            blendedInterfacialModel::interpolate<scalarGeoField>
            (
                blending_->fDisplaced(alphas)
            );
    }

    // Construct displaced coefficients
    tmp<scalarGeoField> fDSum;
    if
    (
        blendedInterfacialModel::valid(modelsGeneralDisplaced_)
     || blendedInterfacialModel::valid(models1DispersedIn2Displaced_)
     || blendedInterfacialModel::valid(models2DispersedIn1Displaced_)
     || blendedInterfacialModel::valid(models1SegregatedWith2Displaced_)
    )
    {
        fDSum = constant(0);

        forAll(alphas, phasei)
        {
            const phaseModel& phaseD = interface_.fluid().phases()[phasei];

            if (interface_.contains(phaseD)) continue;

            // Weight the contribution of a dispersed model for this phase
            // according to the phase fraction in the subset of phases that are
            // not part of this interface. This seems like a reasonable
            // assumption if a stochastic viewpoint is taken. However, it is
            // possible that other forms may be desired in which case this
            // could be abstracted and controlled by the blending method.
            tmp<scalarGeoField> alpha;
            if
            (
                models1DispersedIn2Displaced_.set(phasei)
             || models2DispersedIn1Displaced_.set(phasei)
             || models1SegregatedWith2Displaced_.set(phasei)
             || modelsGeneralDisplaced_.set(phasei)
            )
            {
                alpha =
                    blendedInterfacialModel::interpolate<scalarGeoField>
                    (
                        alphas[phasei]
                       /max
                        (
                            1
                          - alphas[interface_.phase1().index()]
                          - alphas[interface_.phase2().index()],
                            phaseD.residualAlpha()
                        )
                    );
            }

            if (models1DispersedIn2Displaced_.set(phasei))
            {
                f1D2D.set(phasei, alpha()*FD()*F1D2());
                fDSum.ref() += f1D2D[phasei];
            }

            if (models2DispersedIn1Displaced_.set(phasei))
            {
                f2D1D.set(phasei, alpha()*FD()*F2D1());
                fDSum.ref() += f2D1D[phasei];
            }

            if (models1SegregatedWith2Displaced_.set(phasei))
            {
                fSD.set(phasei, alpha()*FD());
                if (f1D2D.set(phasei)) fSD[phasei] -= f1D2D[phasei];
                if (f2D1D.set(phasei)) fSD[phasei] -= f2D1D[phasei];
                fDSum.ref() += fSD[phasei];
            }

            if (modelsGeneralDisplaced_.set(phasei))
            {
                fGD.set(phasei, alpha()*FD());
                if (f1D2D.set(phasei)) fGD[phasei] -= f1D2D[phasei];
                if (f2D1D.set(phasei)) fGD[phasei] -= f2D1D[phasei];
                if (fSD.set(phasei)) fGD[phasei] -= fSD[phasei];
                fDSum.ref() += fGD[phasei];
            }
        }
    }

    // Remove the displaced part from the non-displaced coefficients
    if (fDSum.valid())
    {
        if (f1D2.valid()) f1D2.ref() *= 1 - fDSum();
        if (f2D1.valid()) f2D1.ref() *= 1 - fDSum();
        if (fS.valid()) fS.ref() *= 1 - fDSum();
        if (fG.valid()) fG.ref() *= 1 - fDSum();
    }

    // Flip the sign of the 2 dispersed in 1 models if necessary
    if (subtract)
    {
        auto signedError = [this](const phaseInterface& interface)
        {
            FatalErrorInFunction
                << "A signed quantity was evaluated from the blended "
                << ModelType::typeName << " for " << interface_.name()
                << " but a model was provided for " << interface.name()
                << ". Signed quantities are only possible to evaluate for"
                << " dispersed configurations" << exit(FatalError);
        };

        const phaseModel& phase1 = interface_.phase1();
        const phaseModel& phase2 = interface_.phase2();

        if (fG.valid())
        {
            signedError(phaseInterface(phase1, phase2));
        }
        if (f1D2.valid())
        {
            // Do nothing
        }
        if (f2D1.valid())
        {
            f2D1.ref() *= -1;
        }
        if (fS.valid())
        {
            signedError(segregatedPhaseInterface(phase1, phase2));
        }

        forAll(alphas, phasei)
        {
            const phaseModel& phaseD = interface_.fluid().phases()[phasei];

            if (fGD.set(phasei))
            {
                signedError
                (
                    displacedPhaseInterface(phase1, phase2, phaseD)
                );
            }
            if (f1D2D.set(phasei))
            {
                // Do nothing
            }
            if (f2D1D.set(phasei))
            {
                f2D1.ref() *= -1;
            }
            if (fSD.set(phasei))
            {
                signedError
                (
                    segregatedDisplacedPhaseInterface(phase1, phase2, phaseD)
                );
            }
        }
    }
}


template<class ModelType>
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::correctFixedFluxBCs
(
    GeometricField<Type, PatchField, GeoMesh>& field
) const
{
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    typename typeGeoField::Boundary& fieldBf = field.boundaryFieldRef();

    forAll(fieldBf, patchi)
    {
        if
        (
            (
                !interface_.phase1().stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    interface_.phase1().phi()().boundaryField()[patchi]
                )
            )
         || (
                !interface_.phase2().stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    interface_.phase2().phi()().boundaryField()[patchi]
                )
            )
        )
        {
            fieldBf[patchi] = Zero;
        }
    }
}


template<class ModelType>
void Foam::BlendedInterfacialModel<ModelType>::postProcessBlendingCoefficients
(
    const word& format
) const
{
    using namespace constant::mathematical;

    const phaseSystem& fluid = interface_.fluid();
    const label nPhases = fluid.phases().size();

    // Construct geometry and phase fraction values
    pointField points;
    faceList faces;
    wordList fieldNames(nPhases);
    PtrList<scalarField> fields(nPhases);
    if (nPhases == 1)
    {
        // Single value
        fieldNames[0] = fluid.phases()[0].volScalarField::name();
        fields.set(0, new scalarField(1, 1));
    }
    else if (nPhases == 2)
    {
        // Single axis, using the first phase as the x-coordinate
        static const label nDivisions = 128;
        const scalarField divisionis(scalarList(identityMap(nDivisions)));
        fieldNames[0] = fluid.phases()[0].volScalarField::name();
        fieldNames[1] = fluid.phases()[1].volScalarField::name();
        fields.set(0, new scalarField(divisionis/(nDivisions - 1)));
        fields.set(1, new scalarField(1 - fields[0]));
    }
    else
    {
        // Polygon with as many vertices as there are phases. Each phase
        // fraction equals one at a unique vertex, they vary linearly between
        // vertices, and vary smoothly in the interior of the polygon. The sum
        // of phase fractions is always one throughout the polygon (partition
        // of unity). This is done using Wachspress coordinates. This does not
        // cover the entire N-dimensional phase fraction space for N >= 4 (it's
        // hard to imagine how that could be visualised) but it provides enough
        // to give a good indication of what is going on.

        // Create the nodes of the blending space polygon
        List<point> phaseNodes(nPhases);
        forAll(phaseNodes, phasei)
        {
            const scalar theta = 2*pi*phasei/nPhases;
            phaseNodes[phasei] = point(cos(theta), sin(theta), 0);
        }

        // Create points within the polygon
        static const label nDivisions = 32;
        points.append(point::zero);
        for (label divi = 0; divi < nDivisions; ++ divi)
        {
            const scalar s = scalar(divi + 1)/nDivisions;

            forAll(phaseNodes, phasei)
            {
                for (label i = 0; i < divi + 1; ++ i)
                {
                    const scalar t = scalar(i)/(divi + 1);

                    points.append
                    (
                        s*(1 - t)*phaseNodes[phasei]
                      + s*t*phaseNodes[(phasei + 1) % nPhases]
                    );
                }
            }
        }

        // Create triangles within the polygon
        forAll(phaseNodes, phasei)
        {
            faces.append(triFace(0, phasei + 1, (phasei + 1) % nPhases + 1));
        }
        for (label divi = 1; divi < nDivisions; ++ divi)
        {
            const label pointi0 = nPhases*(divi - 1)*divi/2 + 1;
            const label pointi1 = nPhases*divi*(divi + 1)/2 + 1;

            forAll(phaseNodes, phasei)
            {
                for (label i = 0; i < divi + 1; ++ i)
                {
                    const label pi00 =
                        pointi0
                      + ((phasei*divi + i) % (nPhases*divi));
                    const label pi01 =
                        pointi0
                      + ((phasei*divi + i + 1) % (nPhases*divi));
                    const label pi10 =
                        pointi1
                      + ((phasei*(divi + 1) + i) % (nPhases*(divi + 1)));
                    const label pi11 =
                        pointi1
                      + ((phasei*(divi + 1) + i + 1) % (nPhases*(divi + 1)));

                    faces.append(triFace({pi00, pi10, pi11}));
                    if (i < divi) faces.append(triFace({pi00, pi11, pi01}));
                }
            }
        }

        // Create phase fraction fields
        forAll(fluid.phases(), phasei)
        {
            fieldNames[phasei] = fluid.phases()[phasei].volScalarField::name();
            fields.set(phasei, new scalarField(points.size(), 0));

            const label phasei0 = (phasei + nPhases - 1) % nPhases;
            const label phasei1 = (phasei + 1) % nPhases;

            const point& node0 = phaseNodes[phasei0];
            const point& node = phaseNodes[phasei];
            const point& node1 = phaseNodes[phasei1];

            forAll(points, i)
            {
                const scalar A = mag((node - node0) ^ (node1 - node0));
                const scalar A0 = mag((node - node0) ^ (points[i] - node0));
                const scalar A1 = mag((node1 - node) ^ (points[i] - node));

                if (A0 < rootSmall)
                {
                    fields[phasei][i] =
                        great*mag(points[i] - node0)/mag(node - node0);
                }
                else if (A1 < rootSmall)
                {
                    fields[phasei][i] =
                        great*mag(points[i] - node1)/mag(node - node1);
                }
                else
                {
                    fields[phasei][i] = A/A0/A1;
                }
            }
        }
        forAll(points, i)
        {
            scalar s = 0;
            forAll(fluid.phases(), phasei)
            {
                s += fields[phasei][i];
            }
            forAll(fluid.phases(), phasei)
            {
                fields[phasei][i] /= s;
            }
        }
    }

    // Add the model coefficient fields
    {
        // Create alpha fields on a zero-dimensional one-cell mesh
        const fvMesh mesh(zeroDimensionalFvMesh(fluid.mesh()));
        PtrList<volScalarField> alphas(nPhases);
        forAll(fluid.phases(), phasei)
        {
            alphas.set
            (
                phasei,
                new volScalarField
                (
                    IOobject
                    (
                        fluid.phases()[phasei].volScalarField::name(),
                        mesh.time().name(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar(dimless, fields[phasei][0])
                )
            );
        }

        // Construct blending coefficient fields
        {
            tmp<volScalarField> fG, f1D2, f2D1, fS;
            PtrList<volScalarField> fGD(nPhases);
            PtrList<volScalarField> f1D2D(nPhases);
            PtrList<volScalarField> f2D1D(nPhases);
            PtrList<volScalarField> fSD(nPhases);
            calculateBlendingCoeffs
            (
                alphas.convert<const volScalarField>(),
                fG, f1D2, f2D1, fS,
                fGD, f1D2D, f2D1D, fSD,
                false
            );

            const phaseModel& phase1 = interface_.phase1();
            const phaseModel& phase2 = interface_.phase2();

            auto addField = [&](const phaseInterface& interface)
            {
                fieldNames.append(interface.name());
                fields.append(new scalarField(fields.first().size()));
            };

            if (fG.valid())
            {
                addField(phaseInterface(phase1, phase2));
            }
            if (f1D2.valid())
            {
                addField(dispersedPhaseInterface(phase1, phase2));
            }
            if (f2D1.valid())
            {
                addField(dispersedPhaseInterface(phase2, phase1));
            }
            if (fS.valid())
            {
                addField(segregatedPhaseInterface(phase2, phase1));
            }

            forAll(fluid.phases(), phasei)
            {
                const phaseModel& phaseD = fluid.phases()[phasei];

                if (fGD.set(phasei))
                {
                    addField(displacedPhaseInterface(phase1, phase2, phaseD));
                }
                if (f1D2D.set(phasei))
                {
                    addField
                    (
                        dispersedDisplacedPhaseInterface
                        (
                            phase1,
                            phase2,
                            phaseD
                        )
                    );
                }
                if (f2D1D.set(phasei))
                {
                    addField
                    (
                        dispersedDisplacedPhaseInterface
                        (
                            phase2,
                            phase1,
                            phaseD
                        )
                    );
                }
                if (fSD.set(phasei))
                {
                    addField
                    (
                        segregatedDisplacedPhaseInterface
                        (
                            phase1,
                            phase2,
                            phaseD
                        )
                    );
                }
            }
        }

        // Populate blending coefficient fields
        forAll(fields.first(), i)
        {
            forAll(fluid.phases(), phasei)
            {
                alphas[phasei] = fields[phasei][i];
            }

            tmp<volScalarField> fG, f1D2, f2D1, fS;
            PtrList<volScalarField> fGD(nPhases);
            PtrList<volScalarField> f1D2D(nPhases);
            PtrList<volScalarField> f2D1D(nPhases);
            PtrList<volScalarField> fSD(nPhases);
            calculateBlendingCoeffs
            (
                alphas.convert<const volScalarField>(),
                fG, f1D2, f2D1, fS,
                fGD, f1D2D, f2D1D, fSD,
                false
            );

            label fieldi = nPhases;

            if (fG.valid())
            {
                fields[fieldi ++][i] = fG()[0];
            }
            if (f1D2.valid())
            {
                fields[fieldi ++][i] = f1D2()[0];
            }
            if (f2D1.valid())
            {
                fields[fieldi ++][i] = f2D1()[0];
            }
            if (fS.valid())
            {
                fields[fieldi ++][i] = fS()[0];
            }

            forAll(fluid.phases(), phasei)
            {
                if (fGD.set(phasei))
                {
                    fields[fieldi ++][i] = fGD[phasei][0];
                }
                if (f1D2D.set(phasei))
                {
                    fields[fieldi ++][i] = f1D2D[phasei][0];
                }
                if (f2D1D.set(phasei))
                {
                    fields[fieldi ++][i] = f2D1D[phasei][0];
                }
                if (fSD.set(phasei))
                {
                    fields[fieldi ++][i] = fSD[phasei][0];
                }
            }
        }
    }

    // Write
    const fileName path =
        fluid.mesh().time().globalPath()
       /functionObjects::writeFile::outputPrefix
       /ModelType::typeName;
    Info<< "Writing blending coefficients to " << path/interface_.name()
        << endl;
    if (nPhases <= 2)
    {
        // Strip out the first field and shuffle everything else up
        const autoPtr<scalarField> field0 = fields.set(0, nullptr);
        const word field0Name = fieldNames[0];
        for (label fieldi = 1; fieldi < fields.size(); ++ fieldi)
        {
            fields.set(fieldi - 1, fields.set(fieldi, nullptr).ptr());
            fieldNames[fieldi - 1] = fieldNames[fieldi];
        }
        fields.resize(fields.size() - 1);
        fieldNames.resize(fieldNames.size() - 1);

        setWriter::New
        (
            format,
            IOstream::ASCII,
            IOstream::UNCOMPRESSED
        )->write
        (
            path,
            interface_.name(),
            coordSet(true, field0Name, field0),
            fieldNames,
            fields
        );
    }
    else
    {
        surfaceWriter::New
        (
            format,
            IOstream::ASCII,
            IOstream::UNCOMPRESSED
        )->write
        (
            path,
            interface_.name(),
            points,
            faces,
            true,
            fieldNames,
            fields
        );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ModelType>
template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class ... Args
>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    tmp<GeometricField<Type, PatchField, GeoMesh>>
    (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const bool subtract,
    Args ... args
) const
{
    check();

    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    // Get the blending coefficients
    const label nPhases = interface_.fluid().phases().size();
    tmp<scalarGeoField> fG, f1D2, f2D1, fS;
    PtrList<scalarGeoField> fGD(nPhases);
    PtrList<scalarGeoField> f1D2D(nPhases);
    PtrList<scalarGeoField> f2D1D(nPhases);
    PtrList<scalarGeoField> fSD(nPhases);
    calculateBlendingCoeffs
    (
        interface_.fluid().phases()
       .PtrList<phaseModel>::convert<const volScalarField>(),
        fG, f1D2, f2D1, fS,
        fGD, f1D2D, f2D1D, fSD,
        subtract
    );

    // Construct the result
    tmp<typeGeoField> x =
        typeGeoField::New
        (
            ModelType::typeName + ":"
          + IOobject::groupName(name, interface_.name()),
            interface_.mesh(),
            dimensioned<Type>(dims, Zero)
        );

    // Add the model contributions to the result
    if (modelGeneral_.valid())
    {
        x.ref() += fG*(modelGeneral_().*method)(args ...);
    }
    if (model1DispersedIn2_.valid())
    {
        x.ref() += f1D2*(model1DispersedIn2_().*method)(args ...);
    }
    if (model2DispersedIn1_.valid())
    {
        x.ref() += f2D1*(model2DispersedIn1_().*method)(args ...);
    }
    if (model1SegregatedWith2_.valid())
    {
        x.ref() += fS*(model1SegregatedWith2_().*method)(args ...);
    }
    forAll(interface_.fluid().phases(), phasei)
    {
        if (modelsGeneralDisplaced_.set(phasei))
        {
            x.ref() +=
                fGD[phasei]
               *(modelsGeneralDisplaced_[phasei].*method)(args ...);
        }
        if (models1DispersedIn2Displaced_.set(phasei))
        {
            x.ref() +=
                f1D2D[phasei]
               *(models1DispersedIn2Displaced_[phasei].*method)(args ...);
        }
        if (models2DispersedIn1Displaced_.set(phasei))
        {
            x.ref() +=
                f2D1D[phasei]
               *(models2DispersedIn1Displaced_[phasei].*method)(args ...);
        }
        if (models1SegregatedWith2Displaced_.set(phasei))
        {
            x.ref() +=
                fSD[phasei]
               *(models1SegregatedWith2Displaced_[phasei].*method)(args ...);
        }
    }

    // Correct boundary conditions if necessary
    if (ModelType::correctFixedFluxBCs)
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class ModelType>
template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class ... Args
>
Foam::HashPtrTable<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    HashPtrTable<GeometricField<Type, PatchField, GeoMesh>>
    (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const bool subtract,
    Args ... args
) const
{
    check();

    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    // Get the blending coefficients
    const label nPhases = interface_.fluid().phases().size();
    tmp<scalarGeoField> fG, f1D2, f2D1, fS;
    PtrList<scalarGeoField> fGD(nPhases);
    PtrList<scalarGeoField> f1D2D(nPhases);
    PtrList<scalarGeoField> f2D1D(nPhases);
    PtrList<scalarGeoField> fSD(nPhases);
    calculateBlendingCoeffs
    (
        interface_.fluid().phases()
       .PtrList<phaseModel>::convert<const volScalarField>(),
        fG, f1D2, f2D1, fS,
        fGD, f1D2D, f2D1D, fSD,
        subtract
    );

    // Construct the result
    HashPtrTable<typeGeoField> xs;

    // Add the model contributions to the result
    auto addToXs = [&]
    (
        const scalarGeoField& f,
        const HashPtrTable<typeGeoField>& dxs
    )
    {
        forAllConstIter(typename HashPtrTable<typeGeoField>, dxs, dxIter)
        {
            if (xs.found(dxIter.key()))
            {
                *xs[dxIter.key()] += f**dxIter();
            }
            else
            {
                xs.insert
                (
                    dxIter.key(),
                    typeGeoField::New
                    (
                        ModelType::typeName + ':'
                      + IOobject::groupName
                        (
                            IOobject::groupName(name, dxIter.key()),
                            interface_.name()
                        ),
                        f**dxIter()
                    ).ptr()
                );
            }
        }
    };
    if (modelGeneral_.valid())
    {
        addToXs(fG, (modelGeneral_().*method)(args ...));
    }
    if (model1DispersedIn2_.valid())
    {
        addToXs(f1D2, (model1DispersedIn2_().*method)(args ...));
    }
    if (model2DispersedIn1_.valid())
    {
        addToXs(f2D1, (model2DispersedIn1_().*method)(args ...));
    }
    if (model1SegregatedWith2_.valid())
    {
        addToXs(fS, (model1SegregatedWith2_().*method)(args ...));
    }
    forAll(interface_.fluid().phases(), phasei)
    {
        if (modelsGeneralDisplaced_.set(phasei))
        {
            addToXs
            (
                fGD[phasei],
                (modelsGeneralDisplaced_[phasei].*method)(args ...)
            );
        }
        if (models1DispersedIn2Displaced_.set(phasei))
        {
            addToXs
            (
                f1D2D[phasei],
                (models1DispersedIn2Displaced_[phasei].*method)(args ...)
            );
        }
        if (models2DispersedIn1Displaced_.set(phasei))
        {
            addToXs
            (
                f2D1D[phasei],
                (models2DispersedIn1Displaced_[phasei].*method)(args ...)
            );
        }
        if (models1SegregatedWith2Displaced_.set(phasei))
        {
            addToXs
            (
                fSD[phasei],
                (models1SegregatedWith2Displaced_[phasei].*method)(args ...)
            );
        }
    }

    // Correct boundary conditions if necessary
    if (ModelType::correctFixedFluxBCs)
    {
        forAllIter(typename HashPtrTable<typeGeoField>, xs, xIter)
        {
            correctFixedFluxBCs(*xIter());
        }
    }

    return xs;
}


template<class ModelType>
template<class ... Args>
bool Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    bool (ModelType::*method)(Args ...) const,
    Args ... args
) const
{
    check();

    bool result = false;

    if (modelGeneral_.valid())
    {
        result = result || (modelGeneral_().*method)(args ...);
    }
    if (model1DispersedIn2_.valid())
    {
        result = result || (model1DispersedIn2_().*method)(args ...);
    }
    if (model2DispersedIn1_.valid())
    {
        result = result || (model2DispersedIn1_().*method)(args ...);
    }
    if (model1SegregatedWith2_.valid())
    {
        result = result || (model1SegregatedWith2_().*method)(args ...);
    }

    forAll(interface_.fluid().phases(), phasei)
    {
        if (modelsGeneralDisplaced_.set(phasei))
        {
            result =
                result
             || (modelsGeneralDisplaced_[phasei].*method)(args ...);
        }
        if (models1DispersedIn2Displaced_.set(phasei))
        {
            result =
                result
             || (models1DispersedIn2Displaced_[phasei].*method)(args ...);
        }
        if (models2DispersedIn1Displaced_.set(phasei))
        {
            result =
                result
             || (models2DispersedIn1Displaced_[phasei].*method)(args ...);
        }
        if (models1SegregatedWith2Displaced_.set(phasei))
        {
            result =
                result
              || (models1SegregatedWith2Displaced_[phasei].*method)(args ...);
        }
    }

    return result;
}


template<class ModelType>
template<class ... Args>
Foam::hashedWordList Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    const hashedWordList& (ModelType::*method)(Args ...) const,
    Args ... args
) const
{
    check();

    wordList result;

    if (modelGeneral_.valid())
    {
        result.append((modelGeneral_().*method)(args ...));
    }
    if (model1DispersedIn2_.valid())
    {
        result.append((model1DispersedIn2_().*method)(args ...));
    }
    if (model2DispersedIn1_.valid())
    {
        result.append((model2DispersedIn1_().*method)(args ...));
    }
    if (model1SegregatedWith2_.valid())
    {
        result.append((model1SegregatedWith2_().*method)(args ...));
    }

    forAll(interface_.fluid().phases(), phasei)
    {
        if (modelsGeneralDisplaced_.set(phasei))
        {
            result.append
            (
                (modelsGeneralDisplaced_[phasei].*method)(args ...)
            );
        }
        if (models1DispersedIn2Displaced_.set(phasei))
        {
            result.append
            (
                (models1DispersedIn2Displaced_[phasei].*method)(args ...)
            );
        }
        if (models2DispersedIn1Displaced_.set(phasei))
        {
            result.append
            (
                (models2DispersedIn1Displaced_[phasei].*method)(args ...)
            );
        }
        if (models1SegregatedWith2Displaced_.set(phasei))
        {
            result.append
            (
                 (models1SegregatedWith2Displaced_[phasei].*method)(args ...)
            );
        }
    }

    return hashedWordList(move(result));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::BlendedInterfacialModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, interface.name()),
            interface.fluid().mesh().time().name(),
            interface.fluid().mesh()
        )
    ),
    interface_(interface),
    blending_(),
    modelGeneral_(),
    model1DispersedIn2_(),
    model2DispersedIn1_(),
    model1SegregatedWith2_(),
    modelsGeneralDisplaced_(interface.fluid().phases().size()),
    models1DispersedIn2Displaced_(interface.fluid().phases().size()),
    models2DispersedIn1Displaced_(interface.fluid().phases().size()),
    models1SegregatedWith2Displaced_(interface.fluid().phases().size()),
    checkTimeIndex_(-1)
{
    // Construct blending functions
    const dictionary& blendingDict = interface.fluid().subDict("blending");
    blending_ =
        blendingMethod::New
        (
            ModelType::typeName,
            blendingDict.found(interface.fluid().modelName<ModelType>())
          ? blendingDict.subDict(interface.fluid().modelName<ModelType>())
          : blendingDict.subDict("default"),
            interface
        );

    // Construct the models
    PtrList<phaseInterface> interfaces;
    PtrList<ModelType> models;
    interface.fluid().generateInterfacialModels
    <
        ModelType,
        dispersedDisplacedPhaseInterface,
        segregatedDisplacedPhaseInterface,
        displacedPhaseInterface,
        dispersedPhaseInterface,
        segregatedPhaseInterface,
        phaseInterface
    >
    (
        dict,
        interface,
        interfaces,
        models
    );

    // Unpack the interface and model lists to populate the models used for the
    // different parts of the blending space
    forAll(interfaces, i)
    {
        const phaseInterface& interface = interfaces[i];

        autoPtr<ModelType>* modelPtrPtr;
        PtrList<ModelType>* modelPtrsPtr;

        if (isA<dispersedPhaseInterface>(interface))
        {
            const phaseModel& dispersed =
                refCast<const dispersedPhaseInterface>(interface).dispersed();

            modelPtrPtr =
                interface_.index(dispersed) == 0
              ? &model1DispersedIn2_
              : &model2DispersedIn1_;
            modelPtrsPtr =
                interface_.index(dispersed) == 0
              ? &models1DispersedIn2Displaced_
              : &models2DispersedIn1Displaced_;
        }
        else if (isA<segregatedPhaseInterface>(interface))
        {
            modelPtrPtr = &model1SegregatedWith2_;
            modelPtrsPtr = &models1SegregatedWith2Displaced_;
        }
        else
        {
            modelPtrPtr = &modelGeneral_;
            modelPtrsPtr = &modelsGeneralDisplaced_;
        }

        if (!isA<displacedPhaseInterface>(interface))
        {
            *modelPtrPtr = models.set(i, nullptr);
        }
        else
        {
            const phaseModel& displacing =
                refCast<const displacedPhaseInterface>(interface).displacing();

            modelPtrsPtr->set(displacing.index(), models.set(i, nullptr));
        }
    }

    // Write out the blending space if needed
    if (blendingDict.found("format"))
    {
        const word format = blendingDict.lookup<word>("format");

        const label nPhases = interface_.fluid().phases().size();

        if
        (
            (nPhases <= 2 && format != noSetWriter::typeName)
         || (nPhases > 2 && format != noSurfaceWriter::typeName)
        )
        {
            postProcessBlendingCoefficients(format);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::~BlendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ModelType>
const Foam::phaseInterface&
Foam::BlendedInterfacialModel<ModelType>::interface() const
{
    return interface_;
}


template<class ModelType>
bool Foam::BlendedInterfacialModel<ModelType>::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
