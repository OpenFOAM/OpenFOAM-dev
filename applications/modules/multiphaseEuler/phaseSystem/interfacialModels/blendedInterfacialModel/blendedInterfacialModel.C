/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "blendedInterfacialModel.H"
#include "phaseSystem.H"
#include "surfaceInterpolate.H"
#include "dispersedDisplacedPhaseInterface.H"
#include "segregatedDisplacedPhaseInterface.H"
#include "boolListOps.H"
#include "triFace.H"
#include "zeroDimensionalFvMesh.H"
#include "noSetWriter.H"
#include "noSurfaceWriter.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<>
Foam::tmp<Foam::volScalarField>
Foam::blendedInterfacialModel::interpolate(tmp<volScalarField> f)
{
    return f;
}


template<>
Foam::tmp<Foam::surfaceScalarField>
Foam::blendedInterfacialModel::interpolate(tmp<volScalarField> f)
{
    return fvc::interpolate(f);
}


void Foam::blendedInterfacialModel::check
(
    const word& name,
    const boolSet& valid
) const
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
     && (valid.oneDispersedInTwo || any(valid.oneDispersedInTwoDisplaced))
    )
    {
        WarningInFunction
            << "A " << name << " was provided for "
            << dispersedPhaseInterface(phase1, phase2).name()
            << " but the associated blending does not permit " << phase2.name()
            << " to be continuous so this model will not be used" << endl;
    }

    if
    (
       !can2In1
     && (valid.twoDispersedInOne || any(valid.twoDispersedInOneDisplaced))
    )
    {
        WarningInFunction
            << "A " << name << " was provided for "
            << dispersedPhaseInterface(phase2, phase1).name()
            << " but the associated blending does not permit " << phase1.name()
            << " to be continuous so this model will not be used" << endl;
    }

    if
    (
       !canS
     && (valid.segregated || any(valid.segregatedDisplaced))
    )
    {
        WarningInFunction
            << "A " << name << " was provided for "
            << segregatedPhaseInterface(phase1, phase2).name()
            << " but the associated blending does not permit segregation"
            << " so this model will not be used" << endl;
    }

    if
    (
        valid.general
     && (can1In2 || can2In1 || canS)
     && (!can1In2 || valid.oneDispersedInTwo)
     && (!can2In1 || valid.twoDispersedInOne)
     && (!canS || valid.segregated)
    )
    {
        WarningInFunction
            << "A " << name << " was provided for "
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
            valid.generalDisplaced[phasei]
         && (can1In2 || can2In1 || canS)
         && (!can1In2 || valid.oneDispersedInTwoDisplaced[phasei])
         && (!can2In1 || valid.twoDispersedInOneDisplaced[phasei])
         && (!canS || valid.segregatedDisplaced[phasei])
        )
        {
            WarningInFunction
                << "A " << name << " was provided for "
                << displacedPhaseInterface(phase1, phase2, phaseD).name()
                << " but other displaced and/or segregated models apply"
                << " across the entire phase fraction space so this model"
                << " will not be used" << endl;
        }
    }

    // Warnings associated with gaps in the blending space

    if (!can1In2 && !can2In1 && !canS && !valid.general)
    {
        WarningInFunction
            << "Blending for " << name << "s does not apply "
            << "any configuration-specific modelling, but no general model "
            << "was provided for " << phaseInterface(phase1, phase2).name()
            << ". Consider adding a general model for these phases, or if no "
            << "model is needed then add a \"none\" model to suppress this "
            << "warning." << endl;
    }

    if (can1In2 && !valid.general && !valid.oneDispersedInTwo)
    {
        WarningInFunction
            << "Blending for " << name << "s permits "
            << phase2.name() << " to become continuous, but no model was "
            << "provided for " << dispersedPhaseInterface(phase1, phase2).name()
            << ". Consider adding a model for this configuration (or for "
            << phaseInterface(phase1, phase2).name() << "), or if no model is "
            << "needed then add a \"none\" model to suppress this warning."
            << endl;
    }

    if (can2In1 && !valid.general && !valid.twoDispersedInOne)
    {
        WarningInFunction
            << "Blending for " << name << "s permits "
            << phase1.name() << " to become continuous, but no model was "
            << "provided for " << dispersedPhaseInterface(phase2, phase1).name()
            << ". Consider adding a model for this configuration (or for "
            << phaseInterface(phase1, phase2).name() << "), or if no model is "
            << "needed then add a \"none\" model to suppress this warning."
            << endl;
    }

    if (canS && !valid.general && !valid.segregated)
    {
        WarningInFunction
            << "Blending for " << name << "s permits "
            << "segregation but no model was provided for "
            << segregatedPhaseInterface(phase2, phase1).name()
            << ". Consider adding a model for this configuration (or for "
            << phaseInterface(phase1, phase2).name() << "), or if no model is "
            << "needed then add a \"none\" model to suppress this warning."
            << endl;
    }
}


template<class GeoMesh>
void Foam::blendedInterfacialModel::calculateBlendingCoeffs
(
    const UPtrList<const volScalarField>& alphas,
    const boolSet& valid,
    TmpSet<GeometricField<scalar, GeoMesh>>& fs
) const
{
    typedef GeometricField<scalar, GeoMesh> scalarGeoField;

    tmp<scalarGeoField>& fG = fs.general;
    tmp<scalarGeoField>& f1D2 = fs.oneDispersedInTwo;
    tmp<scalarGeoField>& f2D1 = fs.twoDispersedInOne;
    tmp<scalarGeoField>& fS = fs.segregated;
    PtrList<scalarGeoField>& fGD = fs.generalDisplaced;
    PtrList<scalarGeoField>& f1D2D = fs.oneDispersedInTwoDisplaced;
    PtrList<scalarGeoField>& f2D1D = fs.twoDispersedInOneDisplaced;
    PtrList<scalarGeoField>& fSD = fs.segregatedDisplaced;

    const bool can1In2 = blending_->canBeContinuous(1);
    const bool can2In1 = blending_->canBeContinuous(0);
    const bool canS = blending_->canSegregate();

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
    tmp<scalarGeoField> F1D2, F2D1, FS;
    if (valid.oneDispersedInTwo || any(valid.oneDispersedInTwoDisplaced))
    {
        F1D2 =
            interpolate<scalarGeoField>
            (
                blending_->f1DispersedIn2(alphas)
            );
    }
    if (valid.twoDispersedInOne || any(valid.twoDispersedInOneDisplaced))
    {
        F2D1 =
            interpolate<scalarGeoField>
            (
                blending_->f2DispersedIn1(alphas)
            );
    }
    if (valid.segregated || any(valid.segregatedDisplaced))
    {
        FS =
            interpolate<scalarGeoField>
            (
                blending_->f12Segregated(alphas)
            );
    }

    // Construct non-displaced coefficients
    {
        if (can1In2 && valid.oneDispersedInTwo)
        {
            f1D2 = F1D2().clone();
        }

        if (can2In1 && valid.twoDispersedInOne)
        {
            f2D1 = F2D1().clone();
        }

        if (canS && valid.segregated)
        {
            fS = FS().clone();
        }

        if (valid.general)
        {
            fG = constant(1);
            if (f1D2.valid()) { fG.ref() -= f1D2(); }
            if (f2D1.valid()) { fG.ref() -= f2D1(); }
            if (fS.valid()) { fG.ref() -= fS(); }
        }
    }

    // Construct displaced coefficients
    tmp<scalarGeoField> fDSum;
    if
    (
        any(valid.generalDisplaced)
     || any(valid.oneDispersedInTwoDisplaced)
     || any(valid.twoDispersedInOneDisplaced)
     || any(valid.segregatedDisplaced)
    )
    {
        fDSum = constant(0);

        forAll(alphas, phasei)
        {
            const phaseModel& phaseD = interface_.fluid().phases()[phasei];

            if (interface_.contains(phaseD)) continue;

            // Get the displaced blending functions
            tmp<scalarGeoField> FD =
                interpolate<scalarGeoField>
                (
                    blending_->f12DisplacedBy3(alphas, phasei)
                );

            if (can1In2 && valid.oneDispersedInTwoDisplaced[phasei])
            {
                f1D2D.set(phasei, FD()*F1D2());
                fDSum.ref() += f1D2D[phasei];
            }

            if (can2In1 && valid.twoDispersedInOneDisplaced[phasei])
            {
                f2D1D.set(phasei, FD()*F2D1());
                fDSum.ref() += f2D1D[phasei];
            }

            if (canS && valid.segregatedDisplaced[phasei])
            {
                fSD.set(phasei, FD()*FS());
                fDSum.ref() += fSD[phasei];
            }

            if (valid.generalDisplaced[phasei])
            {
                fGD.set(phasei, FD());
                if (f1D2D.set(phasei)) fGD[phasei] -= f1D2D[phasei];
                if (f2D1D.set(phasei)) fGD[phasei] -= f2D1D[phasei];
                if (fSD.set(phasei)) fGD[phasei] -= fSD[phasei];
                fDSum.ref() += fGD[phasei];
            }
        }
    }

    // Remove the displaced part from the non-displaced coefficients. Maintain
    // the hierarchy by reducing the coefficient of the general model first.
    if (fDSum.valid())
    {
        tmp<scalarGeoField> fRemove(fDSum);

        auto remove = [&fRemove](tmp<scalarGeoField>& f)
        {
            tmp<scalarGeoField> df = min(f(), fRemove());
            f.ref() -= df();
            fRemove.ref() -= df();
        };

        tmp<scalarGeoField> fSumNotG = constant(0);
        if (f1D2.valid()) fSumNotG.ref() += f1D2();
        if (f2D1.valid()) fSumNotG.ref() += f2D1();
        if (fS.valid()) fSumNotG.ref() += fS();

        if (fG.valid())
        {
            remove(fG);
        }
        else
        {
            fRemove.ref() =
                max(fRemove() - (scalar(1) - fSumNotG()), scalar(0));
        }

        tmp<scalarGeoField> fSumNotG0 = fSumNotG().clone();

        remove(fSumNotG);

        tmp<scalarGeoField> factor = fSumNotG/max(fSumNotG0, vSmall);
        if (f1D2.valid()) f1D2.ref() *= factor();
        if (f2D1.valid()) f2D1.ref() *= factor();
        if (fS.valid()) fS.ref() *= factor();
    }
}


template
void Foam::blendedInterfacialModel::calculateBlendingCoeffs<Foam::fvMesh>
(
    const UPtrList<const volScalarField>& alphas,
    const boolSet& valid,
    TmpSet<volScalarField>& fs
) const;


template
void Foam::blendedInterfacialModel::calculateBlendingCoeffs<Foam::surfaceMesh>
(
    const UPtrList<const volScalarField>& alphas,
    const boolSet& valid,
    TmpSet<surfaceScalarField>& fs
) const;


void Foam::blendedInterfacialModel::writeBlendingCoefficients
(
    const word& name,
    const boolSet& valid
) const
{
    check(name, valid);

    const phaseSystem& fluid = interface_.fluid();
    const label nPhases = fluid.phases().size();

    // Get the blending coefficients
    TmpSet<volScalarField> fs(nPhases);
    tmp<volScalarField>& fG = fs.general;
    tmp<volScalarField>& f1D2 = fs.oneDispersedInTwo;
    tmp<volScalarField>& f2D1 = fs.twoDispersedInOne;
    tmp<volScalarField>& fS = fs.segregated;
    PtrList<volScalarField>& fGD = fs.generalDisplaced;
    PtrList<volScalarField>& f1D2D = fs.oneDispersedInTwoDisplaced;
    PtrList<volScalarField>& f2D1D = fs.twoDispersedInOneDisplaced;
    PtrList<volScalarField>& fSD = fs.segregatedDisplaced;
    calculateBlendingCoeffs
    (
        fluid.phases().PtrList<phaseModel>::convert<const volScalarField>(),
        valid,
        fs
    );

    const word prefix(name, ':', interface_.name(), ':');

    Info<< indent << "Writing blending coefficients" << endl;

    auto write = [&name]
    (
        const phaseInterface& interface,
        const word& fName,
        volScalarField& f
    )
    {
        f.rename(word(name, ':', interface.name(), ':', fName));
        f.write();
    };

    if (fG.valid()) write(interface_, "fG", fG.ref());
    if (f1D2.valid()) write(interface_, "f1D2", f1D2.ref());
    if (f2D1.valid()) write(interface_, "f2D1", f2D1.ref());
    if (fS.valid()) write(interface_, "fS", fS.ref());

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phaseD = fluid.phases()[phasei];

        if (interface_.contains(phaseD)) continue;

        const displacedPhaseInterface interfaceD
        (
            interface_.phase1(),
            interface_.phase2(),
            phaseD
        );

        if (fGD.set(phasei)) write(interfaceD, "fG", fGD[phasei]);
        if (f1D2D.set(phasei)) write(interfaceD, "f1D2", f1D2D[phasei]);
        if (f2D1D.set(phasei)) write(interfaceD, "f2D1", f2D1D[phasei]);
        if (fSD.set(phasei)) write(interfaceD, "fS", fSD[phasei]);
    }
}


void Foam::blendedInterfacialModel::postProcessBlendingCoefficients
(
    const word& name,
    const boolSet& valid,
    const entry& blendingEntry
) const
{
    if (blendingEntry.isDict() && blendingEntry.dict().found("format"))
    {
        const word format = blendingEntry.dict().lookup<word>("format");

        postProcessBlendingCoefficients(name, valid, format);
    }
}


void Foam::blendedInterfacialModel::postProcessBlendingCoefficients
(
    const word& name,
    const boolSet& valid,
    const word& format
) const
{
    check(name, valid);

    const phaseSystem& fluid = interface_.fluid();
    const label nPhases = fluid.phases().size();

    // Don't bother if we aren't going to or cannot write anything
    if (nPhases <= 2 && format == noSetWriter::typeName) return;
    if (nPhases > 2 && format == noSurfaceWriter::typeName) return;
    if (!blending_->functionOfAlphas())
    {
        WarningInFunction
            << "Cannot write the blending coefficients for blending "
            << "method " << blending_->type() << " as this method is not "
            << "a pure function of the volume fractions" << endl;
        return;
    }

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
        fieldNames[0] = fluid.phases()[0].volScalarField::name();
        fieldNames[1] = fluid.phases()[1].volScalarField::name();
        fields.set(0, new scalarField(linearSequence01(nDivisions)));
        fields.set(1, new scalarField(1 - fields[0]));
    }
    else
    {
        // Polygon with as many vertices as there are phases. Each phase
        // fraction equals one at a unique vertex, they vary linearly between
        // vertices, and vary smoothly in the interior of the polygon. The sum
        // of phase fractions is always one throughout the polygon (partition of
        // unity). This is done using Wachspress coordinates. This does not
        // cover the entire N-dimensional phase fraction space for N >= 4 (it is
        // hard to imagine how that could be visualised) but it provides enough
        // to give a good indication of what is going on.

        // Create the nodes of the blending space polygon
        List<point> phaseNodes(nPhases);
        forAll(phaseNodes, phasei)
        {
            const scalar theta = 2*constant::mathematical::pi*phasei/nPhases;
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
            TmpSet<volScalarField> fs(nPhases);
            tmp<volScalarField>& fG = fs.general;
            tmp<volScalarField>& f1D2 = fs.oneDispersedInTwo;
            tmp<volScalarField>& f2D1 = fs.twoDispersedInOne;
            tmp<volScalarField>& fS = fs.segregated;
            PtrList<volScalarField>& fGD = fs.generalDisplaced;
            PtrList<volScalarField>& f1D2D = fs.oneDispersedInTwoDisplaced;
            PtrList<volScalarField>& f2D1D = fs.twoDispersedInOneDisplaced;
            PtrList<volScalarField>& fSD = fs.segregatedDisplaced;
            calculateBlendingCoeffs
            (
                alphas.convert<const volScalarField>(),
                valid,
                fs
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

            TmpSet<volScalarField> fs(nPhases);
            tmp<volScalarField>& fG = fs.general;
            tmp<volScalarField>& f1D2 = fs.oneDispersedInTwo;
            tmp<volScalarField>& f2D1 = fs.twoDispersedInOne;
            tmp<volScalarField>& fS = fs.segregated;
            PtrList<volScalarField>& fGD = fs.generalDisplaced;
            PtrList<volScalarField>& f1D2D = fs.oneDispersedInTwoDisplaced;
            PtrList<volScalarField>& f2D1D = fs.twoDispersedInOneDisplaced;
            PtrList<volScalarField>& fSD = fs.segregatedDisplaced;
            calculateBlendingCoeffs
            (
                alphas.convert<const volScalarField>(),
                valid,
                fs
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
       /name;
    Info<< indent <<  "Writing blending coefficients to "
        << path/interface_.name() << endl;
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendedInterfacialModel::blendedInterfacialModel
(
    const word& name,
    const word& typeName,
    const phaseInterface& interface,
    const entry& blendingEntry
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(name, interface.name()),
            interface.fluid().mesh().time().name(),
            interface.fluid().mesh()
        )
    ),
    interface_(interface),
    blending_(blendingMethod::New(typeName, blendingEntry, interface)),
    checkTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendedInterfacialModel::~blendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::phaseInterface&
Foam::blendedInterfacialModel::interface() const
{
    return interface_;
}


bool Foam::blendedInterfacialModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
