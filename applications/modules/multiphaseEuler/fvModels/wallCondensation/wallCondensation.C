/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "wallCondensation.H"

#include "multiphaseEuler.H"

#include "fluidMulticomponentThermo.H"
#include "fluidThermophysicalTransportModel.H"

#include "saturationPressureModel.H"

#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"
#include "wallCondensationPhaseChangeRateFvPatchScalarField.H"
#include "zeroGradientFvPatchFields.H"

#include "addToRunTimeSelectionTable.H"

/*---------------------------------------------------------------------------*\
                Class wallCondensation::properties Declaration
\*---------------------------------------------------------------------------*/

struct Foam::fv::wallCondensation::properties
{
    //- Typedef to shorten the name of the Jayatilleke wall function
    typedef
        compressible::alphatJayatillekeWallFunctionFvPatchScalarField
        alphatJayatillekeWallFunction;

    //- Wall function field
    const wallCondensation& model;

    //- Patch index
    const label patchi;

    //- Liquid volume fraction
    const scalarField& alphaLiquid;

    //- Vapour volume fraction
    const scalarField& alphaVapour;

    //- Liquid thermophysical transport model
    const fluidThermophysicalTransportModel& ttmLiquid;

    //- Vapour thermophysical transport model
    const fluidThermophysicalTransportModel& ttmVapour;

    //- Liquid convective turbulent thermal diffusivity
    const scalarField alphatConvLiquid;

    //- Phase convective turbulent thermal diffusivity
    const scalarField alphatConvVapour;

    //- Patch area by neighbouring cell volume ratio
    const scalarField AbyV;

    //- Vapour heat capacity
    const scalarField& CpVapour;

    //- Cell temperature
    const scalarField TcVapour;

    //- Patch temperature
    const scalarField TwVapour;

    //- Saturation pressure
    const scalarField pSat;

    //- Latent heat
    const scalarField L;

    //- Patch
    inline const fvPatch& patch() const
    {
        return model.mesh().boundary()[patchi];
    }

    //- Constructor
    properties
    (
        const wallCondensation& model,
        const label patchi
    )
    :
        model(model),
        patchi(patchi),
        alphaLiquid(model.liquid_.boundaryField()[patchi]),
        alphaVapour(model.vapour_.boundaryField()[patchi]),
        ttmLiquid
        (
            model.mesh().lookupType<fluidThermophysicalTransportModel>
            (
                model.liquid_.name()
            )
        ),
        ttmVapour
        (
            model.mesh().lookupType<fluidThermophysicalTransportModel>
            (
                model.vapour_.name()
            )
        ),
        alphatConvLiquid
        (
            alphatJayatillekeWallFunction::alphat(ttmLiquid, model.Prt_, patchi)
        ),
        alphatConvVapour
        (
            alphatJayatillekeWallFunction::alphat(ttmVapour, model.Prt_, patchi)
        ),
        AbyV
        (
            patch().magSf()
           /scalarField
            (
                patch().boundaryMesh().mesh().V(),
                patch().faceCells()
            )
        ),
        CpVapour(model.vapour_.thermo().Cp().boundaryField()[patchi]),
        TcVapour
        (
            model.vapour_.thermo().T().boundaryField()[patchi]
           .patchInternalField()
        ),
        TwVapour
        (
            model.vapour_.thermo().T().boundaryField()[patchi]
        ),
        pSat
        (
            model.saturationModelPtr_->pSat
            (
                TwVapour
            )
        ),
        L
        (
            model.L
            (
                patchi,
                TwVapour
            )
        )
    {}
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(wallCondensation, 0);
    addToRunTimeSelectionTable(fvModel, wallCondensation, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::wallCondensation::readCoeffs(const dictionary& dict)
{
    reReadSpecie(dict);

    saturationModelPtr_.reset
    (
        saturationPressureModel::New
        (
            "saturationPressure",
            dict
        ).ptr()
    );

    Prt_ = dict.lookupOrDefault<scalar>("Prt", dimless, 0.85);

    Sct_ = dict.lookupOrDefault<scalar>("Sct", dimless, 0.7);

    specieSemiImplicit_ =
        dict.lookupOrDefault<bool>("specieSemiImplicit", false);
}


void Foam::fv::wallCondensation::correctMDot() const
{
    Info<< type() << ": " << name() << endl << incrIndent;

    //- Reset the phase-change rates in all the near-wall cells
    forAll(mDot_.boundaryField(), patchi)
    {
        if (!isPatchActive(patchi)) continue;

        const labelUList& faceCells = mesh().boundary()[patchi].faceCells();
        forAll(faceCells, i)
        {
            mDot_[faceCells[i]] = scalar(0);
            mDotDy_[faceCells[i]] = scalar(0);
        }
    }

    // Loop the active patches, evaluate the model, and sum the phase change
    // rates into the adjacent cells
    forAll(mDot_.boundaryField(), patchi)
    {
        if (!isPatchActive(patchi)) continue;

        // Access the wall-condensation phase-change patch field for this patch
        wallCondensationPhaseChangeRateFvPatchScalarField& mDot =
            mDotPfRef(patchi);

        scalarField& mDotDy = mDotDy_.boundaryFieldRef()[patchi];

        // Construct properties
        const properties props(*this, patchi);

        const fluidMulticomponentThermo& thermo =
            fluidMulticomponentThermos(true, false)[0];

        const label speciei = specieis()[0];

        // Calculate effective diffusivity
        const scalarField DEff
        (
            props.ttmVapour.D(vapour_.Y(species()[0]), patchi)
          + Prt_/Sct_*props.alphatConvVapour
        );

        const scalarField xc
        (
            vapour_
           .Y(species()[0]).boundaryField()[patchi].patchInternalField()
           /thermo.WiValue(speciei)
           *thermo.W(patchi) // Assuming zeroGradient for species
        );

        const scalarField xw(props.pSat/thermo.p().boundaryField()[patchi]);

        mDot =
          - props.alphaVapour
           *props.AbyV
           *thermo.WiValue(speciei)/thermo.W(patchi)*DEff
           *props.patch().deltaCoeffs()
           *log(max(1 - xc, 0.001)/max(1 - xw, 0.001));

        if (specieSemiImplicit_)
        {
            mDotDy =
                props.alphaVapour
               *props.AbyV
               *DEff
               *props.patch().deltaCoeffs()
               /max(1 - xc, 0.001)*pos(xc - xw);
        }
        else
        {
            mDotDy = Zero;
        }

        // Only allow condensation
        mDot.condensing_ = pos(mDot);
        mDot = max(mDot, scalar(0));

        const scalarField gradT
        (
            props.patch().deltaCoeffs()
           *min(props.TwVapour - props.TcVapour, -rootSmall*props.TcVapour)
        );

        const scalarField q(mDot*props.L/props.AbyV);

        const scalarField alphatCondensingVapour
        (
            q/props.CpVapour/gradT/max(props.alphaVapour, rootSmall)
        );

        mDot.alphatVapour_ = props.alphatConvVapour + alphatCondensingVapour;
        mDot.alphatLiquid_ = props.alphatConvLiquid;

        infoField
        (
            "mDot[" + mesh().boundary()[patchi].name() + "]",
            dimDensity/dimTime,
            mDot
        );

        // Sum the phase change rate into the internal field
        const labelUList& faceCells = mesh().boundary()[patchi].faceCells();
        forAll(faceCells, i)
        {
            mDot_[faceCells[i]] += mDot[i];
            mDotDy_[faceCells[i]] += mDotDy[i];
        }
    }

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::wallCondensation::wallCondensation
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    wallPhaseChange
    (
        name,
        modelType,
        mesh,
        dict,
        readSpecie(coeffs(modelType, dict), true)
    ),
    liquid_(phases().second()),
    vapour_(phases().first()),
    alphatLiquid_(wallPhaseChange::alphats().second()),
    alphatVapour_(wallPhaseChange::alphats().first()),
    p_rgh_
    (
        mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName).p_rgh
    ),
    Prt_(NaN),
    Sct_(NaN),
    saturationModelPtr_(nullptr),
    pressureEquationIndex_(-1),
    specieSemiImplicit_(false),
    mDot_
    (
        IOobject
        (
            name + ":mDot",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, scalar(0)),
        mDotBoundaryTypes
        (
            wallCondensationPhaseChangeRateFvPatchScalarField::typeName
        )
    ),
    mDotDy_
    (
        IOobject
        (
            name + ":mDotDy",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, scalar(0))
    )
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField& Foam::fv::wallCondensation::active
(
    const label patchi
) const
{
    return mDotPfRef(patchi).condensing_;
}


Foam::Pair<const Foam::scalarField&> Foam::fv::wallCondensation::alphats
(
    const label patchi
) const
{
    return
        Pair<const Foam::scalarField&>
        (
            mDotPfRef(patchi).alphatVapour_,
            mDotPfRef(patchi).alphatLiquid_
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallCondensation::Lfraction() const
{
    // Put all the latent heat into the vapour
    return
        volScalarField::Internal::New
        (
            name() + ":Lfraction",
            mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallCondensation::mDot() const
{
    return mDot_.internalField();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallCondensation::mDotDy() const
{
    return mDotDy_.internalField();
}


const Foam::wallCondensationPhaseChangeRateFvPatchScalarField&
Foam::fv::wallCondensation::mDotPf(const label patchi) const
{
    if (!isPatchActive(patchi))
    {
        FatalErrorInFunction
            << "Patch " << mesh().boundary()[patchi].name()
            << " is not a condensation phase change wall" << exit(FatalError);
    }

    return
        refCast<const wallCondensationPhaseChangeRateFvPatchScalarField>
        (
            mDot_.boundaryField()[patchi]
        );
}


Foam::wallCondensationPhaseChangeRateFvPatchScalarField&
Foam::fv::wallCondensation::mDotPfRef(const label patchi) const
{
    if (!isPatchActive(patchi))
    {
        FatalErrorInFunction
            << "Patch " << mesh().boundary()[patchi].name()
            << " is not a condensation phase change wall" << exit(FatalError);
    }

    return
        refCast<wallCondensationPhaseChangeRateFvPatchScalarField>
        (
            mDot_.boundaryFieldRef()[patchi]
        );
}


void Foam::fv::wallCondensation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &liquid_ || &alpha == &vapour_)
     && (&rho == &liquid_.rho() || &rho == &vapour_.rho())
     && &eqn.psi() == &p_rgh_
    )
    {
        // Ensure that the source is up-to date if this is the first call in
        // the current phase loop
        if (pressureEquationIndex_ % 2 == 0) correctMDot();
        pressureEquationIndex_ ++;
    }

    // Let the base class add the actual source
    massTransfer::addSup(alpha, rho, eqn);
}


void Foam::fv::wallCondensation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& heOrYi,
    fvMatrix<scalar>& eqn
) const
{
    const label i = this->index(phaseNames(), alpha.group());

    if (!specieSemiImplicit_ || i != 1)
    {
        return phaseChange::addSup(alpha, rho, heOrYi, eqn);
    }

    const label s = this->sign(phaseNames(), alpha.group());

    const fluidMulticomponentThermo& thermo =
        fluidMulticomponentThermos(true, false)[0];

    const word specieName = heOrYi.member();

    // Mass fraction equation
    if (thermo.containsSpecie(specieName))
    {
        // A non-transferring specie. Do not add a source.
        if (!species().found(specieName)) return;

        // The transferring specie. Add a linearised source.
        tmp<volScalarField::Internal> tmDot = this->mDot();
        tmp<volScalarField::Internal> tmDotDy = this->mDotDy();

        eqn += s*(tmDot() + correction(fvm::Sp(tmDotDy, eqn.psi())));

        return;
    }

    // Something else. Fall back.
    phaseChange::addSup(alpha, rho, heOrYi, eqn);
}


void Foam::fv::wallCondensation::correct()
{
    // Reset the p_rgh equation solution counter
    pressureEquationIndex_ = 0;

    // Correct the total phase change rate
    correctMDot();
}


bool Foam::fv::wallCondensation::read(const dictionary& dict)
{
    if (phaseChange::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
