/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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

#include "comfort.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(comfort, 0);
    addToRunTimeSelectionTable(functionObject, comfort, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::functionObjects::comfort::magU() const
{
    tmp<volScalarField> tmagU = mag(lookupObject<volVectorField>("U"));

    // Switch to use the averaged velocity field in the domain.
    // Consistent with EN ISO 7730 but does not make physical sense
    if (meanVelocity_)
    {
        tmagU.ref() = tmagU->weightedAverage(mesh_.V());
    }

    return tmagU;
}


Foam::dimensionedScalar Foam::functionObjects::comfort::Trad() const
{
    dimensionedScalar Trad(Trad_);

    // The mean radiation is calculated by the mean wall temperatures
    // which are summed and divided by the area | only walls are taken into
    // account. This approach might be correct for a squared room but will
    // definitely be inconsistent for complex room geometries. The norm does
    // not provide any information about the calculation of this quantity.
    if (!TradSet_)
    {
        const volScalarField::Boundary& TBf =
            lookupObject<volScalarField>("T").boundaryField();

        scalar areaIntegral = 0;
        scalar TareaIntegral = 0;

        forAll(TBf, patchi)
        {
            const fvPatchScalarField& pT = TBf[patchi];
            const fvPatch& pTBf = TBf[patchi].patch();
            const scalarField& pSf = pTBf.magSf();

            if (isType<wallFvPatch>(pTBf))
            {
                areaIntegral += gSum(pSf);
                TareaIntegral += gSum(pSf*pT);
            }
        }

        Trad.value() = TareaIntegral/areaIntegral;
    }

    // Bounds based on EN ISO 7730
    if ((Trad.value() < 283.15) || (Trad.value() > 313.15))
    {
        WarningInFunction
            << "The calculated mean wall radiation temperature is out of the\n"
            << "bounds specified in EN ISO 7730:2006\n"
            << "Valid range is 10 degC < T < 40 degC\n"
            << "The actual value is: " << Trad.value() - 273.15 << nl << endl;
    }

    return Trad;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::comfort::pSat() const
{
    static const dimensionedScalar kPaToPa(dimPressure, 1000);
    static const dimensionedScalar A(dimless, 16.6563);
    static const dimensionedScalar B(dimTemperature, 4030.183);
    static const dimensionedScalar C(dimTemperature, -38.15);

    tmp<volScalarField> tpSat(volScalarField::New("pSat", mesh_, pSat_));

    // Calculate the saturation pressure if no user input is given
    if (pSat_.value() == 0)
    {
        const volScalarField& T = lookupObject<volScalarField>("T");

        // Equation based on ISO 7730:2006
        tpSat = kPaToPa*exp(A - B/(T + C));
    }

    return tpSat;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::comfort::Tcloth
(
    volScalarField& hc,
    const dimensionedScalar& metabolicRateSI,
    const dimensionedScalar& extWorkSI,
    const volScalarField& T,
    const dimensionedScalar& Trad
)
{
    const dimensionedScalar factor1(dimTemperature, 308.85);

    const dimensionedScalar factor2
    (
        dimTemperature/metabolicRateSI.dimensions(),
        0.028
    );

    const dimensionedScalar factor3
    (
        dimensionSet(1, 0, -3, -4, 0, 0, 0),
        3.96e-8
    );

    // Heat transfer coefficient based on forced convection [W/m^2/K]
    const volScalarField hcForced
    (
        dimensionedScalar(hc.dimensions()/sqrt(dimVelocity), 12.1)
       *sqrt(magU())
    );

    // Tcl [K] (surface cloth temperature)
    tmp<volScalarField> tTcl
    (
        volScalarField::New
        (
            "Tcl",
            T.mesh(),
            dimTemperature
        )
    );

    volScalarField& Tcl = tTcl.ref();

    // Initial guess
    Tcl = T;

    label i = 0;

    Tcl.storePrevIter();

    // Iterative solving of equation (2)
    do
    {
        Tcl = (Tcl + Tcl.prevIter())/2;
        Tcl.storePrevIter();

        // Heat transfer coefficient based on natural convection
        volScalarField hcNatural
        (
            dimensionedScalar(hc.dimensions()/pow025(dimTemperature), 2.38)
           *pow025(mag(Tcl - T))
        );

        // Set heat transfer coefficient based on equation (3)
        hc =
            pos(hcForced - hcNatural)*hcForced
          + neg0(hcForced - hcNatural)*hcNatural;

        // Calculate surface temperature based on equation (2)
        Tcl =
            factor1
          - factor2*(metabolicRateSI - extWorkSI)
          - Icl_*factor3*fcl_*(pow4(Tcl) - pow4(Trad))
          - Icl_*fcl_*hc*(Tcl - T);

    } while (!converged(Tcl) && i++ < maxClothIter_);

    if (i == maxClothIter_)
    {
        WarningInFunction
            << "The surface cloth temperature did not converge within " << i
            << " iterations\n";
    }

    return tTcl;
}


bool Foam::functionObjects::comfort::converged
(
    const volScalarField& phi
) const
{
    return
        max(mag(phi.primitiveField() - phi.prevIter().primitiveField()))
      < tolerance_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::comfort::comfort
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    clothing_("clothing", dimless, 0),
    metabolicRate_("metabolicRate", dimMass/pow3(dimTime), 0.8),
    extWork_("extWork", dimMass/pow3(dimTime), 0),
    TradSet_(false),
    Trad_("Trad", dimTemperature, 0),
    relHumidity_("relHumidity", dimless, 0.5),
    pSat_("pSat", dimPressure, 0),
    Icl_("Icl", dimensionSet(-1, 0, 3, 1, 0, 0, 0), 0),
    fcl_("fcl", dimless, 0),
    tolerance_(1e-4),
    maxClothIter_(100),
    meanVelocity_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::comfort::~comfort()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::comfort::read(const dictionary& dict)
{
    clothing_.readIfPresent(dict);
    metabolicRate_.readIfPresent(dict);
    extWork_.readIfPresent(dict);
    pSat_.readIfPresent(dict);
    tolerance_ = dict.lookupOrDefault("tolerance", 1e-4);
    maxClothIter_ = dict.lookupOrDefault("maxClothIter", 100);
    meanVelocity_ = dict.lookupOrDefault<Switch>("meanVelocity", false);

    // Read relative humidity if provided and convert from % to fraction
    if (dict.found(relHumidity_.name()))
    {
        relHumidity_.read(dict);
        relHumidity_ /= 100;
    }

    // Read radiation temperature if provided
    if (dict.found(Trad_.name()))
    {
        TradSet_ = true;
        Trad_.read(dict);
    }
    else
    {
        TradSet_ = false;
    }

    Icl_ = dimensionedScalar(Icl_.dimensions(), 0.155)*clothing_;

    fcl_.value() =
        Icl_.value() <= 0.078
      ? 1.0 + 1.290*Icl_.value()
      : 1.05 + 0.645*Icl_.value();

    return true;
}


Foam::wordList Foam::functionObjects::comfort::fields() const
{
    return wordList{"U", "T"};
}


bool Foam::functionObjects::comfort::execute()
{
    const dimensionedScalar Trad(this->Trad());
    const volScalarField pSat(this->pSat());

    const dimensionedScalar metabolicRateSI(58.15*metabolicRate_);
    const dimensionedScalar extWorkSI(58.15*extWork_);

    const volScalarField& T = lookupObject<volScalarField>("T");

    // Heat transfer coefficient [W/m^2/K]
    // This field is updated in Tcloth()
    volScalarField hc
    (
        IOobject
        (
            "hc",
            time_.name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1, 0, -3, -1, 0, 0, 0), 0)
    );

    // Calculate the surface temperature of the cloth by an iterative
    // process using equation (2) from DIN EN ISO 7730 [degC]
    const volScalarField Tcloth
    (
        this->Tcloth
        (
            hc,
            metabolicRateSI,
            extWorkSI,
            T,
            Trad
        )
    );

    // Calculate the PMV quantity
    const dimensionedScalar factor1(dimensionSet(-1, 0, 3, 0, 0, 0, 0), 0.303);
    const dimensionedScalar factor2
    (
        dimless/metabolicRateSI.dimensions(),
        -0.036
    );
    const dimensionedScalar factor3(factor1.dimensions(), 0.028);
    const dimensionedScalar factor4(dimLength/dimTime, 3.05e-3);
    const dimensionedScalar factor5(dimPressure, 5733);
    const dimensionedScalar factor6(dimTime/dimLength, 6.99);
    const dimensionedScalar factor8(metabolicRateSI.dimensions(), 58.15);
    const dimensionedScalar factor9(dimless/dimPressure, 1.7e-5);
    const dimensionedScalar factor10(dimPressure, 5867);
    const dimensionedScalar factor11(dimless/dimTemperature, 0.0014);
    const dimensionedScalar factor12(dimTemperature, 307.15);
    const dimensionedScalar factor13
    (
        dimensionSet(1, 0, -3, -4, 0, 0, 0),
        3.96e-8
    );

    const scalar factor7
    (
        // Special treatment of Term4
        // if metaRate - extWork < factor8, set to zero
        (metabolicRateSI - extWorkSI).value() < factor8.value() ? 0 : 0.42
    );

    Info<< "Calculating the predicted mean vote (PMV)\n";

    // Equation (1)
    tmp<volScalarField> PMV
    (
        volScalarField::New
        (
            "PMV",

            // Term1: Thermal sensation transfer coefficient
            (factor1*exp(factor2*metabolicRateSI) + factor3)
           *(
                (metabolicRateSI - extWorkSI)

                // Term2: Heat loss difference through skin
              - factor4
               *(
                    factor5
                  - factor6*(metabolicRateSI - extWorkSI)
                  - pSat*relHumidity_
                )

                // Term3: Heat loss through sweating
              - factor7*(metabolicRateSI - extWorkSI - factor8)

                // Term4: Heat loss through latent respiration
              - factor9*metabolicRateSI*(factor10 - pSat*relHumidity_)

                // Term5: Heat loss through dry respiration
              - factor11*metabolicRateSI*(factor12 - T)

                // Term6: Heat loss through radiation
              - factor13*fcl_*(pow4(Tcloth) - pow4(Trad))

                // Term7: Heat loss through convection
              - fcl_*hc*(Tcloth - T)
            )
        )
    );

    Info<< "Calculating the predicted percentage of dissatisfaction (PPD)\n";

    // Equation (5)
    tmp<volScalarField> PPD
    (
        volScalarField::New
        (
            "PPD",
            100 - 95*exp(-0.03353*pow4(PMV()) - 0.21790*sqr(PMV()))
        )
    );

    Info<< "Calculating the draught rating (DR)\n";

    const dimensionedScalar Umin("Umin", dimVelocity, 0.05);
    const dimensionedScalar Umax("Umax", dimVelocity, 0.5);
    const dimensionedScalar pre("preU", dimless, 0.37);
    const dimensionedScalar C1("C1", dimVelocity, 3.14);

    // Limit the velocity field to the values given in EN ISO 7733
    volScalarField Umag(mag(lookupObject<volVectorField>("U")));
    Umag.maxMin(Umin, Umax);

    // Calculate the turbulent intensity if turbulent kinetic energy field k
    // exists
    volScalarField TI
    (
        IOobject
        (
            "TI",
            time_.name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    );

    if (foundObject<volScalarField>("k"))
    {
        const volScalarField& k = lookupObject<volScalarField>("k");
        TI = sqrt(2/3*k)/Umag;
    }

    // For unit correctness
    const dimensionedScalar correctUnit
    (
        "correctUnit",
        dimensionSet(0,- 1.62, 1.62, -1, 0, 0, 0),
        1
    );

    // Equation (6)
    tmp<volScalarField> DR
    (
        volScalarField::New
        (
            "DR",
            correctUnit*(factor12 - T)*pow(Umag - Umin, 0.62)*(pre*Umag*TI + C1)
        )
    );

    return
        store(PMV)
     && store(PPD)
     && store(DR);
}


bool Foam::functionObjects::comfort::write()
{
    return
        writeObject("PMV")
     && writeObject("PPD")
     && writeObject("DR");
}


// ************************************************************************* //
