/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "FickianEddyDiffusivity.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "Function2Evaluate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
FickianEddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo,
        false
    ),

    Sct_
    (
        dimensioned<scalar>
        (
            "Sct",
            dimless,
            this->coeffDict_
        )
    ),

    D_(this->thermo().composition().species().size()),
    DT_
    (
        this->coeffDict_.found("DT")
      ? this->thermo().composition().species().size()
      : 0
    )
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if
    (
        unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
        ::read()
    )
    {
        Sct_.readIfPresent(this->coeffDict());

        const speciesTable& species = this->thermo().composition().species();
        const dictionary& Ddict = this->coeffDict_.subDict("D");

        forAll(species, i)
        {
            D_.set(i, Function2<scalar>::New(species[i], Ddict).ptr());
        }

        if (this->coeffDict_.found("DT"))
        {
            const dictionary& DTdict = this->coeffDict_.subDict("DT");

            forAll(species, i)
            {
                DT_.set(i, Function2<scalar>::New(species[i], DTdict).ptr());
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<volScalarField>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition =
        this->thermo().composition();

    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()
       *evaluate
        (
            D_[composition.index(Yi)],
            dimViscosity,
            this->thermo().p(),
            this->thermo().T()
        )
      + (this->Prt_/Sct_)*this->alphat()
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<scalarField>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const basicSpecieMixture& composition =
        this->thermo().composition();

    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *D_[composition.index(Yi)].value
        (
            this->thermo().p().boundaryField()[patchi],
            this->thermo().T().boundaryField()[patchi]
        )
      + this->Prt_.value()/Sct_.value()*this->alphat(patchi);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    if (DT_.size())
    {
        const basicSpecieMixture& composition = this->thermo().composition();
        const volScalarField& p = this->thermo().T();
        const volScalarField& T = this->thermo().T();

        return
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
            j(Yi)
          - fvc::interpolate
            (
                evaluate(DT_[composition.index(Yi)], dimDynamicViscosity, p, T)
            )
           *fvc::snGrad(T)/fvc::interpolate(T);
    }
    else
    {
        return
            unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
            j(Yi);
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->momentumTransport().alphaRhoPhi().group()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();

    if (Y.size())
    {
        surfaceScalarField sumJ
        (
            surfaceScalarField::New
            (
                "sumJ",
                Y[0].mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(sumJ.dimensions()*dimEnergy/dimMass, 0)
            )
        );

        forAll(Y, i)
        {
            if (i != composition.defaultSpecie())
            {
                const volScalarField hi
                (
                    composition.HE(i, this->thermo().p(), this->thermo().T())
                );

                const surfaceScalarField ji(this->j(Y[i]));
                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = composition.defaultSpecie();

            const volScalarField hi
            (
                composition.HE(i, this->thermo().p(), this->thermo().T())
            );

            sumJh -= sumJ*fvc::interpolate(hi);
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
FickianEddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    tmp<fvScalarMatrix> tmpDivq
    (
        fvm::Su
        (
            -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
            he
        )
    );

    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();

    if (!Y.size())
    {
        tmpDivq.ref() -=
            correction(fvm::laplacian(this->alpha()*this->alphaEff(), he));
    }
    else
    {
        tmpDivq.ref() -= fvm::laplacian(this->alpha()*this->alphaEff(), he);

        volScalarField heNew
        (
            volScalarField::New
            (
                "he",
                he.mesh(),
                dimensionedScalar(he.dimensions(), 0)
            )
        );

        surfaceScalarField sumJ
        (
            surfaceScalarField::New
            (
                "sumJ",
                he.mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                he.mesh(),
                dimensionedScalar(sumJ.dimensions()*he.dimensions(), 0)
            )
        );

        forAll(Y, i)
        {
            if (i != composition.defaultSpecie())
            {
                const volScalarField hi
                (
                    composition.HE(i, this->thermo().p(), this->thermo().T())
                );

                heNew += Y[i]*hi;

                const surfaceScalarField ji(this->j(Y[i]));
                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = composition.defaultSpecie();

            const volScalarField hi
            (
                composition.HE(i, this->thermo().p(), this->thermo().T())
            );

            heNew += Y[i]*hi;

            sumJh -= sumJ*fvc::interpolate(hi);
        }

        tmpDivq.ref() +=
            fvc::laplacian(this->alpha()*this->alphaEff(), heNew);

        tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());
    }

    return tmpDivq;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
