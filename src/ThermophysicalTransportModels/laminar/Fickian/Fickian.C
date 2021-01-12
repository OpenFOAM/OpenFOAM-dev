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

#include "Fickian.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "Function2Evaluate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Fickian<BasicThermophysicalTransportModel>::Fickian
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    BasicThermophysicalTransportModel
    (
        type,
        momentumTransport,
        thermo
    ),

    D_(this->thermo().composition().species().size()),

    DT_
    (
        this->coeffDict_.found("DT")
      ? this->thermo().composition().species().size()
      : 0
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool Fickian<BasicThermophysicalTransportModel>::read()
{
    if
    (
        BasicThermophysicalTransportModel::read()
    )
    {
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


template<class BasicThermophysicalTransportModel>
tmp<volScalarField> Fickian<BasicThermophysicalTransportModel>::DEff
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
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> Fickian<BasicThermophysicalTransportModel>::DEff
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
        );
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> Fickian<BasicThermophysicalTransportModel>::q() const
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


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> Fickian<BasicThermophysicalTransportModel>::divq
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


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> Fickian<BasicThermophysicalTransportModel>::j
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
            BasicThermophysicalTransportModel::j(Yi)
          - fvc::interpolate
            (
                evaluate(DT_[composition.index(Yi)], dimDynamicViscosity, p, T)
            )
           *fvc::snGrad(T)/fvc::interpolate(T);
    }
    else
    {
        return BasicThermophysicalTransportModel::j(Yi);
    }
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> Fickian<BasicThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    if (DT_.size())
    {
        const basicSpecieMixture& composition = this->thermo().composition();
        const volScalarField& p = this->thermo().T();
        const volScalarField& T = this->thermo().T();

        return
            BasicThermophysicalTransportModel::divj(Yi)
          - fvc::div
            (
                fvc::interpolate
                (
                    evaluate
                    (
                        DT_[composition.index(Yi)],
                        dimDynamicViscosity,
                        p,
                        T
                    )
                )
               *fvc::snGrad(T)/fvc::interpolate(T)
               *T.mesh().magSf()
            );
    }
    else
    {
        return BasicThermophysicalTransportModel::divj(Yi);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
