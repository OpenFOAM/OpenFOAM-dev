/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "heheuPsiThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he_.internalField();
    const scalarField& heuCells = this->heu_.internalField();
    const scalarField& pCells = this->p_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& TuCells = this->Tu_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);

        TuCells[celli] = this->cellReactants(celli).THE
        (
            heuCells[celli],
            pCells[celli],
            TuCells[celli]
        );
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pTu = this->Tu_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];

        fvPatchScalarField& ph = this->he_.boundaryField()[patchi];
        fvPatchScalarField& pheu = this->heu_.boundaryField()[patchi];

        fvPatchScalarField& pmu_ = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha_ = this->alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha_[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(ph[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha_[facei] = mixture_.alphah(pp[facei], pT[facei]);

                pTu[facei] =
                    this->patchFaceReactants(patchi, facei)
                   .THE(pheu[facei], pp[facei], pTu[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heheuPsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<psiuReactionThermo, MixtureType>(mesh, phaseName),
    Tu_
    (
        IOobject
        (
            "Tu",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    heu_
    (
        IOobject
        (
            MixtureType::thermoType::heName() + 'u',
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->heuBoundaryTypes()
    )
{
    scalarField& heuCells = this->heu_.internalField();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TuCells = this->Tu_.internalField();

    forAll(heuCells, celli)
    {
        heuCells[celli] = this->cellReactants(celli).HE
        (
            pCells[celli],
            TuCells[celli]
        );
    }

    forAll(this->heu_.boundaryField(), patchi)
    {
        fvPatchScalarField& pheu = this->heu_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = this->Tu_.boundaryField()[patchi];

        forAll(pheu, facei)
        {
            pheu[facei] = this->patchFaceReactants(patchi, facei).HE
            (
                pp[facei],
                pTu[facei]
            );
        }
    }

    this->heuBoundaryCorrection(this->heu_);

    calculate();
    this->psi_.oldTime();   // Switch on saving old time
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::~heheuPsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering heheuPsiThermo<BasicPsiThermo, MixtureType>::correct()"
            << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "exiting heheuPsiThermo<BasicPsiThermo, MixtureType>::correct()"
            << endl;
    }
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& p,
    const scalarField& Tu,
    const labelList& cells
) const
{
    tmp<scalarField> theu(new scalarField(Tu.size()));
    scalarField& heu = theu();

    forAll(Tu, celli)
    {
        heu[celli] = this->cellReactants(cells[celli]).HE(p[celli], Tu[celli]);
    }

    return theu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& p,
    const scalarField& Tu,
    const label patchi
) const
{
    tmp<scalarField> theu(new scalarField(Tu.size()));
    scalarField& heu = theu();

    forAll(Tu, facei)
    {
        heu[facei] =
            this->patchFaceReactants(patchi, facei).HE(p[facei], Tu[facei]);
    }

    return theu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::Tb() const
{
    tmp<volScalarField> tTb
    (
        new volScalarField
        (
            IOobject
            (
                "Tb",
                this->T_.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->T_
        )
    );

    volScalarField& Tb_ = tTb();
    scalarField& TbCells = Tb_.internalField();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TCells = this->T_.internalField();
    const scalarField& hCells = this->he_.internalField();

    forAll(TbCells, celli)
    {
        TbCells[celli] = this->cellProducts(celli).THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );
    }

    forAll(Tb_.boundaryField(), patchi)
    {
        fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        const fvPatchScalarField& ph = this->he_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];

        forAll(pTb, facei)
        {
            pTb[facei] =
                this->patchFaceProducts(patchi, facei)
               .THE(ph[facei], pp[facei], pT[facei]);
        }
    }

    return tTb;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::psiu() const
{
    tmp<volScalarField> tpsiu
    (
        new volScalarField
        (
            IOobject
            (
                "psiu",
                this->psi_.time().timeName(),
                this->psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi_.mesh(),
            this->psi_.dimensions()
        )
    );

    volScalarField& psiu = tpsiu();
    scalarField& psiuCells = psiu.internalField();
    const scalarField& TuCells = this->Tu_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(psiuCells, celli)
    {
        psiuCells[celli] =
            this->cellReactants(celli).psi(pCells[celli], TuCells[celli]);
    }

    forAll(psiu.boundaryField(), patchi)
    {
        fvPatchScalarField& ppsiu = psiu.boundaryField()[patchi];

        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = this->Tu_.boundaryField()[patchi];

        forAll(ppsiu, facei)
        {
            ppsiu[facei] =
                this->
                patchFaceReactants(patchi, facei).psi(pp[facei], pTu[facei]);
        }
    }

    return tpsiu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::psib() const
{
    tmp<volScalarField> tpsib
    (
        new volScalarField
        (
            IOobject
            (
                "psib",
                this->psi_.time().timeName(),
                this->psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi_.mesh(),
            this->psi_.dimensions()
        )
    );

    volScalarField& psib = tpsib();
    scalarField& psibCells = psib.internalField();
    const volScalarField Tb_(Tb());
    const scalarField& TbCells = Tb_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(psibCells, celli)
    {
        psibCells[celli] =
            this->cellReactants(celli).psi(pCells[celli], TbCells[celli]);
    }

    forAll(psib.boundaryField(), patchi)
    {
        fvPatchScalarField& ppsib = psib.boundaryField()[patchi];

        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        forAll(ppsib, facei)
        {
            ppsib[facei] =
                this->patchFaceReactants
                (patchi, facei).psi(pp[facei], pTb[facei]);
        }
    }

    return tpsib;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::muu() const
{
    tmp<volScalarField> tmuu
    (
        new volScalarField
        (
            IOobject
            (
                "muu",
                this->T_.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& muu_ = tmuu();
    scalarField& muuCells = muu_.internalField();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TuCells = this->Tu_.internalField();

    forAll(muuCells, celli)
    {
        muuCells[celli] = this->cellReactants(celli).mu
        (
            pCells[celli],
            TuCells[celli]
        );
    }

    forAll(muu_.boundaryField(), patchi)
    {
        fvPatchScalarField& pMuu = muu_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = this->Tu_.boundaryField()[patchi];

        forAll(pMuu, facei)
        {
            pMuu[facei] = this->patchFaceReactants(patchi, facei).mu
            (
                pp[facei],
                pTu[facei]
            );
        }
    }

    return tmuu;
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::mub() const
{
    tmp<volScalarField> tmub
    (
        new volScalarField
        (
            IOobject
            (
                "mub",
                this->T_.time().timeName(),
                this->T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& mub_ = tmub();
    scalarField& mubCells = mub_.internalField();
    const volScalarField Tb_(Tb());
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TbCells = Tb_.internalField();

    forAll(mubCells, celli)
    {
        mubCells[celli] = this->cellProducts(celli).mu
        (
            pCells[celli],
            TbCells[celli]
        );
    }

    forAll(mub_.boundaryField(), patchi)
    {
        fvPatchScalarField& pMub = mub_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        forAll(pMub, facei)
        {
            pMub[facei] = this->patchFaceProducts(patchi, facei).mu
            (
                pp[facei],
                pTb[facei]
            );
        }
    }

    return tmub;
}


// ************************************************************************* //
