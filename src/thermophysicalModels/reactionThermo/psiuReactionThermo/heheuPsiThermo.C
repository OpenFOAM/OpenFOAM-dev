/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    const scalarField& hCells = this->he_;
    const scalarField& heuCells = this->heu_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& TuCells = this->Tu_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

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

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& TuBf =
        this->Tu_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& heuBf =
        this->heu().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pTu = TuBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pheu = heuBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);

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
    scalarField& heuCells = this->heu_.primitiveFieldRef();
    const scalarField& pCells = this->p_;
    const scalarField& TuCells = this->Tu_;

    forAll(heuCells, celli)
    {
        heuCells[celli] = this->cellReactants(celli).HE
        (
            pCells[celli],
            TuCells[celli]
        );
    }

    volScalarField::Boundary& heuBf = heu_.boundaryFieldRef();

    forAll(heuBf, patchi)
    {
        fvPatchScalarField& pheu = heuBf[patchi];
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
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
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
    scalarField& heu = theu.ref();

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
    scalarField& heu = theu.ref();

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
                IOobject::NO_WRITE,
                false
            ),
            this->T_
        )
    );

    volScalarField& Tb_ = tTb.ref();
    scalarField& TbCells = Tb_.primitiveFieldRef();
    const scalarField& pCells = this->p_;
    const scalarField& TCells = this->T_;
    const scalarField& hCells = this->he_;

    forAll(TbCells, celli)
    {
        TbCells[celli] = this->cellProducts(celli).THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );
    }

    volScalarField::Boundary& TbBf = Tb_.boundaryFieldRef();

    forAll(TbBf, patchi)
    {
        fvPatchScalarField& pTb = TbBf[patchi];

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
                IOobject::NO_WRITE,
                false
            ),
            this->psi_.mesh(),
            this->psi_.dimensions()
        )
    );

    volScalarField& psiu = tpsiu.ref();
    scalarField& psiuCells = psiu.primitiveFieldRef();
    const scalarField& TuCells = this->Tu_;
    const scalarField& pCells = this->p_;

    forAll(psiuCells, celli)
    {
        psiuCells[celli] =
            this->cellReactants(celli).psi(pCells[celli], TuCells[celli]);
    }

    volScalarField::Boundary& psiuBf = psiu.boundaryFieldRef();

    forAll(psiuBf, patchi)
    {
        fvPatchScalarField& ppsiu = psiuBf[patchi];

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
                IOobject::NO_WRITE,
                false
            ),
            this->psi_.mesh(),
            this->psi_.dimensions()
        )
    );

    volScalarField& psib = tpsib.ref();
    scalarField& psibCells = psib.primitiveFieldRef();
    const volScalarField Tb_(Tb());
    const scalarField& TbCells = Tb_;
    const scalarField& pCells = this->p_;

    forAll(psibCells, celli)
    {
        psibCells[celli] =
            this->cellReactants(celli).psi(pCells[celli], TbCells[celli]);
    }

    volScalarField::Boundary& psibBf = psib.boundaryFieldRef();

    forAll(psibBf, patchi)
    {
        fvPatchScalarField& ppsib = psibBf[patchi];

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
                IOobject::NO_WRITE,
                false
            ),
            this->T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& muu_ = tmuu.ref();
    scalarField& muuCells = muu_.primitiveFieldRef();
    const scalarField& pCells = this->p_;
    const scalarField& TuCells = this->Tu_;

    forAll(muuCells, celli)
    {
        muuCells[celli] = this->cellReactants(celli).mu
        (
            pCells[celli],
            TuCells[celli]
        );
    }

    volScalarField::Boundary& muuBf = muu_.boundaryFieldRef();

    forAll(muuBf, patchi)
    {
        fvPatchScalarField& pMuu = muuBf[patchi];
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
                IOobject::NO_WRITE,
                false
            ),
            this->T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& mub_ = tmub.ref();
    scalarField& mubCells = mub_.primitiveFieldRef();
    const volScalarField Tb_(Tb());
    const scalarField& pCells = this->p_;
    const scalarField& TbCells = Tb_;

    forAll(mubCells, celli)
    {
        mubCells[celli] = this->cellProducts(celli).mu
        (
            pCells[celli],
            TbCells[celli]
        );
    }

    volScalarField::Boundary& mubBf = mub_.boundaryFieldRef();

    forAll(mubBf, patchi)
    {
        fvPatchScalarField& pMub = mubBf[patchi];
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
