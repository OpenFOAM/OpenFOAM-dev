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

#include "momentumSurfaceFilm.H"
#include "surfaceInterpolate.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::dimensionedScalar&
Foam::momentumSurfaceFilm::deltaSmall() const
{
    return deltaSmall_;
}


inline const Foam::word&
Foam::momentumSurfaceFilm::phaseName() const
{
    return phaseName_;
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::rho() const
{
    return thermo_->rho();
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::mu() const
{
    return mu_;
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::delta() const
{
    return delta_;
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::alpha() const
{
    return alpha_;
}


inline const Foam::volVectorField&
Foam::momentumSurfaceFilm::U() const
{
    return U_;
}


inline Foam::tmp<Foam::volVectorField::Internal>
Foam::momentumSurfaceFilm::Us() const
{
    return momentumTransport_->Us();
}


inline const Foam::surfaceScalarField&
Foam::momentumSurfaceFilm::phi() const
{
    return phi_;
}


inline const Foam::surfaceScalarField&
Foam::momentumSurfaceFilm::phiU() const
{
    return phiU_;
}


inline const Foam::volScalarField::Internal&
Foam::momentumSurfaceFilm::continuityErr() const
{
    return continuityErr_;
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::coverage() const
{
    return coverage_;
}


inline Foam::volVectorField&
Foam::momentumSurfaceFilm::USpPrimary()
{
    return USpPrimary_;
}


inline Foam::volScalarField&
Foam::momentumSurfaceFilm::pSpPrimary()
{
    return pSpPrimary_;
}


inline Foam::volScalarField&
Foam::momentumSurfaceFilm::rhoSpPrimary()
{
    return rhoSpPrimary_;
}


inline Foam::volVectorField::Internal&
Foam::momentumSurfaceFilm::USp()
{
    return USp_;
}


inline Foam::volScalarField::Internal&
Foam::momentumSurfaceFilm::pSp()
{
    return pSp_;
}


inline Foam::volScalarField::Internal&
Foam::momentumSurfaceFilm::rhoSp()
{
    return rhoSp_;
}


inline const Foam::volVectorField::Internal&
Foam::momentumSurfaceFilm::USp() const
{
    return USp_;
}


inline const Foam::volScalarField::Internal&
Foam::momentumSurfaceFilm::pSp() const
{
    return pSp_;
}


inline const Foam::volScalarField::Internal&
Foam::momentumSurfaceFilm::rhoSp() const
{
    return rhoSp_;
}


inline const Foam::volVectorField&
Foam::momentumSurfaceFilm::UPrimary() const
{
    return UPrimary_;
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::pPrimary() const
{
    return thermo_->p();
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::rhoPrimary() const
{
    return rhoPrimary_;
}


inline const Foam::volScalarField&
Foam::momentumSurfaceFilm::muPrimary() const
{
    return muPrimary_;
}


inline const Foam::rhoThermo&
Foam::momentumSurfaceFilm::thermo() const
{
    return thermo_();
}


inline Foam::surfaceFilmModels::ejectionModelList&
Foam::momentumSurfaceFilm::ejection()
{
    return ejection_;
}


inline Foam::surfaceFilmModels::transferModelList&
Foam::momentumSurfaceFilm::transfer()
{
    return transfer_;
}


inline const Foam::surfaceFilmModels::momentumTransportModel&
Foam::momentumSurfaceFilm::momentumTransport() const
{
    return momentumTransport_();
}


inline Foam::tmp<Foam::volScalarField::Internal>
Foam::momentumSurfaceFilm::mass() const
{
    return rho()()*delta_()*magSf();
}


inline Foam::tmp<Foam::volScalarField::Internal>
Foam::momentumSurfaceFilm::deltaMass() const
{
    return rhoSp_*magSf()*time().deltaT();
}


// ************************************************************************* //
