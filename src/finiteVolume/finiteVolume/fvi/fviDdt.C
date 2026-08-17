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

#include "fviDdt.H"
#include "fvMesh.H"
#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvi
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolInternalField<Type>>
ddt(const dimensioned<Type> dt, const fvMesh& mesh)
{
    return fv::ddtScheme<Type>::New
    (
        mesh,
        mesh.schemes().ddt("ddt(" + dt.name() + ')')
    ).ref().fviDdt(dt);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt(const VolInternalField<Type>& vf)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + vf.name() + ')')
    ).ref().fviDdt(vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt(const VolField<Type>& vf)
{
    return fvi::ddt(vf());
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const dimensionedScalar& rho,
    const VolInternalField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().fviDdt(rho, vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    return fvi::ddt(rho, vf());
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const volInternalScalarField& rho,
    const VolInternalField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().fviDdt(rho, vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fvi::ddt(rho(), vf());
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const one&,
    const VolInternalField<Type>& vf
)
{
    return fvi::ddt(vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const one&,
    const VolField<Type>& vf
)
{
    return fvi::ddt(vf());
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const volInternalScalarField& alpha,
    const volInternalScalarField& rho,
    const VolInternalField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt
        (
            "ddt("
          + alpha.name() + ','
          + rho.name() + ','
          + vf.name() + ')'
        )
    ).ref().fviDdt(alpha, rho, vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fvi::ddt(alpha(), rho(), vf());
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const one&,
    const one&,
    const VolInternalField<Type>& vf
)
{
    return fvi::ddt(vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const one&,
    const one&,
    const VolField<Type>& vf
)
{
    return fvi::ddt(vf());
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const one&,
    const volInternalScalarField& rho,
    const VolInternalField<Type>& vf
)
{
    return fvi::ddt(rho, vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const one&,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fvi::ddt(rho(), vf());
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const volInternalScalarField& alpha,
    const one&,
    const VolInternalField<Type>& vf
)
{
    return fvi::ddt(alpha, vf);
}


template<class Type>
tmp<VolInternalField<Type>>
ddt
(
    const volScalarField& alpha,
    const one&,
    const VolField<Type>& vf
)
{
    return fvi::ddt(alpha(), vf());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
