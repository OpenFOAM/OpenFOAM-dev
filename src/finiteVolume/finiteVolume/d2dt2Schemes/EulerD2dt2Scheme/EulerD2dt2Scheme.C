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

#include "EulerD2dt2Scheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
EulerD2dt2Scheme<Type>::fvcD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT2 =
        4.0/sqr(mesh().time().deltaT() + mesh().time().deltaT0());

    IOobject d2dt2IOobject
    (
        "d2dt2("+vf.name()+')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = mesh().time().deltaTValue();
    scalar deltaT0 = mesh().time().deltaT0Value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = rDeltaT2.value()/2.0;

        scalarField VV0 = mesh().V() + mesh().V0();
        scalarField V0V00 = mesh().V0() + mesh().V00();

        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                mesh(),
                rDeltaT2.dimensions()*vf.dimensions(),
                halfRdeltaT2*
                (
                    coefft*VV0*vf.primitiveField()

                  - (coefft*VV0 + coefft00*V0V00)
                   *vf.oldTime().primitiveField()

                  + (coefft00*V0V00)*vf.oldTime().oldTime().primitiveField()
                )/mesh().V(),
                rDeltaT2.value()*
                (
                    coefft*vf.boundaryField()
                  - coefft0*vf.oldTime().boundaryField()
                  + coefft00*vf.oldTime().oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                rDeltaT2*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
EulerD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT2 =
        4.0/sqr(mesh().time().deltaT() + mesh().time().deltaT0());

    IOobject d2dt2IOobject
    (
        "d2dt2("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = mesh().time().deltaTValue();
    scalar deltaT0 = mesh().time().deltaT0Value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2.value();
        scalar quarterRdeltaT2 = 0.25*rDeltaT2.value();

        const scalarField VV0rhoRho0
        (
            (mesh().V() + mesh().V0())
          * (rho.primitiveField() + rho.oldTime().primitiveField())
        );

        const scalarField V0V00rho0Rho00
        (
            (mesh().V0() + mesh().V00())
          * (
                rho.oldTime().primitiveField()
              + rho.oldTime().oldTime().primitiveField()
            )
        );

        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                mesh(),
                rDeltaT2.dimensions()*rho.dimensions()*vf.dimensions(),
                quarterRdeltaT2*
                (
                    coefft*VV0rhoRho0*vf.primitiveField()

                  - (coefft*VV0rhoRho0 + coefft00*V0V00rho0Rho00)
                   *vf.oldTime().primitiveField()

                  + (coefft00*V0V00rho0Rho00)
                   *vf.oldTime().oldTime().primitiveField()
                )/mesh().V(),
                halfRdeltaT2*
                (
                    coefft
                   *(rho.boundaryField() + rho.oldTime().boundaryField())
                   *vf.boundaryField()

                  - (
                        coefft
                       *(
                           rho.boundaryField()
                         + rho.oldTime().boundaryField()
                        )
                      + coefft00
                       *(
                           rho.oldTime().boundaryField()
                         + rho.oldTime().oldTime().boundaryField()
                        )
                    )*vf.oldTime().boundaryField()

                  + coefft00
                   *(
                       rho.oldTime().boundaryField()
                     + rho.oldTime().oldTime().boundaryField()
                    )*vf.oldTime().oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        dimensionedScalar halfRdeltaT2 = 0.5*rDeltaT2;

        const volScalarField rhoRho0(rho + rho.oldTime());
        const volScalarField rho0Rho00(rho.oldTime() +rho.oldTime().oldTime());

        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                halfRdeltaT2*
                (
                    coefft*rhoRho0*vf
                  - (coefft*rhoRho0 + coefft00*rho0Rho00)*vf.oldTime()
                  + coefft00*rho0Rho00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type>>
EulerD2dt2Scheme<Type>::fvmD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar deltaT = mesh().time().deltaTValue();
    scalar deltaT0 = mesh().time().deltaT0Value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);
    scalar coefft0  = coefft + coefft00;

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = rDeltaT2/2.0;

        const scalarField VV0(mesh().V() + mesh().V0());
        const scalarField V0V00(mesh().V0() + mesh().V00());

        fvm.diag() = (coefft*halfRdeltaT2)*VV0;

        fvm.source() = halfRdeltaT2*
        (
            (coefft*VV0 + coefft00*V0V00)
           *vf.oldTime().primitiveField()

          - (coefft00*V0V00)*vf.oldTime().oldTime().primitiveField()
        );
    }
    else
    {
        fvm.diag() = (coefft*rDeltaT2)*mesh().V();

        fvm.source() = rDeltaT2*mesh().V()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
EulerD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol
            /dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar deltaT = mesh().time().deltaTValue();
    scalar deltaT0 = mesh().time().deltaT0Value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2;

        const scalarField VV0(mesh().V() + mesh().V0());
        const scalarField V0V00(mesh().V0() + mesh().V00());

        fvm.diag() = rho.value()*(coefft*halfRdeltaT2)*VV0;

        fvm.source() = halfRdeltaT2*rho.value()*
        (
            (coefft*VV0 + coefft00*V0V00)
           *vf.oldTime().primitiveField()

          - (coefft00*V0V00)*vf.oldTime().oldTime().primitiveField()
        );
    }
    else
    {
        fvm.diag() = (coefft*rDeltaT2)*mesh().V()*rho.value();

        fvm.source() = rDeltaT2*mesh().V()*rho.value()*
        (
            (coefft + coefft00)*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
EulerD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol
            /dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar deltaT = mesh().time().deltaTValue();
    scalar deltaT0 = mesh().time().deltaT0Value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar quarterRdeltaT2 = 0.25*rDeltaT2;

        const scalarField VV0rhoRho0
        (
            (mesh().V() + mesh().V0())
           *(rho.primitiveField() + rho.oldTime().primitiveField())
        );

        const scalarField V0V00rho0Rho00
        (
            (mesh().V0() + mesh().V00())
           *(
               rho.oldTime().primitiveField()
             + rho.oldTime().oldTime().primitiveField()
            )
        );

        fvm.diag() = (coefft*quarterRdeltaT2)*VV0rhoRho0;

        fvm.source() = quarterRdeltaT2*
        (
            (coefft*VV0rhoRho0 + coefft00*V0V00rho0Rho00)
           *vf.oldTime().primitiveField()

          - (coefft00*V0V00rho0Rho00)
           *vf.oldTime().oldTime().primitiveField()
        );
    }
    else
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2;

        const scalarField rhoRho0
        (
            rho.primitiveField()
          + rho.oldTime().primitiveField()
        );

        const scalarField rho0Rho00
        (
            rho.oldTime().primitiveField()
          + rho.oldTime().oldTime().primitiveField()
        );

        fvm.diag() = (coefft*halfRdeltaT2)*mesh().V()*rhoRho0;

        fvm.source() = halfRdeltaT2*mesh().V()*
        (
            (coefft*rhoRho0 + coefft00*rho0Rho00)
           *vf.oldTime().primitiveField()

          - (coefft00*rho0Rho00)
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
