/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "SLTSDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void SLTSDdtScheme<Type>::relaxedDiag
(
    scalarField& rD,
    const surfaceScalarField& phi
) const
{
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();
    scalarField diag(rD.size(), 0.0);

    forAll(owner, faceI)
    {
        if (phi[faceI] > 0.0)
        {
            diag[owner[faceI]] += phi[faceI];
            rD[neighbour[faceI]] += phi[faceI];
        }
        else
        {
            diag[neighbour[faceI]] -= phi[faceI];
            rD[owner[faceI]] -= phi[faceI];
        }
    }

    forAll(phi.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pphi = phi.boundaryField()[patchi];
        const labelUList& faceCells = pphi.patch().patch().faceCells();

        forAll(pphi, patchFacei)
        {
            if (pphi[patchFacei] > 0.0)
            {
                diag[faceCells[patchFacei]] += pphi[patchFacei];
            }
            else
            {
                rD[faceCells[patchFacei]] -= pphi[patchFacei];
            }
        }
    }

    rD += (1.0/alpha_ - 2.0)*diag;
}


template<class Type>
tmp<volScalarField> SLTSDdtScheme<Type>::SLrDeltaT() const
{
    const surfaceScalarField& phi =
        mesh().objectRegistry::template
            lookupObject<surfaceScalarField>(phiName_);

    const dimensionedScalar& deltaT = mesh().time().deltaT();

    tmp<volScalarField> trDeltaT
    (
        new volScalarField
        (
            IOobject
            (
                "rDeltaT",
                phi.instance(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("rDeltaT", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& rDeltaT = trDeltaT();

    relaxedDiag(rDeltaT, phi);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        rDeltaT.internalField() = max
        (
            rDeltaT.internalField()/mesh().V(),
            scalar(1)/deltaT.value()
        );
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
            mesh().objectRegistry::template lookupObject<volScalarField>
            (
                rhoName_
            ).oldTime();

        rDeltaT.internalField() = max
        (
            rDeltaT.internalField()/(rho.internalField()*mesh().V()),
            scalar(1)/deltaT.value()
        );
    }
    else
    {
        FatalErrorIn("SLTSDdtScheme<Type>::CofrDeltaT() const")
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);
    }

    rDeltaT.correctBoundaryConditions();

    return trDeltaT;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt().internalField() =
            rDeltaT.internalField()*dt.value()*(1.0 - mesh().V0()/mesh().V());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*rho.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()
               *alpha.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    alpha.internalField()
                   *rho.internalField()
                   *vf.internalField()

                  - alpha.oldTime().internalField()
                   *rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.boundaryField()*
                (
                    alpha.boundaryField()
                   *rho.boundaryField()
                   *vf.boundaryField()

                  - alpha.oldTime().boundaryField()
                   *rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT
               *(
                   alpha*rho*vf
                 - alpha.oldTime()*rho.oldTime()*vf.oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
SLTSDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT(SLrDeltaT()().internalField());

    Info<< "SLTSDdtScheme<Type>::fvmDdt: max/min rDeltaT "
        << gMax(rDeltaT) << " " << gMin(rDeltaT) << endl;

    fvm.diag() = rDeltaT*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
SLTSDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT(SLrDeltaT()().internalField());

    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
SLTSDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT(SLrDeltaT()().internalField());

    fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
SLTSDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT(SLrDeltaT()().internalField());

    fvm.diag() = rDeltaT*alpha.internalField()*rho.internalField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().internalField()
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().internalField()
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    IOobject ddtIOobject
    (
        "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    const surfaceScalarField rDeltaT(fvc::interpolate(SLrDeltaT()));

    fluxFieldType phiCorr
    (
        mesh().Sf() & (Uf.oldTime() - fvc::interpolate(U.oldTime()))
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            this->fvcDdtPhiCoeff
            (
                U.oldTime(),
                (mesh().Sf() & Uf.oldTime()),
                phiCorr
            )
           *rDeltaT*phiCorr
        )
    );
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtCorr(" + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    const surfaceScalarField rDeltaT(fvc::interpolate(SLrDeltaT()));

    fluxFieldType phiCorr
    (
        phi.oldTime() - (mesh().Sf() & fvc::interpolate(U.oldTime()))
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
           *rDeltaT*phiCorr
        )
    );
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    IOobject ddtIOobject
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    const surfaceScalarField rDeltaT(fvc::interpolate(SLrDeltaT()));

    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiCorr
        (
            mesh().Sf() & (Uf.oldTime() - fvc::interpolate(rhoU0))
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff
                (
                    rhoU0,
                    mesh().Sf() & Uf.oldTime(),
                    phiCorr
                )
               *rDeltaT*phiCorr
            )
        );
    }
    else if
    (
        U.dimensions() == dimDensity*dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        fluxFieldType phiCorr
        (
             mesh().Sf() & (Uf.oldTime() - fvc::interpolate(U.oldTime()))
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    mesh().Sf() & Uf.oldTime(),
                    phiCorr
                )
               *rDeltaT*phiCorr
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "SLTSDdtScheme<Type>::fvcDdtPhiCorr"
        )   << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiCorr
        (
            phi.oldTime() - (mesh().Sf() & fvc::interpolate(rhoU0))
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(rhoU0, phi.oldTime(), phiCorr)
               *rDeltaT*phiCorr
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        fluxFieldType phiCorr
        (
            phi.oldTime() - (mesh().Sf() & fvc::interpolate(U.oldTime()))
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
               *rDeltaT*phiCorr
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "SLTSDdtScheme<Type>::fvcDdtPhiCorr"
        )   << "dimensions of phi are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<surfaceScalarField> SLTSDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("0", dimVolume/dimTime, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
