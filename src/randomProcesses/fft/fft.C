/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "fft.H"
#include "fftRenumber.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#define TWOPI 6.28318530717959

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void fft::transform
(
    complexField& field,
    const labelList& nn,
    transformDirection isign
)
{
    forAll(nn, idim)
    {
        // Check for power of two
        unsigned int dimCount = nn[idim];
        if (!dimCount || (dimCount & (dimCount - 1)))
        {
            FatalErrorIn
            (
                 "fft::transform(complexField&, const labelList&, "
                 "transformDirection)"
            )   << "number of elements in direction " << idim
                << " is not a power of 2" << endl
                << "    Number of elements in each direction = " << nn
                << abort(FatalError);
        }
    }

    const label ndim = nn.size();

    label i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
    label ibit, k1, k2, n, nprev, nrem, idim;
    scalar tempi, tempr;
    scalar theta, wi, wpi, wpr, wr, wtemp;
    scalar* data = reinterpret_cast<scalar*>(field.begin()) - 1;


    // if inverse transform : renumber before transform

    if (isign == REVERSE_TRANSFORM)
    {
        fftRenumber(field, nn);
    }


    label ntot = 1;
    forAll(nn, idim)
    {
        ntot *= nn[idim];
    }


    nprev = 1;

    for (idim=ndim; idim>=1; idim--)
    {
        n = nn[idim-1];
        nrem = ntot/(n*nprev);
        ip1 = nprev << 1;
        ip2 = ip1*n;
        ip3 = ip2*nrem;
        i2rev = 1;

        for (i2=1; i2<=ip2; i2+=ip1)
        {
            if (i2 < i2rev)
            {
                for (i1=i2; i1<=i2 + ip1 - 2; i1+=2)
                {
                    for (i3=i1; i3<=ip3; i3+=ip2)
                    {
                        i3rev = i2rev + i3 - i2;
                        SWAP(data[i3], data[i3rev]);
                        SWAP(data[i3 + 1], data[i3rev + 1]);
                    }
                }
            }

            ibit = ip2 >> 1;
            while (ibit >= ip1 && i2rev > ibit)
            {
                i2rev -= ibit;
                ibit >>= 1;
            }

            i2rev += ibit;
        }

        ifp1 = ip1;

        while (ifp1 < ip2)
        {
            ifp2 = ifp1 << 1;
            theta = isign*TWOPI/(ifp2/ip1);
            wtemp = sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;

            for (i3 = 1; i3 <= ifp1; i3 += ip1)
            {
                for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2)
                {
                    for (i2 = i1; i2 <= ip3; i2 += ifp2)
                    {
                        k1 = i2;
                        k2 = k1 + ifp1;
                        tempr = scalar(wr*data[k2]) - scalar(wi*data[k2 + 1]);
                        tempi = scalar(wr*data[k2 + 1]) + scalar(wi*data[k2]);
                        data[k2] = data[k1] - tempr;
                        data[k2 + 1] = data[k1 + 1] - tempi;
                        data[k1] += tempr;
                        data[k1 + 1] += tempi;
                    }
                }

                wr = (wtemp = wr)*wpr - wi*wpi + wr;
                wi = wi*wpr + wtemp*wpi + wi;
            }
            ifp1 = ifp2;
        }
        nprev *= n;
    }


    // if forward transform : renumber after transform

    if (isign == FORWARD_TRANSFORM)
    {
        fftRenumber(field, nn);
    }


    // scale result (symmetric scale both forward and inverse transform)

    scalar recRootN = 1.0/sqrt(scalar(ntot));

    forAll(field, i)
    {
        field[i] *= recRootN;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef SWAP
#undef TWOPI

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<complexField> fft::forwardTransform
(
    const tmp<complexField>& tfield,
    const labelList& nn
)
{
    tmp<complexField> tfftField(new complexField(tfield));

    transform(tfftField(), nn, FORWARD_TRANSFORM);

    tfield.clear();

    return tfftField;
}


tmp<complexField> fft::reverseTransform
(
    const tmp<complexField>& tfield,
    const labelList& nn
)
{
    tmp<complexField> tifftField(new complexField(tfield));

    transform(tifftField(), nn, REVERSE_TRANSFORM);

    tfield.clear();

    return tifftField;
}


tmp<complexVectorField> fft::forwardTransform
(
    const tmp<complexVectorField>& tfield,
    const labelList& nn
)
{
    tmp<complexVectorField> tfftVectorField
    (
        new complexVectorField
        (
            tfield().size()
        )
    );

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        tfftVectorField().replace
        (
            cmpt,
            forwardTransform(tfield().component(cmpt), nn)
        );
    }

    tfield.clear();

    return tfftVectorField;
}


tmp<complexVectorField> fft::reverseTransform
(
    const tmp<complexVectorField>& tfield,
    const labelList& nn
)
{
    tmp<complexVectorField> tifftVectorField
    (
        new complexVectorField
        (
            tfield().size()
        )
    );

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        tifftVectorField().replace
        (
            cmpt,
            reverseTransform(tfield().component(cmpt), nn)
        );
    }

    tfield.clear();

    return tifftVectorField;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
