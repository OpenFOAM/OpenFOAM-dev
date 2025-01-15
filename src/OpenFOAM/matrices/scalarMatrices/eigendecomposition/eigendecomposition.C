/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "eigendecomposition.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::eigendecomposition::tred2()
{
    // This is derived from the Algol procedures tred2 by
    // Bowdler, Martin, Reinsch, and Wilkinson,
    // Handbook for Automatic Computation: Volume II: Linear Algebra

    const label n = V_.n();

    for (label j = 0; j < n; j++)
    {
        d_[j] = V_(n-1, j);
    }

    // Householder reduction to tridiagonal form

    for (label i = n-1; i > 0; i--)
    {
        // Scale to avoid under/overflow

        scalar scale = 0;
        scalar h = 0;

        for (label k = 0; k < i; k++)
        {
            scale = scale + mag(d_[k]);
        }

        if (scale == 0)
        {
            e_[i] = d_[i-1];

            for (label j = 0; j < i; j++)
            {
                d_[j] = V_(i-1, j);
                V_(i, j) = 0;
                V_(j, i) = 0;
            }
        }
        else
        {
            // Generate Householder vector

            for (label k = 0; k < i; k++)
            {
                d_[k] /= scale;
                h += d_[k]*d_[k];
            }

            scalar f = d_[i-1];
            scalar g = sqrt(h);

            if (f > 0)
            {
                g = -g;
            }

            e_[i] = scale*g;
            h = h - f*g;
            d_[i-1] = f - g;

            for (label j = 0; j < i; j++)
            {
                e_[j] = 0;
            }

            // Apply similarity transformation to remaining columns

            for (label j = 0; j < i; j++)
            {
                f = d_[j];
                V_(j, i) = f;
                g = e_[j] + V_(j, j)*f;
                for (label k = j+1; k <= i-1; k++)
                {
                    g += V_(k, j)*d_[k];
                    e_[k] += V_(k, j)*f;
                }
                e_[j] = g;
            }

            f = 0;

            for (label j = 0; j < i; j++)
            {
                e_[j] /= h;
                f += e_[j]*d_[j];
            }

            const scalar hh = f/(h + h);

            for (label j = 0; j < i; j++)
            {
                e_[j] -= hh*d_[j];
            }

            for (label j = 0; j < i; j++)
            {
                f = d_[j];
                g = e_[j];

                for (label k = j; k <= i-1; k++)
                {
                    V_(k, j) -= (f*e_[k] + g*d_[k]);
                }

                d_[j] = V_(i-1, j);
                V_(i, j) = 0;
            }
        }
        d_[i] = h;
    }

    // Accumulate transformations

    for (label i = 0; i < n-1; i++)
    {
        V_(n-1, i) = V_(i, i);
        V_(i, i) = 1;

        const scalar h = d_[i+1];

        if (h != 0)
        {
            for (label k = 0; k <= i; k++)
            {
                d_[k] = V_(k, i+1)/h;
            }

            for (label j = 0; j <= i; j++)
            {
                scalar g = 0;

                for (label k = 0; k <= i; k++)
                {
                    g += V_(k, i+1)*V_(k, j);
                }

                for (label k = 0; k <= i; k++)
                {
                    V_(k, j) -= g*d_[k];
                }
            }
        }

        for (label k = 0; k <= i; k++)
        {
            V_(k, i+1) = 0;
        }
    }

    for (label j = 0; j < n; j++)
    {
        d_[j] = V_(n-1, j);
        V_(n-1, j) = 0;
    }

    V_(n-1, n-1) = 1;
    e_[0] = 0;
}


void Foam::eigendecomposition::tql2()
{
    // This is derived from the Algol procedures tql2, by
    // Bowdler, Martin, Reinsch, and Wilkinson,
    // Handbook for Automatic Computation: Volume II: Linear Algebra

    const label n = V_.n();

    for (label i = 1; i < n; i++)
    {
        e_[i-1] = e_[i];
    }
    e_[n-1] = 0;

    scalar f = 0;
    scalar tst1 = 0;

    for (label l = 0; l < n; l++)
    {
        // Find small subdiagonal element

        tst1 = max(tst1, mag(d_[l]) + mag(e_[l]));

        label m = l;

        // Original while-loop from Java code
        while (m < n)
        {
            if (mag(e_[m]) <= small*tst1)
            {
                break;
            }
            m++;
        }

        // If m == l, d_[l] is an eigenvalue,
        // otherwise, iterate

        if (m > l)
        {
            label iter = 0;

            do
            {
                iter = iter + 1;  // Could check iteration count here

                // Compute implicit shift

                scalar g = d_[l];
                scalar p = (d_[l+1] - g)/(2*e_[l]);
                scalar r = hypot(p, 1);

                if (p < 0)
                {
                    r = -r;
                }

                d_[l] = e_[l]/(p + r);
                d_[l+1] = e_[l]*(p + r);
                const scalar dl1 = d_[l+1];
                scalar h = g - d_[l];

                for (label i = l+2; i < n; i++)
                {
                    d_[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation

                p = d_[m];
                scalar c = 1;
                scalar c2 = c;
                scalar c3 = c;
                scalar el1 = e_[l+1];
                scalar s = 0;
                scalar s2 = 0;

                for (label i = m-1; i >= l; i--)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c*e_[i];
                    h = c*p;
                    r = hypot(p, e_[i]);
                    e_[i+1] = s*r;
                    s = e_[i]/r;
                    c = p/r;
                    p = c*d_[i] - s*g;
                    d_[i+1] = h + s*(c*g + s*d_[i]);

                    // Accumulate transformation

                    for (label k = 0; k < n; k++)
                    {
                        h = V_(k, i+1);
                        V_(k, i+1) = s*V_(k, i) + c*h;
                        V_(k, i) = c*V_(k, i) - s*h;
                    }
                }

                p = -s*s2*c3*el1*e_[l]/dl1;
                e_[l] = s*p;
                d_[l] = c*p;

                // Check for convergence

            } while (mag(e_[l]) > small*tst1);
        }
        d_[l] = d_[l] + f;
        e_[l] = 0;
    }

    // Sort eigenvalues and corresponding vectors

    for (label i = 0; i < n-1; i++)
    {
        label k = i;
        scalar p = d_[i];

        for (label j = i+1; j < n; j++)
        {
            if (d_[j] < p)
            {
                k = j;
                p = d_[j];
            }
        }

        if (k != i)
        {
            d_[k] = d_[i];
            d_[i] = p;

            for (label j = 0; j < n; j++)
            {
                p = V_(j, i);
                V_(j, i) = V_(j, k);
                V_(j, k) = p;
            }
        }
    }
}


void Foam::eigendecomposition::orthes()
{
    // This is derived from the Algol procedures orthes and ortran,
    // by Martin and Wilkinson,
    // Handbook for Automatic Computation: Volume II: Linear Algebra

    const label n = V_.n();
    const label low = 0;
    const label high = n-1;

    for (label m = low+1; m <= high-1; m++)
    {
        // Scale column

        scalar scale = 0;
        for (label i = m; i <= high; i++)
        {
            scale = scale + mag(H_(i, m-1));
        }

        if (scale != 0)
        {
            // Compute Householder transformation

            scalar h = 0;
            for (label i = high; i >= m; i--)
            {
                ort_[i] = H_(i, m-1)/scale;
                h += ort_[i]*ort_[i];
            }

            scalar g = sqrt(h);

            if (ort_[m] > 0)
            {
                g = -g;
            }

            h = h - ort_[m]*g;
            ort_[m] = ort_[m] - g;

            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)

            for (label j = m; j < n; j++)
            {
                scalar f = 0;
                for (label i = high; i >= m; i--)
                {
                    f += ort_[i]*H_(i, j);
                }

                f = f/h;

                for (label i = m; i <= high; i++)
                {
                    H_(i, j) -= f*ort_[i];
                }
            }

            for (label i = 0; i <= high; i++)
            {
                scalar f = 0;

                for (label j = high; j >= m; j--)
                {
                    f += ort_[j]*H_(i, j);
                }

                f = f/h;

                for (label j = m; j <= high; j++)
                {
                    H_(i, j) -= f*ort_[j];
                }
            }

            ort_[m] = scale*ort_[m];
            H_(m, m-1) = scale*g;
        }
    }

    // Accumulate transformations (Algol's ortran)

    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            V_(i, j) = (i == j ? 1 : 0);
        }
    }

    for (label m = high-1; m >= low+1; m--)
    {
        if (H_(m, m-1) != 0)
        {
            for (label i = m+1; i <= high; i++)
            {
                ort_[i] = H_(i, m-1);
            }

            for (label j = m; j <= high; j++)
            {
                scalar g = 0;

                for (label i = m; i <= high; i++)
                {
                    g += ort_[i]*V_(i, j);
                }

                // Double division avoids possible underflow
                g = (g/ort_[m])/H_(m, m-1);

                for (label i = m; i <= high; i++)
                {
                    V_(i, j) += g*ort_[i];
                }
            }
        }
    }
}


inline void Foam::eigendecomposition::cdiv
(
    scalar& cdivr,
    scalar& cdivi,
    const scalar xr,
    const scalar xi,
    const scalar yr,
    const scalar yi
)
{
    if (mag(yr) > mag(yi))
    {
        const scalar r = yi/yr;
        const scalar d = yr + r*yi;
        cdivr = (xr + r*xi)/d;
        cdivi = (xi - r*xr)/d;
    }
    else
    {
        const scalar r = yr/yi;
        const scalar d = yi + r*yr;
        cdivr = (r*xr + xi)/d;
        cdivi = (r*xi - xr)/d;
    }
}


void Foam::eigendecomposition::hqr2()
{
    // This is derived from the Algol procedure hqr2,
    // by Martin and Wilkinson,
    // Handbook for Automatic Computation: Volume II: Linear Algebra

    // Initialise

    const label nn = V_.n();

    label n = nn-1;
    const label low = 0;
    const label high = nn-1;

    scalar exshift = 0;
    scalar p = 0;
    scalar q = 0;
    scalar r = 0;
    scalar s = 0;
    scalar z = 0;
    scalar t, w, x, y;

    // Store roots isolated by balanc and compute matrix norm

    scalar norm = 0;
    for (label i = 0; i < nn; i++)
    {
        if ((i < low) || (i > high))
        {
            d_[i] = H_(i, i);
            e_[i] = 0;
        }
        for (label j = max(i-1, 0); j < nn; j++)
        {
            norm = norm + mag(H_(i, j));
        }
    }

    // Outer loop over eigenvalue index

    label iter = 0;
    while (n >= low)
    {
        // Look for single small sub-diagonal element

        label l = n;
        while (l > low)
        {
            s = mag(H_(l-1, l-1)) + mag(H_(l, l));
            if (s == 0)
            {
                s = norm;
            }
            if (mag(H_(l, l-1)) < small*s)
            {
                break;
            }
            l--;
        }

        // Check for convergence
        // One root found

        if (l == n)
        {
            H_(n, n) = H_(n, n) + exshift;
            d_[n] = H_(n, n);
            e_[n] = 0;
            n--;
            iter = 0;

            // Two roots found
        }
        else if (l == n-1)
        {
            w = H_(n, n-1)*H_(n-1, n);
            p = (H_(n-1, n-1) - H_(n, n))/2;
            q = p*p + w;
            z = sqrt(mag(q));
            H_(n, n) = H_(n, n) + exshift;
            H_(n-1, n-1) = H_(n-1, n-1) + exshift;
            x = H_(n, n);

            // scalar pair

            if (q >= 0)
            {
                if (p >= 0)
                {
                    z = p + z;
                }
                else
                {
                    z = p - z;
                }

                d_[n-1] = x + z;
                d_[n] = d_[n-1];

                if (z != 0)
                {
                    d_[n] = x - w/z;
                }

                e_[n-1] = 0;
                e_[n] = 0;
                x = H_(n, n-1);
                s = mag(x) + mag(z);
                p = x/s;
                q = z/s;
                r = sqrt(p*p+q*q);
                p = p/r;
                q = q/r;

                // Row modification
                for (label j = n-1; j < nn; j++)
                {
                    z = H_(n-1, j);
                    H_(n-1, j) = q*z + p*H_(n, j);
                    H_(n, j) = q*H_(n, j) - p*z;
                }

                // Column modification
                for (label i = 0; i <= n; i++)
                {
                    z = H_(i, n-1);
                    H_(i, n-1) = q*z + p*H_(i, n);
                    H_(i, n) = q*H_(i, n) - p*z;
                }

                // Accumulate transformations
                for (label i = low; i <= high; i++)
                {
                    z = V_(i, n-1);
                    V_(i, n-1) = q*z + p*V_(i, n);
                    V_(i, n) = q*V_(i, n) - p*z;
                }

                // Complex pair
            }
            else
            {
                d_[n-1] = x + p;
                d_[n] = x + p;
                e_[n-1] = z;
                e_[n] = -z;
            }

            n = n - 2;
            iter = 0;

            // No convergence yet

        }
        else
        {
            // Form shift

            x = H_(n, n);
            y = 0;
            w = 0;

            if (l < n)
            {
                y = H_(n-1, n-1);
                w = H_(n, n-1)*H_(n-1, n);
            }

            // Wilkinson's original ad hoc shift
            if (iter == 10)
            {
                exshift += x;
                for (label i = low; i <= n; i++)
                {
                    H_(i, i) -= x;
                }
                s = mag(H_(n, n-1)) + mag(H_(n-1, n-2));
                x = y = 0.75*s;
                w = -0.4375*s*s;
            }

            // MATLAB's new ad hoc shift
            if (iter == 30)
            {
                s = (y - x)/2;
                s = s*s + w;
                if (s > 0)
                {
                    s = sqrt(s);

                    if (y < x)
                    {
                        s = -s;
                    }

                    s = x - w/((y - x)/2 + s);

                    for (label i = low; i <= n; i++)
                    {
                        H_(i, i) -= s;
                    }

                    exshift += s;
                    x = y = w = 0.964;
                }
            }

            iter = iter + 1; // Could check iteration count here

            // Look for two consecutive small sub-diagonal elements

            label m = n-2;
            while (m >= l)
            {
                z = H_(m, m);
                r = x - z;
                s = y - z;
                p = (r*s - w)/H_(m+1, m) + H_(m, m+1);
                q = H_(m+1, m+1) - z - r - s;
                r = H_(m+2, m+1);
                s = mag(p) + mag(q) + mag(r);
                p = p/s;
                q = q/s;
                r = r/s;

                if (m == l)
                {
                    break;
                }

                if
                (
                    mag(H_(m, m-1))*(mag(q) + mag(r))
                  < small
                   *(mag(p)*(mag(H_(m-1, m-1)) + mag(z) + mag(H_(m+1, m+1))))
                )
                {
                    break;
                }

                m--;
            }

            for (label i = m+2; i <= n; i++)
            {
                H_(i, i-2) = 0;
                if (i > m+2)
                {
                    H_(i, i-3) = 0;
                }
            }

            // Double QR step involving rows l:n and columns m:n

            for (label k = m; k <= n-1; k++)
            {
                const label notlast = (k != n-1);

                if (k != m)
                {
                    p = H_(k, k-1);
                    q = H_(k+1, k-1);
                    r = (notlast ? H_(k+2, k-1) : 0);
                    x = mag(p) + mag(q) + mag(r);
                    if (x != 0)
                    {
                        p = p/x;
                        q = q/x;
                        r = r/x;
                    }
                }

                if (x == 0)
                {
                    break;
                }

                s = sqrt(p*p + q*q + r*r);

                if (p < 0)
                {
                    s = -s;
                }

                if (s != 0)
                {
                    if (k != m)
                    {
                        H_(k, k-1) = -s*x;
                    }
                    else if (l != m)
                    {
                        H_(k, k-1) = -H_(k, k-1);
                    }
                    p = p + s;
                    x = p/s;
                    y = q/s;
                    z = r/s;
                    q = q/p;
                    r = r/p;

                    // Row modification
                    for (label j = k; j < nn; j++)
                    {
                        p = H_(k, j) + q*H_(k+1, j);
                        if (notlast)
                        {
                            p = p + r*H_(k+2, j);
                            H_(k+2, j) = H_(k+2, j) - p*z;
                        }
                        H_(k, j) = H_(k, j) - p*x;
                        H_(k+1, j) = H_(k+1, j) - p*y;
                    }

                    // Column modification
                    for (label i = 0; i <= min(n, k+3); i++)
                    {
                        p = x*H_(i, k) + y*H_(i, k+1);
                        if (notlast)
                        {
                            p = p + z*H_(i, k+2);
                            H_(i, k+2) = H_(i, k+2) - p*r;
                        }
                        H_(i, k) = H_(i, k) - p;
                        H_(i, k+1) = H_(i, k+1) - p*q;
                    }

                    // Accumulate transformations
                    for (label i = low; i <= high; i++)
                    {
                        p = x*V_(i, k) + y*V_(i, k+1);
                        if (notlast)
                        {
                            p = p + z*V_(i, k+2);
                            V_(i, k+2) = V_(i, k+2) - p*r;
                        }
                        V_(i, k) = V_(i, k) - p;
                        V_(i, k+1) = V_(i, k+1) - p*q;
                    }
                }  // (s != 0)
            }  // k loop
        }  // check convergence
    }  // while (n >= low)

    // Backsubstitute to find vectors of upper triangular form

    if (norm == 0)
    {
        return;
    }

    for (n = nn-1; n >= 0; n--)
    {
        p = d_[n];
        q = e_[n];

        // scalar vector

        if (q == 0)
        {
            label l = n;
            H_(n, n) = 1;

            for (label i = n-1; i >= 0; i--)
            {
                w = H_(i, i) - p;
                r = 0;
                for (label j = l; j <= n; j++)
                {
                    r = r + H_(i, j)*H_(j, n);
                }

                if (e_[i] < 0)
                {
                    z = w;
                    s = r;
                }
                else
                {
                    l = i;

                    if (e_[i] == 0)
                    {
                        if (w != 0)
                        {
                            H_(i, n) = -r/w;
                        }
                        else
                        {
                            H_(i, n) = -r/(small*norm);
                        }

                        // Solve real equations
                    }
                    else
                    {
                        x = H_(i, i+1);
                        y = H_(i+1, i);
                        q = (d_[i] - p)*(d_[i] - p) + e_[i]*e_[i];
                        t = (x*s - z*r)/q;
                        H_(i, n) = t;

                        if (mag(x) > mag(z))
                        {
                            H_(i+1, n) = (-r - w*t)/x;
                        }
                        else
                        {
                            H_(i+1, n) = (-s - y*t)/z;
                        }
                    }

                    // Overflow control
                    t = mag(H_(i, n));
                    if ((small*t)*t > 1)
                    {
                        for (label j = i; j <= n; j++)
                        {
                            H_(j, n) = H_(j, n)/t;
                        }
                    }
                }
            }

            // Complex vector
        }
        else if (q < 0)
        {
            label l = n-1;
            scalar cdivr, cdivi;

            // Last vector component imaginary so matrix is triangular
            if (mag(H_(n, n-1)) > mag(H_(n-1, n)))
            {
                H_(n-1, n-1) = q/H_(n, n-1);
                H_(n-1, n) = -(H_(n, n) - p)/H_(n, n-1);
            }
            else
            {
                cdiv(cdivr, cdivi, 0, -H_(n-1, n), H_(n-1, n-1)-p, q);
                H_(n-1, n-1) = cdivr;
                H_(n-1, n) = cdivi;
            }

            H_(n, n-1) = 0;
            H_(n, n) = 1;

            for (label i = n-2; i >= 0; i--)
            {
                scalar ra = 0;
                scalar sa = 0;
                scalar vr, vi;

                for (label j = l; j <= n; j++)
                {
                    ra = ra + H_(i, j)*H_(j, n-1);
                    sa = sa + H_(i, j)*H_(j, n);
                }
                w = H_(i, i) - p;

                if (e_[i] < 0)
                {
                    z = w;
                    r = ra;
                    s = sa;
                }
                else
                {
                    l = i;
                    if (e_[i] == 0)
                    {
                        cdiv(cdivr, cdivi, -ra, -sa, w, q);
                        H_(i, n-1) = cdivr;
                        H_(i, n) = cdivi;
                    }
                    else
                    {
                        // Solve complex equations

                        x = H_(i, i+1);
                        y = H_(i+1, i);
                        vr = (d_[i] - p)*(d_[i] - p) + e_[i]*e_[i] - q*q;
                        vi = (d_[i] - p)*2*q;

                        if ((vr == 0) && (vi == 0))
                        {
                            vr = small*norm*(mag(w) + mag(q) +
                            mag(x) + mag(y) + mag(z));
                        }

                        cdiv
                        (
                            cdivr, cdivi,
                            x*r-z*ra+q*sa, x*s-z*sa-q*ra, vr, vi
                        );

                        H_(i, n-1) = cdivr;
                        H_(i, n) = cdivi;

                        if (mag(x) > (mag(z) + mag(q)))
                        {
                            H_(i+1, n-1) = (-ra - w*H_(i, n-1) + q*H_(i, n))/x;
                            H_(i+1, n) = (-sa - w*H_(i, n) - q*H_(i, n-1))/x;
                        }
                        else
                        {
                            cdiv
                            (
                                cdivr, cdivi,
                                -r-y*H_(i, n-1), -s-y*H_(i, n), z, q
                            );
                            H_(i+1, n-1) = cdivr;
                            H_(i+1, n) = cdivi;
                        }
                    }

                    // Overflow control

                    t = max(mag(H_(i, n-1)), mag(H_(i, n)));
                    if ((small*t)*t > 1)
                    {
                        for (label j = i; j <= n; j++)
                        {
                            H_(j, n-1) = H_(j, n-1)/t;
                            H_(j, n) = H_(j, n)/t;
                        }
                    }
                }
            }
        }
    }

    // Vectors of isolated roots
    for (label i = 0; i < nn; i++)
    {
        if (i < low || i > high)
        {
            for (label j = i; j < nn; j++)
            {
                V_(i, j) = H_(i, j);
            }
        }
    }

    // Back transformation to get eigenvectors of original matrix
    for (label j = nn-1; j >= low; j--)
    {
        for (label i = low; i <= high; i++)
        {
            z = 0;
            for (label k = low; k <= min(j, high); k++)
            {
                z = z + V_(i, k)*H_(k, j);
            }
            V_(i, j) = z;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eigendecomposition::eigendecomposition(const scalarSquareMatrix& A)
:
    symmetric_(false),
    d_(A.n()),
    e_(A.n()),
    V_(A.n())
{
    const label n = A.n();

    for (label j = 0; (j < n) && symmetric_; j++)
    {
        for (label i = 0; (i < n) && symmetric_; i++)
        {
            symmetric_ = (A(i, j) == A(j, i));
        }
    }

    if (symmetric_)
    {
        for (label i = 0; i < n; i++)
        {
            for (label j = 0; j < n; j++)
            {
                V_(i, j) = A(i, j);
            }
        }

        // Tridiagonalise
        tred2();

        // Diagonalise
        tql2();
    }
    else
    {
        H_.setSize(n);
        ort_.setSize(n);

        for (label j = 0; j < n; j++)
        {
            for (label i = 0; i < n; i++)
            {
                H_(i, j) = A(i, j);
            }
        }

        // Reduce to Hessenberg form
        orthes();

        // Reduce Hessenberg to real Schur form
        hqr2();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eigendecomposition::D(scalarSquareMatrix& D_) const
{
    const label n = V_.n();
    D_.setSize(n);
    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            D_(i, j) = 0;
        }
        D_(i, i) = d_[i];
        if (e_[i] > 0)
        {
            D_(i, i+1) = e_[i];
        }
        else if (e_[i] < 0)
        {
            D_(i, i-1) = e_[i];
        }
    }
}


// ************************************************************************* //
