### ADAPTIVE GAUSS-LOBATTO FROM GANDER AND GAUTSCHI (2000)

### constants from section 4 of the paper
const α = sqrt(2/3)
const β = 1/sqrt(5)

const X1 = 0.94288241569547971905635175843185720232
const X2 = 0.64185334234578130578123554132903188354
const X3 = 0.23638319966214988028222377349205292599

const A = 0.015827191973480183087169986733305510591
const B = 0.094273840218850045531282505077108171960
const C = 0.155071987336585396253963597980210298680
const D = 0.18882157396018245442000533937297167125
const E = 0.199773405226858526792068022066048840246
const F = 0.22492646533333952701601768799639508076
const G = 0.24261107190140773379964095790325635233




"""
    adaptlob(f,a,b; tol::Real=1e-10, trace::Bool=false)

Numerically integrate `f` over the closed real interval `[a,b]` to within tolerance `tol` using the adaptive Gauss-Lobatto method from Gander and Gautschi (2000).

Passing `trace=true` prints out information about where the function is being sampled.

`f` must be defined at `a` and `b`.

```
using QuadGG
adaptlob(sqrt,0,1)
```
"""
function adaptlob(f,a,b; tol::Real=1e-10, trace::Bool=false)
    m = (a+b)/2
    h = (b-a)/2

    x = (a,m-X1*h,m-α*h,m-X2*h,m-β*h,m-X3*h,m,m+X3*h,m+β*h,m+X2*h,m+α*h,m+X1*h,b)
    y = f.(x)

    
    fa,fb = y[1],y[end]

    _i2 = (h/6) * (y[1] + y[13] + 5*(y[5]+y[9]))
    _i1 = (h/1470) * (77*(y[1]+y[13]) + 432*(y[3]+y[11]) + 625 * (y[2]+y[12]) + 672*y[7])

    _is = h * (A*(y[1]+y[13]) + B*(y[2]+y[12]) + C*(y[3]+y[11]) + D*(y[4]+y[10]) + E*(y[5]+y[9]) + F*(y[6]+y[8]) + G*y[7])

    sgn = sign(_is)
    sgn = ifelse(iszero(sgn),1,sgn)
    erri1 = abs(_i1 - _is)
    erri2 = abs(_i2 - _is)
    
    R = ifelse(!iszero(erri2), erri1/erri2, 1)

    if 0 < R < 1
        tol = tol/R
    end

    _is_new = sgn * abs(_is) * tol / eps()

    if iszero(_is_new)
        _is_new = b-a
    end

    flag = WarnFlag(false)
    Q = adaptlobstp(f,a,b,fa,fb,_is_new,trace,flag)
    return Q
end


function adaptlobstp(f,a,b,fa,fb,_is,trace,flag)
    m = (a+b)/2
    h = (b-a)/2

    mll = m-α*h
    ml = m-β*h
    mr = m+β*h
    mrr = m+α*h

    fmll = f(mll)
    fml = f(ml)
    fm = f(m)
    fmr = f(mr)
    fmrr = f(mrr)

    i2=(h/6)*(fa+fb+5*(fml+fmr))
    i1=(h/1470)*(77*(fa+fb)+432*(fmll+fmrr)+625*(fml+fmr)+672*fm)

    if ≈(_is+i1-i2, _is, atol=1e-16) || mll <= a || b <= mrr || !isfinite(_is)
        if (m <= a || b <= m) && !flag.nowarn
            @warn "Interval contains no more machine numbers. Required tolerance may not be met."
            flag.nowarn = true
        end
        Q = i1
        trace && @printf "%20.15f    %20.15f    %20.15f\n" a b-a Q
    else
        Q = adaptlobstp(f,a,mll,fa,fmll,_is,trace,flag) +
            adaptlobstp(f,mll,ml,fmll,fml,_is,trace,flag) +
            adaptlobstp(f,ml,m,fml,fm,_is,trace,flag) +
            adaptlobstp(f,m,mr,fm,fmr,_is,trace,flag) +
            adaptlobstp(f,mr,mrr,fmr,fmrr,_is,trace,flag) +
            adaptlobstp(f,mrr,b,fmrr,fb,_is,trace,flag)
    end
    return Q
end