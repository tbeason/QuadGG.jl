### ADAPTIVE SIMPSON RULE FROM GANDER AND GAUTSCHI (2000)


"""
    adaptsim(f,a,b; tol::Real=1e-10, trace::Bool=false)

Numerically integrate `f` over the closed real interval `[a,b]` to within tolerance `tol` using the adaptive Simpson's rule from Gander and Gautschi (2000).

Passing `trace=true` prints out information about where the function is being sampled.

`f` must be defined at `a` and `b`.

```
using QuadGG
adaptsim(sqrt,0,1)
```
"""
function adaptsim(f,a,b; tol::Real=1e-10, trace::Bool=false)
    m = (a+b)/2

    fa = f(a)
    fm = f(m)
    fb = f(b)

    sumyy = sum(f(v) for v in a .+ [0.9501,0.2311,0.6068,0.4860,0.8913].*(b-a))

    _is = (b-a) / 8 * (fa + fm + fb + sumyy)

    if iszero(_is)
        _is = b-a
    end

    _is_new = _is * tol / eps()

    flag = WarnFlag(false)
    Q = adaptsimstp(f,a,b,fa,fm,fb,_is_new,trace,flag)
    return Q
end


function adaptsimstp(f,a,b,fa,fm,fb,_is,trace,flag)
    m = (a+b)/2
    h = (b-a)/4

    fml = f(a+h)
    fmr = f(b-h)
    
    i1 = h/1.5 * (fa + 4*fm + fb)
    i2 = h/3 * (fa + 4*(fml + fmr) + 2*fm + fb)
    i1 = (16*i2 - i1)/15

    if ≈(_is+i1-i2, _is, atol=1e-16) || m <= a || b <= m|| !isfinite(_is)
        if (m <= a || b <= m) && !flag.nowarn
            @warn "Interval contains no more machine numbers. Required tolerance may not be met."
            flag.nowarn = true
        end
        Q = i1
        trace && @printf "%20.15f    %20.15f    %20.15f\n" a b-a Q
    else
        Q = adaptsimstp(f,a,m,fa,fml,fm,_is,trace,flag) + adaptsimstp(f,m,b,fm,fmr,fb,_is,trace,flag)
    end
    return Q
end