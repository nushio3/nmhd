set xrange [0.01:1]
set xlabel 'B_0'
set ylabel 'omega'

kva(b) = 2 * pi * b / (5.0/3.0)**2
kva2(b) = kva(b)**2
Omega = 1.0
Omega2 = Omega**2

omega2(b) = 0.5*(Omega2 + 2*kva2(b) - Omega * sqrt(Omega2 + 16*kva2(b)))

omega_stable(b) = sqrt(omega2(b))
omega_unstable(b) = sqrt(-omega2(b))

set samples 1e4
plot omega_stable(x), omega_unstable(x) lt 3
