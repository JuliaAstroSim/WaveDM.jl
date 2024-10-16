RAR(a_b, a0) = a_b / (1 - exp(-sqrt(abs(a_b)/a0)))

Milgrom(::Missing, ::Any, ::Any) = missing
Milgrom(::Any, ::Missing, ::Any) = missing
Milgrom(a_b, a0, nuIndex) = nu(abs(a_b)/a0, nuIndex) * a_b

RAR_inv(a_b, a_MOND) = a_MOND > a_b ? a_b / (log(1-a_b/a_MOND))^2 : missing
Milgrom_inv(a_b, a_MOND, nuIndex) = a_MOND > a_b ? a_b * (((2*(a_MOND/a_b)^nuIndex - 1)^2 - 1) / 4)^(1/nuIndex) : missing

# Lsq fitting models
# modelRAR(x, p) = x ./ (1 .- exp.(-sqrt.(abs.(x)/abs(p[1]))))
modelRAR(x, p) = x ./ (1 .- exp.(-sqrt.(abs.(x)/p[1])))