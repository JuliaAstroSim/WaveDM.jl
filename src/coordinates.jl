function vec_cartesian_to_spherical(Ax, Ay, Az, xxx, yyy, zzz, rrr)
    sqrt_xy = sqrt.(xxx.^2 + yyy.^2)
    Ar = (Ax .* xxx + Ay .* yyy + Az .* zzz) ./ rrr
    Aθ = (Ax .* xxx .* zzz + Ay .* yyy .* zzz - Az .* (xxx.^2 + yyy.^2)) ./ (rrr .* sqrt_xy)
    Aϕ = (-Ax .* yyy + Ay .* xxx) ./ sqrt_xy
    return Ar, Aθ, Aϕ
end

function vec_cartesian_to_cylindrical(Ax, Ay, Az, xxx, yyy, RRR)
    A_R = (Ax .* xxx + Ay .* yyy) ./ RRR
    A_phi = (-Ax .* yyy + Ay .* xxx) ./ RRR
    A_z = Az
    return A_R, A_phi, A_z
end
