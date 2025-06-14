using .HybridCableModel

function ode_function!(dX, X, mp, t)
    ode_delegate!(mt, dX, X, mp, t)
end

function ode_delegate!(mt::FullModel, dX, X, mp, t)
    x = X[1:(mp.N+mp.pulley_num+1)]
    dx = X[(mp.N+mp.pulley_num+2):end]

    _s = s(mt, mp, x)
    _ds = ds(mt, mp, x, dx)
    _M = M(mt, mp, _s, _ds)
    _F = F(mt, mp, _s, _ds)
    _G = G(mt, mp, x, dx)
    _C = C(mt, mp, x, dx)

    id_M = I(mp.N + mp.pulley_num + 1)
    zero_c = zeros(mp.N + mp.pulley_num + 1, 1)

    _A = _M * [id_M; _G]
    _B = _F - _M * [zero_c; _C]
    # write(io, "current time: " * string(t) * "\n")
    # write(io, "condition number of _A: " * string(cond(_A)) * "\n")
    # write(io, "condition number of _M: " * string(cond(_M)) * "\n")
    # write(io, "condition number of [id_M; _G]: " * string(cond([id_M; _G])) * "\n")
    # write(io, "diff: cond(A)-cond(M)*cond(G) = " * string(cond(_A) - cond(_M) * cond([id_M; _G])) * "\n\n")
    ddx = inv(_A) * _B
    dX[1:mp.N+mp.pulley_num+1] = dx
    dX[mp.N+mp.pulley_num+2:end] = ddx
    nothing
end

function ode_delegate!(mt::PartialModel, dX, X, mp, t)
    x = X[1:(mp.N+mp.pulley_num+1)]
    dx = X[(mp.N+mp.pulley_num+2):end]

    _s = s(mt, mp, x)
    _ds = ds(mt, mp, x, dx)
    _M = M(mt, mp, _s, _ds)
    _F = F(mt, mp, _s, _ds)
    _G = G(mt, mp, x, dx)
    _C = C(mt, mp, x, dx)

    id_M = I(mp.N + mp.pulley_num + 1)
    zero_c = zeros(mp.N + mp.pulley_num + 1, 1)

    _A = _M * [id_M; _G]
    _B = _F - _M * [zero_c; _C]
    # write(io, "current time: " * string(t) * "\n")
    # write(io, "condition number of _A: " * string(cond(_A)) * "\n")
    # write(io, "condition number of _M: " * string(cond(_M)) * "\n")
    # write(io, "condition number of [id_M; _G]: " * string(cond([id_M; _G])) * "\n")
    # write(io, "diff: cond(A)-cond(M)*cond(G) = " * string(cond(_A) - cond(_M) * cond([id_M; _G])) * "\n\n")
    ddx = inv(_A) * _B
    dX[1:mp.N+mp.pulley_num+1] = dx
    dX[mp.N+mp.pulley_num+2:end] = ddx
    nothing
end

function ode_delegate!(mt::SimpleModel, dX, X, mp, t)
    x = X[1:(mp.N+1)]
    dx = X[(mp.N+2):end]

    _M = M(mt, mp, x, dx)
    _F = F(mt, mp, x, dx)
    _K = K(mt, mp, x, dx)

    ddx = inv(_M) * (_F - _K * x)
    dX[1:mp.N+1] = dx
    dX[mp.N+2:end] = ddx
    nothing
end

# =============== FULL MODEL ====================

function M(mt::FullModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+mp.pulley_num+1]
    xe = x[end]
    θ = s[mp.N+mp.pulley_num+2]
    θₚ = s[mp.N+mp.pulley_num+3:end]

    dx = ds[1:mp.N+mp.pulley_num+1]
    dxe = dx[end]
    dθ = ds[mp.N+mp.pulley_num+2]
    dθₚ = ds[mp.N+mp.pulley_num+3:end]

    M11 = (-mp.ρ * xe) .* [∫Φ(mp, xe, i) for i in 1:(mp.N+mp.pulley_num)]'
    M12 = mp.m + mp.ρ * xe
    for i in 1:(mp.N+mp.pulley_num)
        M12 += mp.ρ * x[i] * ∫x∂ₓΦ(mp, xe, i)
    end
    M13 = -(mp.Id + mp.ρ * mp.rd^3 * θ) / mp.rd
    M14 = zeros(1, mp.pulley_num)
    for i in 1:mp.pulley_num
        M14[i] = mp.Ip[i] / mp.rp[i] * (∂ₓₑuₚ(mt, mp, i, x) - 1)
    end

    M21 = zeros(mp.N + mp.pulley_num, mp.N + mp.pulley_num)
    for i in 1:mp.N+mp.pulley_num
        for j in 1:mp.N+mp.pulley_num
            M21[i, j] = ∫ΦΦ(mp, xe, i, j)
        end
    end

    M22 = zeros(mp.N + mp.pulley_num, 1)
    M22_1 = zeros(mp.N + mp.pulley_num, 1)
    M22_2 = zeros(mp.N + mp.pulley_num, 1)
    for i in 1:mp.N+mp.pulley_num
        M22_1[i] = 0
        for j in 1:mp.N+mp.pulley_num
            M22_1[i] += ∫x∂ₓΦΦ(mp, xe, j, i) * x[j]
        end
    end
    M22_1 = -(1 / xe) * M22_1

    for i in 1:mp.N+mp.pulley_num
        M22_2[i] = -∫Φ(mp, xe, i)
    end

    M22 = M22_1 + M22_2

    M23 = zeros(mp.N + mp.pulley_num, 1)

    M24 = zeros(mp.N + mp.pulley_num, mp.pulley_num)
    for i in 1:mp.N+mp.pulley_num
        for j in 1:mp.pulley_num
            M24[i, j] = mp.Ip[j] / (mp.ρ * mp.rp[j]) * Φ(mp, xe, i, 1 - mp.lp[j] / xe)
        end
    end

    _M = [M11 M12 M13 M14;
        M21 M22 M23 M24
    ]
    return _M
end

function F(mt::FullModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+mp.pulley_num+1]
    xe = x[end]
    θ = s[mp.N+mp.pulley_num+2]
    θₚ = s[mp.N+mp.pulley_num+3:end]

    dx = ds[1:mp.N+mp.pulley_num+1]
    dxe = dx[end]
    dθ = ds[mp.N+mp.pulley_num+2]
    dθₚ = ds[mp.N+mp.pulley_num+3:end]
    a = sqrt(mp.A * mp.E / mp.ρ)

    Fnon = mp.ρ * dxe^2 + mp.k * (xe - mp.L) - fm(dxe, mp.Tm, mp.Cm)
    Fnon += -1 / mp.rd * (1 / 2 * mp.ρ * mp.rd^3 * dθ^2 + f(dθ, mp.Td, mp.Cd) - mp.T)
    for i in 1:mp.pulley_num
        Fnon += (∂ₓₑuₚ(mt, mp, i, x) - 1) / mp.rp[i] * f(dθₚ[i], mp.Tp[i], mp.Cp[i])
    end

    Fs1 = zeros(mp.N + mp.pulley_num + 1, 1)
    Fs1[1] = -Fnon + mp.ρ * xe * ∫ζ(mt, mp, x, dx)
    for i in 1:mp.N+mp.pulley_num
        _tmp = 0
        for j in 1:mp.N+mp.pulley_num
            _tmp += ∫∂ₓΦ∂ₓΦ(mp, xe, j, i) * x[j]
        end
        Fs1[i+1] = -a^2 / xe^2 * _tmp
    end

    Fs2 = zeros(mp.N + mp.pulley_num + 1, 1)
    for i in 1:mp.N+mp.pulley_num
        _tmp = 0
        for j in 1:mp.pulley_num
            _tmp += f(dθₚ[j], mp.Tp[j], mp.Cp[j]) / (mp.ρ * mp.rp[j]) * Φ(mp, xe, i, 1 - mp.lp[j] / xe)
        end
        Fs2[i+1] = -(_tmp + ∫ζΦ(mt, mp, i, x, dx))
    end

    _F = Fs1 + Fs2
    return _F
end

function G(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _G = zeros(mp.pulley_num + 1, mp.N + mp.pulley_num + 1)
    for i in 1:mp.N+mp.pulley_num
        _G[1, i] = 0
    end
    _G[1, mp.N+mp.pulley_num+1] = -1 / mp.rd

    for i in 1:mp.pulley_num
        for j in 1:mp.N+mp.pulley_num
            _G[i+1, j] = Φ(mp, xe, j, 1 - mp.lp[i] / xe) / mp.rp[i]
        end
        _G[i+1, mp.N+mp.pulley_num+1] = (∂ₓu(mt, mp, x, xe - mp.lp[i]) - 1) / mp.rp[i]
        _tmp = 0
        for j in 1:mp.N+mp.pulley_num
            _tmp += ∂ₓΦ(mp, xe, j, 1 - mp.lp[i] / xe) * x[j]
        end
        _G[i+1, mp.N+mp.pulley_num+1] += -(1 - mp.lp[i] / xe) / (mp.rp[i] * xe) * _tmp
    end

    return _G
end

function C(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _C = zeros(mp.pulley_num + 1, 1)
    _C[1] = 0
    for i in 1:mp.pulley_num
        _C[i+1] = 2 / mp.rp[i] * ∂ₓₜu(mt, mp, x, dx, xe - mp.lp[i]) * dxe
        _C[i+1] += 1 / mp.rp[i] * ∂ₓₓu(mt, mp, x, xe - mp.lp[i]) * dxe^2
        _C[i+1] += 1 / mp.rp[i] * ζ(mt, mp, x, dx, 1 - mp.lp[i] / xe)
    end

    return _C
end


# =============== PARTIAL MODEL ====================

function M(mt::PartialModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+mp.pulley_num+1]
    xe = x[end]
    θ = s[mp.N+mp.pulley_num+2]
    θₚ = s[mp.N+mp.pulley_num+3:end]

    dx = ds[1:mp.N+mp.pulley_num+1]
    dxe = dx[end]
    dθ = ds[mp.N+mp.pulley_num+2]
    dθₚ = ds[mp.N+mp.pulley_num+3:end]

    M11 = (-mp.ρ * xe) .* [∫Φ(mp, xe, i) for i in 1:(mp.N+mp.pulley_num)]'
    M12 = mp.m + mp.ρ * xe
    M13 = -(mp.Id + mp.ρ * mp.rd^3 * θ) / mp.rd
    M14 = zeros(1, mp.pulley_num)
    for i in 1:mp.pulley_num
        M14[i] = mp.Ip[i] / mp.rp[i] * (∂ₓₑuₚ(mt, mp, i, x) - 1)
    end

    M21 = zeros(mp.N + mp.pulley_num, mp.N + mp.pulley_num)
    for i in 1:mp.N+mp.pulley_num
        for j in 1:mp.N+mp.pulley_num
            M21[i, j] = ∫ΦΦ(mp, xe, i, j)
        end
    end

    M22 = zeros(mp.N + mp.pulley_num, 1)
    M22_2 = zeros(mp.N + mp.pulley_num, 1)

    for i in 1:mp.N+mp.pulley_num
        M22_2[i] = -∫Φ(mp, xe, i)
    end

    M22 = M22_2

    M23 = zeros(mp.N + mp.pulley_num, 1)

    M24 = zeros(mp.N + mp.pulley_num, mp.pulley_num)
    for i in 1:mp.N+mp.pulley_num
        for j in 1:mp.pulley_num
            M24[i, j] = mp.Ip[j] / (mp.ρ * mp.rp[j]) * Φ(mp, xe, i, 1 - mp.lp[j] / xe)
        end
    end

    _M = [M11 M12 M13 M14;
        M21 M22 M23 M24
    ]
    return _M
end

function F(mt::PartialModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+mp.pulley_num+1]
    xe = x[end]
    θ = s[mp.N+mp.pulley_num+2]
    θₚ = s[mp.N+mp.pulley_num+3:end]

    dx = ds[1:mp.N+mp.pulley_num+1]
    dxe = dx[end]
    dθ = ds[mp.N+mp.pulley_num+2]
    dθₚ = ds[mp.N+mp.pulley_num+3:end]
    a = sqrt(mp.A * mp.E / mp.ρ)

    Fnon = mp.ρ * dxe^2 + mp.k * (xe - mp.L) - fm(dxe, mp.Tm, mp.Cm)
    Fnon += -1 / mp.rd * (1 / 2 * mp.ρ * mp.rd^3 * dθ^2 + f(dθ, mp.Td, mp.Cd) - mp.T)
    for i in 1:mp.pulley_num
        Fnon += (∂ₓₑuₚ(mt, mp, i, x) - 1) / mp.rp[i] * f(dθₚ[i], mp.Tp[i], mp.Cp[i])
    end

    Fs1 = zeros(mp.N + mp.pulley_num + 1, 1)
    Fs1[1] = -Fnon
    for i in 1:mp.N+mp.pulley_num
        _tmp = 0
        for j in 1:mp.N+mp.pulley_num
            _tmp += ∫∂ₓΦ∂ₓΦ(mp, xe, j, i) * x[j]
        end
        Fs1[i+1] = -a^2 / xe^2 * _tmp
    end

    Fs2 = zeros(mp.N + mp.pulley_num + 1, 1)
    for i in 1:mp.N+mp.pulley_num
        _tmp = 0
        for j in 1:mp.pulley_num
            _tmp += f(dθₚ[j], mp.Tp[j], mp.Cp[j]) / (mp.ρ * mp.rp[j]) * Φ(mp, xe, i, 1 - mp.lp[j] / xe)
        end
        Fs2[i+1] = -_tmp
    end

    _F = Fs1 + Fs2
    return _F
end

function G(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _G = zeros(mp.pulley_num + 1, mp.N + mp.pulley_num + 1)
    for i in 1:mp.N+mp.pulley_num
        _G[1, i] = 0
    end
    _G[1, mp.N+mp.pulley_num+1] = -1 / mp.rd

    for i in 1:mp.pulley_num
        for j in 1:mp.N+mp.pulley_num
            _G[i+1, j] = Φ(mp, xe, j, 1 - mp.lp[i] / xe) / mp.rp[i]
        end
        _G[i+1, mp.N+mp.pulley_num+1] = (∂ₓu(mt, mp, x, xe - mp.lp[i]) - 1) / mp.rp[i]
    end

    return _G
end

function C(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _C = zeros(mp.pulley_num + 1, 1)
    _C[1] = 0
    for i in 1:mp.pulley_num
        _C[i+1] = 2 / mp.rp[i] * ∂ₓₜu(mt, mp, x, dx, xe - mp.lp[i]) * dxe
        _C[i+1] += 1 / mp.rp[i] * ∂ₓₓu(mt, mp, x, xe - mp.lp[i]) * dxe^2
    end

    return _C
end

# =============== SIMPLE MODEL ====================

function M(mt::SimpleModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    M11 = zeros(mp.N, mp.N)
    for i in 1:mp.N
        for j in 1:mp.N
            M11[i, j] = ∫ΦΦ(mp, xe, i, j)
        end
    end
    M12 = -[∫Φ(mp, xe, i) for i in 1:mp.N]
    M21 = M12'
    M22 = 1 / xe * (mp.m / mp.ρ + mp.Id / (mp.ρ * mp.rd^2) + mp.L)
    _M = [M11 M12; M21 M22]
    return _M
end


function K(mt::SimpleModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]
    a = sqrt(mp.A * mp.E / mp.ρ)

    K11 = zeros(mp.N, mp.N)
    for i in 1:mp.N
        for j in 1:mp.N
            K11[i, j] = ∫∂ₓΦ∂ₓΦ(mp, xe, j, i)
        end
    end
    K11 = a^2 / xe^2 * K11

    K12 = zeros(mp.N, 1)
    K21 = K12'
    K22 = 0.0
    _K = [K11 K12; K21 K22]
    return _K
end

function F(mt::SimpleModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]
    _dθ = dθ(mt, mp, dxe)
    F1 = zeros(mp.N, 1)
    F2 = -1 / 2 * mp.ρ * dxe^2 - (mp.T - f(_dθ, mp.Td, mp.Cd)) / mp.rd
    F2 = F2 + mp.k * (mp.L - xe) + fm(dxe, mp.Tm, mp.Cm)
    F2 = F2 / (mp.ρ * xe)
    _F = [F1; F2]
    return _F
end