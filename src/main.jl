begin
    using DifferentialEquations
    using Sundials
    using GLMakie
    using LinearAlgebra
    # using DelimitedFiles
    using ForwardDiff
    using DataFrames
    using PrettyTables
    using BenchmarkTools, Profile
    using MAT
    include("hybrid_cable_model.jl")
    import .HybridCableModel
end
begin
    include("hybrid_cable_model.jl")
    include("util.jl")
    include("ode.jl")
    include("bc.jl")
    include("post.jl")
end

begin
    N = 3
    tspan = (0.0, 6.0)
    tol = 1e-8

    L = 0.42 + 1.04 * 1
    pulley_num = 1

    rd = 1.5e-2
    Id = 1.7485e-5 + 1.11e-3
    # REAL PARAMETERS
    Td = 0.155 * rd
    Cd = 0.020 * rd
    # Td = 0.0 * rd
    # Cd = 0.0 * rd

    rp = 1.5e-2
    Ip = 1.7485e-5
    Tp = Td
    Cp = Cd

    T = 48.37 * rd

    k = 1290.44
    m = 0.65

    mt = FullModel()
    # mt = PartialModel()
    # mt = SimpleModel()

    if isa(mt, SimpleModel)
        Cd = Cd + pulley_num * Cp
        Td = Td + pulley_num * Tp
        Id = Id + pulley_num * Ip
    end

    mp = ModelParam(
        N=N,
        tspan=tspan,
        tol=tol,
        L=L,
        rd=rd,
        Id=Id,
        k=k,
        T=T,
        Td=Td,
        Cd=Cd,
        pulley_num=pulley_num,
        rp=fill(rp, pulley_num),
        lp=[0.42 + 1.04 * (i - 1) for i = 1:pulley_num],
        Ip=fill(Ip, pulley_num),
        Tp=fill(Tp, pulley_num),
        Cp=fill(Cp, pulley_num)
    )
    display(mt)
    show_model_param(mp)
end

begin
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    # @time sol = solve(prob, CVODE_Adams(), reltol=mp.tol, abstol=mp.tol)  # non stiff
    # @time sol = solve(prob, CVODE_BDF(), reltol=mp.tol, abstol=mp.tol)  # stiff
    display("solver finished!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 40 seconds = 80
end

post_sol(mt, mp, sol)

# generate dispalcement video for all time steps
# begin
#     frame_rate = 24
#     sol_t_range = LinRange(sol.t[1], sol.t[end], frame_rate * 5)
#     length_size = 1000

#     function get_all_deformation(sol_t, sol, mp, mt, length_size=1000)
#         sol_u_all = zeros(length_size)
#         x = sol(sol_t)[1:mp.N+mp.pulley_num+1]
#         ux = LinRange(0, x[end], length_size)
#         for i in 1:length_size
#             sol_u_all[i] = u(mt, mp, x, ux[i]) * 1000 # unit: mm
#         end
#         return sol_u_all
#     end

#     function get_all_strains(sol_t, sol, mp, mt, length_size=1000)
#         sol_uₓ_all = zeros(length_size)
#         x = sol(sol_t)[1:mp.N+mp.pulley_num+1]
#         ux = LinRange(0, x[end], length_size)
#         for i in 1:length_size
#             sol_uₓ_all[i] = ∂ₓu(mt, mp, x, ux[i]) * 1000 # unit: mm
#         end
#         return sol_uₓ_all
#     end

#     fig = Figure()
#     ax1 = Axis(fig[1, 1:2], xlabel="coordinate(x/xe)", ylabel="deformation(mm)", title="t=" * string(0))
#     ax2 = Axis(fig[2, 1:2], xlabel="coordinate(x/xe)", ylabel="strain")
#     ux = LinRange(0, 1, length_size)
#     # record(fig, "ideal_1pulley_deformation_with_fixlimit_strain.mp4", sol_t_range; framerate=frame_rate) do t
#     #     display(t)
#     #     sol_u_all = get_all_deformation(t, sol, mp, mt)
#     #     sol_uₓ_all = get_all_strains(t, sol, mp, mt)
#     #     empty!(ax1)
#     #     empty!(ax2)

#     #     lines!(ax1, ux, sol_u_all, label="u_all", color=:tomato)
#     #     lines!(ax2, ux, sol_uₓ_all, label="u_all", color=:tomato)

#     #     ax1.title = "t=" * string(t)

#     #     ylims!(ax1, -4e-1, 1e-5)

#     #     display(fig)
#     # end

#     record(fig, "ideal_1pulley_deformation_with_autolimit_strain.mp4", sol_t_range; framerate=frame_rate) do t
#         display(t)
#         sol_u_all = get_all_deformation(t, sol, mp, mt)
#         sol_uₓ_all = get_all_strains(t, sol, mp, mt)
#         empty!(ax1)
#         empty!(ax2)

#         lines!(ax1, ux, sol_u_all, label="u_all", color=:tomato)
#         lines!(ax2, ux, sol_uₓ_all, label="u_all", color=:tomato)
#         autolimits!(ax1)
#         autolimits!(ax2)
#         ax1.title = "t=" * string(t)

#         display(fig)
#     end
# end

# FullModel or PartialModel 3D deformation data
begin
    length_size = 300

    sol_t_range = LinRange(sol.t[1], sol.t[end], length_size)
    sol_x_range = LinRange(0, 1, length_size)
    # x axis is the coordinate of the cable, y axis is the time, z axis is the deformation/strain
    sol_u_all = zeros(length_size, length(sol_t_range))
    sol_uₓ_all = zeros(length_size, length(sol_t_range))

    function get_all_deformation(sol_t, sol, mp, mt, length_size=1000)
        sol_u_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+mp.pulley_num+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_u_all[i] = u(mt, mp, x, ux[i]) * 1000 # unit: mm
        end
        return sol_u_all
    end

    function get_all_strains(sol_t, sol, mp, mt, length_size=1000)
        sol_uₓ_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+mp.pulley_num+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_uₓ_all[i] = ∂ₓu(mt, mp, x, ux[i])
        end
        return sol_uₓ_all
    end

    for id in eachindex(sol_t_range)
        t = sol_t_range[id]
        sol_u_all[:, id] = get_all_deformation(t, sol, mp, mt, length_size)
        sol_uₓ_all[:, id] = get_all_strains(t, sol, mp, mt, length_size)
    end

    GLMakie.activate!()
    fig = surface(sol_x_range, sol_t_range, sol_u_all,
        axis=(type=Axis3, azimuth=pi / 4,
            zlabel="deformation", xlabel="x̃", ylabel="t"))

    display(GLMakie.Screen(), fig)

    GLMakie.activate!()
    fig = surface(sol_x_range, sol_t_range, sol_uₓ_all,
        axis=(type=Axis3, azimuth=pi / 4,
            zlabel="strain", xlabel="x̃", ylabel="t"))
    display(GLMakie.Screen(), fig)
    display("3D deformation data generated!")
end

# SimpleModel 3D deformation data
begin
    length_size = 300

    sol_t_range = LinRange(sol.t[1], sol.t[end], length_size)
    sol_x_range = LinRange(0, 1, length_size)
    # x axis is the coordinate of the cable, y axis is the time, z axis is the deformation/strain
    sol_u_all = zeros(length_size, length(sol_t_range))
    sol_uₓ_all = zeros(length_size, length(sol_t_range))

    function get_all_deformation(sol_t, sol, mp, mt, length_size=1000)
        sol_u_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_u_all[i] = u(mt, mp, x, ux[i]) * 1000 # unit: mm
        end
        return sol_u_all
    end

    function get_all_strains(sol_t, sol, mp, mt, length_size=1000)
        sol_uₓ_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_uₓ_all[i] = ∂ₓu(mt, mp, x, ux[i])
        end
        return sol_uₓ_all
    end

    for id in eachindex(sol_t_range)
        t = sol_t_range[id]
        sol_u_all[:, id] = get_all_deformation(t, sol, mp, mt, length_size)
        sol_uₓ_all[:, id] = get_all_strains(t, sol, mp, mt, length_size)
    end

    GLMakie.activate!()
    fig = surface(sol_x_range, sol_t_range, sol_u_all,
        axis=(type=Axis3, azimuth=pi / 4,
            zlabel="deformation", xlabel="x̃", ylabel="t"))

    display(GLMakie.Screen(), fig)

    GLMakie.activate!()
    fig = surface(sol_x_range, sol_t_range, sol_uₓ_all,
        axis=(type=Axis3, azimuth=pi / 4,
            zlabel="strain", xlabel="x̃", ylabel="t"))
    display(GLMakie.Screen(), fig)
    display("3D deformation data generated!")
end

# plot the 3D deformation and strain field
# surface(sol_x_range, sol_t_range, sol_u_all)
# surface(sol_x_range, sol_t_range, sol_uₓ_all)

# export 3D deformation data
# begin
# include("data/3D strain/ideal_0pulley_deformation_N3_3D.jl")
#     # include("data/3D strain/ideal_0pulley_deformation_N10_3D.jl")
#     # include("data/3D strain/ideal_1pulley_deformation_N3_3D.jl")
#     # include("data/3D strain/ideal_1pulley_deformation_N10_3D.jl")
# end

# exprt 2D deformation data
# begin
#     # export post data
#     # include("data/2D strain/ideal_0pulley_deformation_N3.jl")
#     # include("data/2D strain/ideal_0pulley_deformation_N10.jl")
#     # include("data/2D strain/ideal_1pulley_deformation_N3.jl")
#     # include("data/2D strain/ideal_1pulley_deformation_N10.jl")
#     # display("2d deformation and strain data exported!")
# end
