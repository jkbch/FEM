using LinearAlgebra
using SparseArrays
using Plots
using Polynomials
using SymRCM
using AMD

function xy(
    x0::Float64, 
    y0::Float64, 
    L1::Float64, 
    L2::Float64, 
    noelms1::Int64, 
    noelms2::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}
    VX = repeat(collect(LinRange(x0, x0+L1, noelms1+1)), inner=noelms2+1)
    VY = repeat(collect(LinRange(y0+L2, y0, noelms2+1)), noelms1+1)
    return VX, VY
end

function conelmtab(noelms1::Int64, noelms2::Int64)::Matrix{Int64}
    k = [i for i in 1:(noelms1*(noelms2+1)) if i % (noelms2+1) != 0]

    return [
        k (2 + noelms2 .+ k) (1 + noelms2 .+ k);
        k (1 .+ k) (2 + noelms2 .+ k)
    ]
end

function basfun(VX, VY, EToV)
    xjs = VX[EToV[:, [2,3,1]]]
    yjs = VY[EToV[:, [2,3,1]]]

    xks = VX[EToV[:, [3,1,2]]]
    yks = VY[EToV[:, [3,1,2]]]

    as = xjs .* yks - xks .* yjs
    bs = yjs - yks
    cs = xks - xjs

    return as, bs, cs
end

function constructBeds(
    VX::Vector{Float64},
    VY::Vector{Float64},
    EToV::Matrix{Int64},
    tol::Float64,
    fd::Function,
)::Matrix{Int64}
    xc = (VX[EToV] + VX[EToV[:, [2, 3, 1]]]) ./ 2
    yc = (VY[EToV] + VY[EToV[:, [2, 3, 1]]]) ./ 2
    return getindex.(findall(abs.(fd.(xc, yc)) .<= tol), [1 2])
end

function constructBnodes(
    VX::Vector{Float64},
    VY::Vector{Float64},
    tol::Float64,
    fd::Function,
)::Vector{Int64}
    return findall(abs.(fd.(VX, VY)) .<= tol)
end

function dirbc(
    bnodes::Vector{Int64}, 
    f::Vector{Float64}, 
    A::SparseMatrixCSC{Float64, Int64}, 
    b::Vector{Float64}
)::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    for (i, k) in enumerate(bnodes)
        b[k] = f[i]

        indices1 = findall(A[1:k-1, k] .!= 0) 
        indices2 = findall(A[k, 1+k:end] .!= 0) .+ k

        b[indices1] -= A[indices1,k] .* f[i]
        b[indices2] -= A[k,indices2] .* f[i]

        A[indices1, k] .= 0
        A[k, indices2] .= 0
        A[k, k] = 1
    end

    return A, b
end

function edgeIndices(
    EToV::Matrix{Int64},
    beds::Matrix{Int64},
)::Tuple{Vector{Int64}, Vector{Int64}}
    n = beds[:, 1]
    r = beds[:, 2]
    s = r .% 3 .+ 1

    i = EToV[CartesianIndex.(n, r)]
    j = EToV[CartesianIndex.(n, s)]

    return i,j
end

function neubc(
    VX::Vector{Float64},
    VY::Vector{Float64},
    EToV::Matrix{Int64},
    beds::Matrix{Int64},
    q::Vector{Float64},
    b::Vector{Float64}
)::Vector{Float64}
    i, j = edgeIndices(EToV, beds)
    q1 = q .* sqrt.((VX[j] - VX[i]).^2 + (VY[j] - VY[i]).^2) ./ 2

    b[i] -= q1
    b[j] -= q1

    return b
end


function assembly(
    VX::Vector{Float64},
    VY::Vector{Float64},
    EToV::Matrix{Int64},
    lam1::Float64,
    lam2::Float64,
    qt::Vector{Float64}
)::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    N = size(EToV)[1]
    M = length(VX)
    A = spzeros(M, M)
    b = zeros(M)

    as, bs, cs = basfun(VX, VY, EToV)
    deltas = sum(as, dims=2) ./ 2
    qs = abs.(deltas) .* sum(qt[EToV], dims=2) / 9

    r = [1,1,1,2,2,3]
    s = [1,2,3,2,3,3]
    i = EToV[:,r]
    j = EToV[:,s]
    ks = (lam1 .* bs[:,r] .* bs[:,s] .+ lam2 .* cs[:,r] .* cs[:,s]) ./ (4 .* abs.(deltas))
    idx = CartesianIndex.(min.(i, j), max.(i, j))

    for n in 1:N
        A[idx[n, :]] += ks[n, :]
        b[EToV[n, :]] .+= qs[n]
    end

    return A, b
end

function assembly2(
    VX::Vector{Float64},
    VY::Vector{Float64},
    EToV::Matrix{Int64},
    lam1::Float64,
    lam2::Float64,
    qt::Vector{Float64}
)::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    N = size(EToV)[1]
    M = length(VX)
    A = spzeros(M, M)
    b = zeros(M)

    as, bs, cs = basfun(VX, VY, EToV)
    deltas = sum(as, dims=2) ./ 2
    qs = abs.(deltas) .* sum(qt[EToV], dims=2) / 9

    for r in 1:3
        i = EToV[:,r]

        for n in 1:N
            b[i[n]] += qs[n]
        end

        for s in r:3
            j = EToV[:,s]
            ks = (lam1 .* bs[:,r] .* bs[:,s] + lam2 .* cs[:,r] .* cs[:,s]) ./ (4 .* abs.(deltas))
            idx = CartesianIndex.(min.(i, j), max.(i, j))

            for n in 1:N
                A[idx[n]] += ks[n]
            end
        end
    end

    return A, b
end

function assembly3(
    VX::Vector{Float64},
    VY::Vector{Float64},
    EToV::Matrix{Int64},
    lam1::Float64,
    lam2::Float64,
    qt::Vector{Float64}
)::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    N = size(EToV)[1]
    M = length(VX)

    A = spzeros(M, M)
    B = zeros(M)

    as, bs, cs = basfun(VX, VY, EToV)
    deltas = sum(as, dims=2) ./ 2
    qs = abs.(deltas) .* sum(qt[EToV], dims=2) / 9

    for n in 1:N
        delta = deltas[n]
        q = qs[n]
        b = bs[n, :]
        c = cs[n, :]

        for r in 1:3
            i = EToV[n,r]
            B[i] += q

            for s in r:3
                j = EToV[n,s]
                kn = (lam1*b[r]*b[s] + lam2*c[r]*c[s]) / (4 * abs(delta))
                A[min(i, j), max(i, j)] += kn
            end
        end
    end

    return A, B
end

function solveNDBVP(
    VX::Vector{Float64},
    VY::Vector{Float64},
    EToV::Matrix{Int64},
    lam1::Float64,
    lam2::Float64,
    qt::Function,
    q::Function,
    f::Function,
    fd_gamma1::Function,
    fd_gamma2::Function,
    tol::Float64
)::Vector{Float64}
    A, b = assembly3(VX, VY, EToV, lam1, lam2, qt.(VX, VY))

    beds = constructBeds(VX, VY, EToV, tol, fd_gamma1)
    i, j = edgeIndices(EToV, beds)
    b = neubc(VX, VY, EToV, beds, q.(VX[i], VY[i], VX[j], VY[j]), b)

    bnodes = constructBnodes(VX, VY, tol, fd_gamma2)
    A, b = dirbc(bnodes, f.(VX[bnodes], VY[bnodes]), A, b)

    p = symamd(A)
    ip = similar(p)
    ip[p] = 1:length(p)

    A = Symmetric(A)

    uhat = A[p,p] \ b[p]
    uhat = uhat[ip]

    return uhat
end

function Driver28b(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, f, qt, q)
    fd_gamma1(x, y) = min(x - x0, y - y0)
    fd_gamma2(x, y) = min(x0 + L1 - x, y0 + L2 - y)
    tol = 0.0001

    VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
    EToV = conelmtab(noelms1, noelms2)
    uhat = solveNDBVP(VX, VY, EToV, lam1, lam2, qt, q, f, fd_gamma1, fd_gamma2, tol)
    return VX, VY, EToV, uhat
end