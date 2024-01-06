function errorestimate(
    xc::Vector{Float64},
    xf::Vector{Float64},
    uc::Vector{Float64},
    uf::Vector{Float64},
    idxMarked::Vector{Int64}
)::Vector{Float64}
    K = length(idxMarked)
    err = zeros(K)

    for k in 1:K
        i = idxMarked[k]
        j = i + k - 1
        h = xc[i+1] - xc[i]

        err[i] = sqrt(3)*sqrt(2)*sqrt(h*(2*uf[j+1]^2 + (uf[j] - 3*uc[i] + uf[j+2] - 3*uc[i+1])*uf[j+1] + 2*uc[i+1]^2 + (-uf[j]/2 + 2*uc[i] - (5*uf[j+2])/2)*uc[i+1] + uf[j]^2 - (5*uc[i]*uf[j])/2 + 2*uc[i]^2 - uc[i]*uf[j+2]/2 + uf[j+2]^2))/6
    end

    return err
end

function refine_marked(xc::Vector{Float64}, idxMarked::Vector{Int64})::Vector{Float64}
    M = length(xc)
    K = length(idxMarked)

    xf = zeros(M + K)

    i = 1
    for j in 1:M-1
        xf[j+i-1] = xc[j]

        if i <= K && idxMarked[i] == j
            xf[j+i] = (xc[j+1] + xc[j])/2
            i += 1
        end
    end

    xf[end] = xc[end]

    return xf
end

function BVP1Drhs(c::Float64, d::Float64, x::Vector{Float64}, f::Function) :: Vector{Float64}
    M = length(x)
    
    # Algorithm 1
    A = spzeros(M, M)
    b = zeros(M)

    for i in 1:M-1
        h = x[i+1] - x[i]
        k1 = 1/h + h/3
        k2 = -1/h + h/6

        A[i, i] += k1
        A[i, i+1] = k2
        A[i+1, i] = k2
        A[i+1, i+1] = k1
    end

    for i in 2:M-1
        h1 = x[i] - x[i-1]
        h2 = x[i+1] - x[i]
        b[i] = -f(x[i-1])*h1/6 - f(x[i])*(h1 + h2)/3 - f(x[i+1])*h2/6
    end
    
    # Algorithm 2
    b[1] = c
    b[2] -= A[1,2] * c
    b[M-1] -= A[M-1,M] * d
    b[M] = d
    
    A[1,1] = 1
    A[1,2] = 0
    A[2,1] = 0
    A[M, M] = 1
    A[M-1, M] = 0
    A[M, M-1] = 0

    u = A \ b

    return u
end

function DriverAMR17(
    L::Int64, 
    c::Float64, 
    d::Float64, 
    xc::Vector{Float64}, 
    func::Function, 
    tol::Float64, 
    maxit::Int64
)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    iter = 0
    for _ in 1:maxit
        idxMarked = collect(1:length(xc)-1)
        xf = refine_marked(xc, idxMarked)

        uc = BVP1Drhs(c, d, xc, func)
        uf = BVP1Drhs(c, d, xf, func)

        err = errorestimate(xc, xf, uc, uf, idxMarked)

        idxMarked = findall(err .>= tol)
        if isempty(idxMarked)
            break
        end

        xc = refine_marked(xc, idxMarked)
        iter += 1
    end

    uc = BVP1Drhs(c, d, xc, func)
    
    return xc, uc, iter
end