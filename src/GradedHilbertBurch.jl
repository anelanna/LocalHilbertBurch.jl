function W_matrix(m::Vector{Int64},a::Int,b::Int)
    m[1] == 0 || prepend!(m,0)
    t = length(m)-1
    W = zero_matrix(QQ, t+1, t)
    for i in 1:t+1
        for j in 1:t
            W[i,j] = m[j+1]-m[i]+(a//b)*(i-j)
        end
    end
    return W
end

function W_degree_bounds(m::Vector{Int64},a::Int, b::Int)
    m[1] == 0 || prepend!(m,0)
    t = length(m) - 1
    W = W_matrix(m,a,b)
    WDB = zero_matrix(ZZ,t+1, t)
    for i in 1:t+1, j in 1:t
        if i<=j
            WDB[i,j] = numerator(ceil(W[i,j])-1)
        elseif i>j
            WDB[i,j] = numerator(floor(W[i,j]))
        end
    end
    return WDB
end

function graded_ypolywithcoeffs_matrix(m::Vector{Int64}, a::Int, b::Int, v)
    m[1] == 0 || prepend!(m,0)
    t = length(m)-1
    DB = degree_bounds(m)
    WDB = W_degree_bounds(m,a,b)
    R = parent(v[1])
    result = zero_matrix(R, t+1, t)
    ccount = 0
    for i in 1:t+1, j in 1:t
        result[i,j],ccount = ypolywithcoeffs(0, min(DB[i,j],WDB[i,j]), ccount, v)
    end
    return result, ccount
end

struct GradedCell <: CellType
    a::Int
    b::Int
    m::Vector{Int64}
    d::Vector{Int64}
    hilb::Vector{fmpq}
    H::MatElem
    E::Ideal
    W::fmpq_mat
    N::MatElem
    M::MatElem
    dim::Int64
    I
end

function GradedCell(m::Vector{Int64}, a, b, R, Q)
    d = copy(m)
    d = m_to_d(d)
    HQ = hbm_wrap(d, Q)
    HR = hbm_wrap(d, R)
    E = ideal(minors(HQ,ncols(HQ)))
    W = W_matrix(copy(m),a,b)
    N, dim = graded_ypolywithcoeffs_matrix(copy(m),a,b, gens(R))
    M = HR+N
    I = ideal(minors(M,ncols(M)))
    hilb = hilbert_function_as_vector(sum(m), Q, E)
    return GradedCell(a, b, m, d, hilb, HQ, E, W, N, M, dim, I)
end

function graded_sorted_celllist(n::Int64,a::Int,b::Int)
    result = Dict{Int64, Vector{CellType}}()
    for i in 0:2n
        result[i] = Vector{CellType}()
    end
    R,_ = PolynomialRing(QQ, [["x","y"];["c_"*string(i) for i in 1:2n]])
    Q,_ = GradedPolynomialRing(QQ, ["x", "y"])
    for partition in Generic.partitions(n)
        c = GradedCell(reverse(partition), a, b, R, Q)
        push!(result[dim(c)],c)
    end
    return result
end


export W_matrix,W_degree_bounds,graded_ypolywithcoeffs_matrix, GradedCell, graded_sorted_celllist
