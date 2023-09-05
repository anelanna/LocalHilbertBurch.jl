module LocalHilbertBurch

abstract type CellType end

using Oscar
using Markdown
# Write your package code here.

@doc Markdown.doc"""
    hilbert_burch_matrix(d::Vector{Int64}, x, y)

Compute the Hilbert-Burch matrix from a given degree vector `d`.

# Examples
```jldoctest
julia> using Oscar

julia> d = [1,1,2]
3-element Vector{Int64}:
 1
 1
 2

julia> R,(x,y) = PolynomialRing(QQ,2)
(Multivariate Polynomial Ring in x1, x2 over Rational Field, Nemo.fmpq_mpoly[x1, x2])

julia> hilbert_burch_matrix(d, x, y)
[ x2     0      0]
[-x1    x2      0]
[  0   -x1   x2^2]
[  0     0    -x1]
```
"""
function hilbert_burch_matrix(d::Vector{Int64}, x, y)
    t = length(d)
    H = zero_matrix(parent(x), t+1, t)
    for i in 1:t
        H[i,i] = y^d[i]
        H[i+1,i] = -x
    end
    return H
end

function U_matrix(m::Vector{Int64})
    m[1] == 0 || prepend!(m,0)
    t = length(m)-1
    U = zero_matrix(ZZ, t+1, t)
    for i in 1:t+1
        for j in 1:t
            U[i,j] = m[j+1]-m[i]+i-j
        end
    end
    return U
end

function m_to_d(m::Vector{Int64})
    m[1] == 0 || prepend!(m,0)
    [m[i+1] - m[i] for i in 1:(length(m)-1)]
end

function order_bounds(m::Vector{Int64})
    m[1] == 0 || prepend!(m,0)
    t = length(m) - 1
    U = U_matrix(m)
    OB = zero_matrix(ZZ, t+1, t)
    for i in 1:t+1, j in 1:t
        s = i<=j
        if i<=j
            s = 1
        end
        OB[i,j] = max(U[i,j]+s,0)
    end
    return OB
end

function degree_bounds(m::Vector{Int64})
    m[1] == 0 || prepend!(m,0)
    d = m_to_d(m)
    t = length(m)-1
    DB = zero_matrix(ZZ, t+1, t)
    for i in 1:t+1, j in 1:t
        DB[i,j] = d[min(i,j)]-1
    end
    return DB
end

function ypolywithcoeffs(k::fmpz, l::fmpz, ccount::Int64, v)
    result = zero(parent(v[1]))
    y = v[2]
    for i in k:l
        ccount += 1
        result += v[ccount+2]*(y^i)
    end
    return result,ccount
end
ypolywithcoeffs(k::Int, l::fmpz, ccount::Int64, v) = ypolywithcoeffs(ZZ(k), l, ccount, v)

function ypolywithcoeffs_matrix(m::Vector{Int64}, v)
    m[1] == 0 || prepend!(m,0)
    t = length(m)-1
    OB = order_bounds(m)
    DB = degree_bounds(m)
    R = parent(v[1])
    result = zero_matrix(R, t+1, t)
    ccount = 0
    for i in 1:t+1, j in 1:t
        result[i,j],ccount = ypolywithcoeffs(OB[i,j], DB[i,j], ccount, v)
    end
    return result, ccount
end

function homypolywithcoeffs_matrix(m::Vector{Int64}, v)
    m[1] == 0 || prepend!(m,0)
    t = length(m)-1
    U = U_matrix(m)
    DB = degree_bounds(m)
    R = parent(v[1])
    result = zero_matrix(R, t+1, t)
    ccount = 0
    for j in 1:t, i in j:t+1 
        if (0 <= U[i,j] && U[i,j] <= DB[i,j])
            result[i,j],ccount = ypolywithcoeffs(U[i,j], U[i,j], ccount, v)
        end
    end
    return result, ccount
end

struct Cell <: CellType
    m::Vector{Int64}
    d::Vector{Int64}
    hilb::Vector{fmpq}
    H::MatElem
    E::Ideal
    U::fmpz_mat
    N::MatElem
    M::MatElem
    dim::Int64
    N_hom::MatElem
    M_hom::MatElem
    dim_hom::Int64
    I
end

function hbm_wrap(d::Vector{Int64}, R)
    v = gens(R)
    x,y = v[1],v[2]
    return hilbert_burch_matrix(d, x, y)
end

function hilbert_function_as_vector(n,Q,E)
    R,_ = quo(Q, E)
    s = 0
    while hilbert_function(R,s) != 0
        s = s+1
    end
    return [hilbert_function(R,i) for i in 0:s-1]
end

function Cell(m::Vector{Int64}, R, Q)
    d = copy(m)
    d = m_to_d(d)
    HQ = hbm_wrap(d, Q)
    HR = hbm_wrap(d, R)
    E = ideal(minors(HQ,ncols(HQ)))
    U = U_matrix(copy(m))
    N, dim = ypolywithcoeffs_matrix(copy(m), gens(R))
    M = HR+N
    I = ideal(minors(M,ncols(M)))
    N_hom, dim_hom = homypolywithcoeffs_matrix(copy(m), gens(R))
    M_hom = HR+N_hom
    hilb = hilbert_function_as_vector(sum(m), Q, E)
    return Cell(m, d, hilb, HQ, E, U, N, M, dim, N_hom, M_hom, dim_hom, I)
end


function dimm(m::Vector{Int64})
    m[1] == 0 || prepend!(m,0)
    t = length(m)-1
    d = copy(m)
    d = m_to_d(d)
    OB = order_bounds(m)
    N = 0
    for i in 1:(t+1), j in 1:t
      if OB[i,j] <= d[min(i,j)]-1
         N = N + d[min(i,j)] - OB[i,j]
      end
    end
    return N
end

function dim_homm(m::Vector{Int64})
    m[1] == 0 || prepend!(m,0)
    t = length(m)-1
    d = copy(m)
    d = m_to_d(d)
    U = U_matrix(m)
    DB = degree_bounds(m)
    N = 0
    for j in 1:t, i in j+1:t+1
      if (0 <= U[i,j] && U[i,j] <= DB[i,j])
         N = N + 1
      end
    end
    return N
end

function dim(c::CellType)
   return c.dim
end


function dim_hom(c::CellType)
   return c.dim_hom
end

function sorted_celllist(ResultType::Type{T}, n::Int64) where T
    result = Dict{Int64, Vector{T}}()
    for i in 0:n-1
        result[i] = Vector{T}()
    end
    R,_ = PolynomialRing(QQ, [["x","y"];["c_"*string(i) for i in 1:n]])
    Q,_ = GradedPolynomialRing(QQ, ["x", "y"])
    for partition in Generic.partitions(n)
        c = T(reverse(partition), R, Q)
        push!(result[c.dim],c)
    end
    return result
end
sorted_celllist(n::Int) = sorted_celllist(Cell, n)


struct SmallCell
    m::Vector{Int64}
    hilb::Vector{fmpq}
    dim::Int64
    dim_hom::Int64
end


function SmallCell(m::Vector{Int64}, R, Q)
    d = copy(m)
    d = m_to_d(d)
    HQ = hbm_wrap(d, Q)
    E = ideal(minors(HQ,ncols(HQ)))
    U = U_matrix(copy(m))
    hilb = hilbert_function_as_vector(sum(m), Q, E)
    dim = dimm(m)
    dim_hom = dim_homm(m)
    return SmallCell(m, hilb, dim, dim_hom)
end

function sorted_celllist_by_HF(ResultType::Type{T}, n::Int64) where T
    result = Dict{Vector{fmpq}, Vector{T}}()
    R,_ = PolynomialRing(QQ, [["x","y"];["c_"*string(i) for i in 1:n]])
    Q,_ = GradedPolynomialRing(QQ, ["x", "y"])
    for partition in Generic.partitions(n)
        c = T(reverse(partition), R, Q)
        hf = c.hilb
        if !haskey(result, hf)
            result[hf] = Vector{T}()
        end
        push!(result[hf],c)
    end
    return result
end
sorted_celllist_by_HF(n::Int) = sorted_celllist_by_HF(Cell, n)

function check_dim_diffs(n::Int64)
   l = sorted_celllist_by_HF(SmallCell, n)
   for hf in keys(l)
      selected = l[hf]
      t = selected[1].dim - selected[1].dim_hom
      for sc in selected
         tt = sc.dim - sc.dim_hom
         if t != tt
            println(hf)
            println(sc.m)
            println(sc.dim)
            println(sc.dim_hom)
            println(t)
            return false
         end
      end
   end
   return true
end

export hilbert_burch_matrix, U_matrix, order_bounds, m_to_d, degree_bounds, hilbert_function_ad_vector, Cell, sorted_celllist, sorted_celllist_by_HF, dim, dim_hom, dimm, dim_homm, SmallCell, check_dim_diffs
include("GradedHilbertBurch.jl")

end
