include("directions.jl")

#=
The following code is based on the Sobol sequence generator by Frances
Y. Kuo and Stephen Joe. The license terms are provided below.

Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

=#

"""
    sobol_sequence(N::Int, D::Int)

Generate `N` by `D` array of Sobol sequence samples
"""
function sobol_sequence(N::Int, D::Int)

    # Generate (N x D) numpy array of Sobol sequence samples
    scale = 31
    result = zeros(N, D)

    if D > length(directions) + 1
        error("Error in Sobol sequence: not enough dimensions")
    end
    
    L = Int(ceil(log(N) / log(2)))
    if L > scale
        error("Error in Sobol sequence: not enough bits")
    end

    for i in 1:D
        V = zeros(UInt64, L+1)

        if i == 1
            for j in 2:L+1
                V[j] = 1 << (scale - (j-1)) # all m's = 1
            end
        else
            m = convert(Array{Int}, directions[i - 1])
            a = m[1]
            s = length(m) - 1

            # The following code discards the first row of the ``m`` array
            # because it has floating point errors, e.g. values of 2.24e-314

            if L <= s
                for j in 2:L+1
                    V[j] = m[j] << (scale - (j-1))
                end
            else
                for j in 2:s+1
                    V[j] = m[j] << (scale - (j-1))
                end
                for j in s+2:L+1
                    V[j] = xor(V[j-s], (V[j - s] >> s))
                    for k in 1:s
                        V[j] = xor(V[j], ((a >> (s - 1 - k)) & 1) * V[j - k])
                    end
                end
            end
        end

        X = UInt64(0)
        for j in 2:N
            X = xor(X, V[index_of_least_significant_zero_bit((j-1) - 1)])
            result[j,i] = Float64(X / (2 ^ scale))
        end
    end

    return result
end


"""
    index_of_least_significant_zero_bit(value::Int)

Return the 1-indexed index of the least significant zero bit of `value`
"""
function index_of_least_significant_zero_bit(value::Int)
    index = 2
    while (value & 1) != 0
        value >>= 1
        index +=1
    end

    return index
end
