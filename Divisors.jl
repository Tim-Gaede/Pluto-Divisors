### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ 2efdbfe0-f079-11ea-07b3-9faee3ed637f
using Pkg

# ╔═╡ f113ed44-f078-11ea-3d4a-0953e793934f
using Formatting

# ╔═╡ 3b96db7e-f079-11ea-3a02-4541d747e74e
Pkg.add("Formatting")

# ╔═╡ a6e9d700-f074-11ea-04aa-5ba864c236ae
n = 72

# ╔═╡ 0adf222a-f074-11ea-116b-670c55753921
function divisors_brute(n)
	
	dvrs = [1]	
	for d in (2 : n-1)
		if n % d == 0
			push!(dvrs, d)
		end
	end
	push!(dvrs, n)
	
	dvrs
end 

# ╔═╡ c1bb534e-f077-11ea-1cdf-7140e0d84c23
divisors_brute(n)

# ╔═╡ a127c5b0-f0ca-11ea-2e61-836220b27a63
function divisors_brutish(n)
	
	root_down = convert(Int64, floor(√n)) # square root of n rounded down
	
	dvrs = [1]
	for d in (2 : root_down)
		if n % d == 0
			push!(dvrs, d) 
		end
	end
	
	
	# Now we will iterate backwards through the divisors
	# to generate new divisors.
	
	# If n is square         start with penultimate          else with last  
	root_down^2 == n    ?    iₕᵢ = length(dvrs) - 1    :    iₕᵢ = length(dvrs)	 
	
	for i in (iₕᵢ : -1 : 2)
		push!(dvrs, n ÷ dvrs[i])
	end	
	push!(dvrs, n)
	
	
	dvrs
end

# ╔═╡ 738ca0e2-f0cc-11ea-1130-15a514db81ac
divisors_brutish(n)

# ╔═╡ 5f038976-f075-11ea-3a68-85e1a963077a
# Implementation of the Sieve of Eratosthenes (era-TOSS-the-knees)
# where we take advantage of the option to iterate a for loop by more than 1.
#
# Multiplying two numbers to
# get an odd would imply that both of those numbers were also odd.
function primesTo(n::Int)
    if n < 2;    return []; end
    primes = Int64[2]
    sizehint!(primes, convert( Int64, floor( n / log(n) ) ))
    oddsAlive = trues((n-1) ÷ 2) # oddsAlive[i] represents 2i + 1

    i_sqrt = (convert( Int64, floor(√n) ) - 1) ÷ 2
    for i in (1 : i_sqrt)
        if oddsAlive[i] # It's prime.  Kill odd multiples of it
            push!(primes, 2i + 1)
            Δi = 2i + 1
            for iₓ = i+Δi : Δi : length(oddsAlive);   oddsAlive[iₓ] = false; end
        end
    end
    for i in (i_sqrt + 1 : length(oddsAlive)) # Surviving odds also prime
        if oddsAlive[i];    push!(primes, 2i + 1); end
    end

    primes
end 

# ╔═╡ 719f54de-f075-11ea-1601-fb01690e54c0
primesTo(60)

# ╔═╡ a3ba5f5e-f075-11ea-0071-cb03b010b38c
 
function primeFactors_and_powers!(n::Int, primes::Array{Int,1})
# Will extend the array of primes if it is possibly inadequate for factorization

	if n < 1    throw("n must be a positive integer.")  end
	if n == 1  return [1]  end



	sqrt_n = convert(Int64, floor(√n))
	if sqrt_n ≥ last(primes)
		primesNew = primesUpTo(2sqrt_n)
		len = length(primes)
		for i in (len + 1 : length(primesNew))
			push!(primes, primesNew[i])
		end
	end

		
	fctrs_and_pwrs = []
	rem = n
	root = convert(Int64, floor(√rem + 0.5))
	i = 1
	while rem !== 1    &&    primes[i] ≤ root
		if rem % primes[i] == 0

			pwr = 0
			while rem % primes[i] == 0
				rem ÷= primes[i]
				pwr += 1
			end
			push!(fctrs_and_pwrs, (primes[i], pwr))
		end
		i += 1
		root = convert(Int64, floor(√rem + 0.5))
	end
	if rem !== 1;    push!(fctrs_and_pwrs, (rem, 1));    end

		
	fctrs_and_pwrs
end
	
 


# ╔═╡ 6b552c38-f076-11ea-2b28-b7bd1f4e3503
 

function divisors(fctrs_and_pwrs)	 

	# Use the prime factors and powers to generate all divisors
	result = Int64[]
	numDvrs = Int64(1)
	for f_and_p in fctrs_and_pwrs
		numDvrs *= (f_and_p[2] + 1)
	end

	codeMax = numDvrs - 1
	for code in (0 : codeMax)
		rem = code
		dvr = Int64(1)
		i = 1

		while rem > 0
			exponent = rem % (fctrs_and_pwrs[i][2] + 1)
			dvr *= fctrs_and_pwrs[i][1]^exponent
			rem ÷= (fctrs_and_pwrs[i][2] + 1)
			i += 1
		end

		push!(result, dvr)
	end



	result
end
    
	 

# ╔═╡ c4e40f54-f075-11ea-316e-2992a7c70f2f
primes = primesTo(999)

# ╔═╡ d46a9b5a-f075-11ea-2c4d-b31c7d5699a7
 

fctrs_and_pwrs = primeFactors_and_powers!(n, primes)
 

# ╔═╡ 30b8fce6-f076-11ea-148d-61b9a29982c0
sort(divisors(fctrs_and_pwrs))

# ╔═╡ 622c75ce-f078-11ea-3d79-7958142ac435
 
 
function superscript(n::Int)
    if n == 0;    return "⁰";    end

    res = ""
    n < 0    ?    rem = -n    :    rem = n

    sups = split("¹²³⁴⁵⁶⁷⁸⁹⁰", "")

    while rem > 0
        d = rem % 10 # digit
        rem ÷= 10

        d == 0    ?    res *= string(sups[10])    :    res *= string(sups[d])
    end
    if n < 0;    res *= "⁻";   end

    reverse(res) # returned
end
 
	 


# ╔═╡ 77afb8ca-f078-11ea-3511-d9b5287c8d38
 
 
function underscoresL(n::Number)
# Places uderscores left of the decimal point    

    n_commas = format(n, commas=true)

        
    n_fmd = ""
    for i in (1 :  length(n_commas))
        n_commas[i] == ','    ?    n_fmd *= '_'    :    n_fmd *= n_commas[i]
    end

    n_fmd
end  
    
    

# ╔═╡ bd19eef0-f078-11ea-398e-7138faf17d6b
function formatted_factors(a)
    if length(a) == 0;    return "";    end

    res = underscoresL(a[1][1])  * superscript(a[1][2])

    for i in (2 : length(a))
        
        item = " × " * underscoresL(a[i][1]) * superscript(a[i][2])
        res *= item
        
    end

    res
end

# ╔═╡ 1025a144-f079-11ea-2e40-c90147f94060
formatted_factors(fctrs_and_pwrs)

# ╔═╡ cb02d6a6-f078-11ea-1448-11eabff24393
function main()
    println("\n"^2)

    primes = primesTo(10_000)

    for n in (2 : 999)
        fctrs = primeFactors_and_powers!(n, primes)
        println(lpad(n, 5), " = ", formatted_factors(fctrs) )
    end
 

    println()
 
    println("\n"^3, "▬"^displaysize(stdout)[2])

end

# ╔═╡ cff3a91a-f078-11ea-0d7f-2ba4acf7868a
main()

# ╔═╡ Cell order:
# ╠═2efdbfe0-f079-11ea-07b3-9faee3ed637f
# ╠═3b96db7e-f079-11ea-3a02-4541d747e74e
# ╠═f113ed44-f078-11ea-3d4a-0953e793934f
# ╠═a6e9d700-f074-11ea-04aa-5ba864c236ae
# ╠═0adf222a-f074-11ea-116b-670c55753921
# ╠═c1bb534e-f077-11ea-1cdf-7140e0d84c23
# ╠═a127c5b0-f0ca-11ea-2e61-836220b27a63
# ╠═738ca0e2-f0cc-11ea-1130-15a514db81ac
# ╠═5f038976-f075-11ea-3a68-85e1a963077a
# ╠═719f54de-f075-11ea-1601-fb01690e54c0
# ╠═a3ba5f5e-f075-11ea-0071-cb03b010b38c
# ╠═6b552c38-f076-11ea-2b28-b7bd1f4e3503
# ╠═c4e40f54-f075-11ea-316e-2992a7c70f2f
# ╠═d46a9b5a-f075-11ea-2c4d-b31c7d5699a7
# ╠═30b8fce6-f076-11ea-148d-61b9a29982c0
# ╠═622c75ce-f078-11ea-3d79-7958142ac435
# ╠═77afb8ca-f078-11ea-3511-d9b5287c8d38
# ╠═bd19eef0-f078-11ea-398e-7138faf17d6b
# ╠═1025a144-f079-11ea-2e40-c90147f94060
# ╠═cb02d6a6-f078-11ea-1448-11eabff24393
# ╠═cff3a91a-f078-11ea-0d7f-2ba4acf7868a
