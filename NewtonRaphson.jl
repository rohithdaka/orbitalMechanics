#inputs E, e , M, tolerance allowed and iteration check to keep track of 
# recursion
 
function NewtonRaphson(E,e,M,tolerance,iteration)
    F = E - e*sin(E) - M # defining F
    if abs(F) < tolerance
        return (E, iteration) # If F is within tolerance Efinal. so return
    end
    dF = 1 - e*cos(E) # else find dF/dE 
    Enew = E - F/dF
    (Ea , I) = NewtonRaphson(Enew,e,M,tolerance,iteration+1)
end