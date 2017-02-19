function C(z)
# takes z as inputs and returns C functions 
    if z>0 
        return (1 - cos(sqrt(z)))/z;
    elseif z <0
        return (cosh(sqrt(-z)) -1 )/(-z);
    else
        return 1/6;
    end
end

function S(z)
# takes z as inputs and returns S functions 
    if z>0 
        return (sqrt(z) - sin(sqrt(z)))/z^(1.5) ; 
    elseif z <0
        return (sinh(sqrt(-z)) - sqrt(-z))/(-z)^(1.5) ;
    else
        return 1/6;
    end
end

# f, g, fdot and gdot functions 
function f(Rₒ,Vₒ,X,ΔT,α,μ)
    rₒ = norm(Rₒ)
    z = α*X^2
    return (1 - X^2/rₒ*C(z)) 
end

function g(Rₒ,Vₒ,X,ΔT,α,μ)
    z = α*X^2
    return (ΔT- X^3/sqrt(μ)*S(z))
end

function fᵈ(R,Rₒ,Vₒ,X,ΔT,α,μ)
    rₒ = norm(Rₒ) # Distance of satellite from central body
    r = norm(R)
    z = α*X^2
    return sqrt(μ)/r/rₒ*(z*S(z) - 1)*X
end

function gᵈ(R,Rₒ,Vₒ,X,ΔT,α,μ)
    r = norm(R)
    z = α*X^2    
    return (1 - X^2/r*C(z))
end


function NewtonRaphsonUniversal(X,M,μ,rₒ,vₒʳ,α,tolerance,iteration)    
    σₒ = rₒ*vₒʳ/sqrt(μ)
    z = α*X^2 
    F = σₒ*X^2*C(z) + (1-α*rₒ)*X^3*S(z) + rₒ*X - M/abs(α)^1.5;  # M = n*dt = sqrt(u/a^3)*dt 
     
    if (abs(F) < tolerance) || (iteration > 30) 
        return X# If F is within tolerance Xfinal is reached. so return
    end
         
    dF = σₒ*X*(1-α*X^2*S(z))+ (1-α*rₒ)*X^2*C(z)+rₒ; # else find dF/dX 
     
    Xnew = X - F/dF;
    return NewtonRaphsonUniversal(Xnew,M,μ,rₒ,vₒʳ,α,tolerance,iteration+1); # loop again  
end

function OrbitalStateAfterTime(Rₒ,Vₒ,μ,ΔT)
    # Rₒ and Vₒ are Observed Position and Velocity Vectors orbiting celestial body with
    # μ as standard gravitational parameter. 
    # This function calculates Position and Velocity of satellite 
    # after ΔT seconds
    rₒ = norm(Rₒ) # Distance of satellite from central body
    vₒ = norm(Vₒ) # Speed of satellite with respect to central body.

    #From energy equation, inverse of semi major axis can be calculated
    α = -(vₒ^2/μ  - 2/rₒ) 
    a = 1/α # semi major axis of the conic section that belongs to the satellite under consideration
    # Mean Anomaly after ΔT seconds 
    M = sqrt(μ*abs(α)^3)*ΔT

    vₒʳ= sum(Rₒ.*Vₒ)/rₒ # Radial velocity of satellite 

    Xᵢ = M/sqrt(abs(α)) # Initial Guess for Universal NewtonRaphson method 
    tolerance = 1e-9 # Acceptable tolerance for the final answer.
    X = NewtonRaphsonUniversal(Xᵢ,M,μ,rₒ,vₒʳ,α,tolerance,1)

    R = f(Rₒ,Vₒ,X,ΔT,α,μ)*Rₒ + g(Rₒ,Vₒ,X,ΔT,α,μ)*Vₒ
    V = fᵈ(R,Rₒ,Vₒ,X,ΔT,α,μ)*Rₒ + gᵈ(R,Rₒ,Vₒ,X,ΔT,α,μ)*Vₒ

    return [R V]
end


